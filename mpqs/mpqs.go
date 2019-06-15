package mpqs

import (
	"log"
	"math/rand"
	"runtime"
	"sort"

	"github.com/randall77/factorlib/big"
	"github.com/randall77/factorlib/linear"
	"github.com/randall77/factorlib/math"
	"github.com/randall77/factorlib/primepower"
	"github.com/randall77/factorlib/primes"
	"github.com/randall77/factorlib/sieve"
)

// mpqs = Multiple Polynomial Quadratic Sieve
//
// First, read the description of the regular Quadratic Sieve in ../qs/qs.go.
//
// define f(x) = (ax+b)^2 - n
//
// f(x) = a^2x^2+2abx+b^2-n
//      = a(ax^2+2bx+c)        where c=(b^2-n)/a
//
// We choose a to be a product of primes in the factor base.
// We choose b such that b^2==n mod a so the division for computing c works.
//
// Then we sieve ax^2+2bx+c to find x such that ax^2+2bx+c factors over
// the factor base.
//
// We sieve around x0 = (sqrt(n)-b)/a, because that is where f(x) is small.
//
// How do we choose a?  For larger a we can extract a larger known
// factor from f(x).  But a larger a also mean f(x) grows faster as x
// diverges from x0.  Pick a so that f(x0+m)/a is minimal when sieving
// [-m,m].
//
// f(x0+m)/a = (a^2(x0+m)^2 + 2ab(x0+m) + b^2-n)/a
//           = a(x0+m)^2 + 2b(x0+m) + c
//           = ax0^2+2ax0m+am^2+2bx0+2bm+c
//           = ax0^2+2bx0+c + am^2+2ax0m+2bm
//           = 0 + (am+2ax0+2b)m
//           = (am + 2sqrt(n) - 2b + 2b) m
//           = (am + 2sqrt(n)) m
// choose a <= 2sqrt(n)/m, after that f(x0+m)/a starts growing linearly with a.
//
//

func Factor(n big.Int, rnd *rand.Rand, logger *log.Logger) ([]big.Int, error) {
	// mpqs does not work for powers of a single prime.  Check that first.
	if f, err := primepower.Factor(n, rnd, logger); err == nil {
		return f, nil
	}

	// Pick a factor base
	fb, a := makeFactorBase(n)
	if a != 0 {
		return []big.Int{big.Int64(a), n.Div64(a)}, nil
	}

	maxp := fb[len(fb)-1]

	// Figure out maximum possible a we want.
	amax := n.SqrtCeil().Lsh(1).Div64(sieverange)
	// Point to stop adding more factors to a.
	// It is ok if a is a bit small.
	amin := amax.Div64(maxp)
	logger.Printf("a range: [%d,%d]\n", amin, amax)

	// matrix is used to do gaussian elimination on mod 2 exponents.
	m := linear.NewMatrix(uint(len(fb)))

	// pair up large primes that we find using this table
	type largerecord struct {
		x big.Int
		f []uint
	}
	largeprimes := map[int64]largerecord{}

	// Set up worker goroutines.
	workers := runtime.NumCPU() // TODO: set up as a parameter somehow?
	logger.Printf("workers: %d\n", workers)
	tasks := make(chan task, workers+1)
	results := make(chan result, workers+1)
	for i := 0; i < workers; i++ {
		go worker(tasks, results)
	}

	// Send initial tasks.
	var id int64 // Identify tasks with a sequence number.
	for i := 0; i < workers+1; i++ {
		tasks <- task{
			id:   id,
			n:    n,
			amin: amin,
			fb:   fb,
			seed: rnd.Int63(),
		}
		id++
	}

	// Process results.
	pending := map[int64]result{}
	var next int64
	for {
		// Wait for a result.
		res := <-results
		// Queue up another task.
		tasks <- task{
			id:   id,
			n:    n,
			amin: amin,
			fb:   fb,
			seed: rnd.Int63(),
		}
		id++
		// Record the result.
		pending[res.id] = res

		// Process pending results.
		// We want to process these results in ID order, so that our results
		// are deterministic even when the underlying Go scheduler is not.
		for {
			res, ok := pending[next]
			if !ok {
				break
			}
			delete(pending, next)
			next++

			for _, r := range res.r {
				x := r.X
				factors := r.Factors
				remainder := r.Remainder
				if remainder != 1 {
					// try to find another record with the same largeprime
					lr, ok := largeprimes[remainder]
					if !ok {
						// haven't seen this large prime yet.  Save record for later
						largeprimes[remainder] = largerecord{x, factors}
						continue
					}
					// combine current equation with other largeprime equation
					// x1^2 === prod(f1) * largeprime
					// x2^2 === prod(f2) * largeprime
					//fmt.Printf("  largeprime %d match\n", remainder)
					x = x.Mul(lr.x).Mod(n).Mul(big.Int64(remainder).ModInv(n)).Mod(n) // TODO: could remainder divide n?
					factors = append(factors, lr.f...)
				}

				// Add equation to the matrix
				idlist := m.AddRow(factors, eqn{x, factors})
				if idlist == nil {
					if m.Rows()%100 == 0 {
						logger.Printf("%d/%d falsepos=%d largeprimes=%d\n", m.Rows(), len(fb), sieve.FalsePos(), len(largeprimes))
						sieve.ClearFalsePos()
					}
					continue
				}

				// We found a set of equations with all even powers.
				// Compute a and b where a^2 === b^2 mod n
				a := big.One
				b := big.One
				odd := make([]bool, len(fb))
				for _, id := range idlist {
					e := id.(eqn)
					a = a.Mul(e.x).Mod(n)
					for _, i := range e.f {
						if !odd[i] {
							// first occurrence of this factor
							odd[i] = true
							continue
						}
						// second occurrence of this factor
						b = b.Mul64(fb[i]).Mod(n)
						odd[i] = false
					}
				}
				for _, o := range odd {
					if o {
						panic("gauss elim failed")
					}
				}

				if a.Cmp(b) == 0 {
					// trivial equation, ignore it
					logger.Println("triv A")
					continue
				}
				if a.Add(b).Cmp(n) == 0 {
					// trivial equation, ignore it
					logger.Println("triv B")
					continue
				}
				f := a.Add(b).GCD(n)
				res := []big.Int{f, n.Div(f)}

				// Tell workers to stop.
				close(tasks)

				return res, nil
			}
		}
	}
}

type task struct {
	id   int64
	n    big.Int
	amin big.Int
	fb   []int64
	seed int64
}
type result struct {
	id int64
	r  []sieve.Result
}

func worker(tasks chan task, results chan result) {
	for {
		t, ok := <-tasks
		if !ok {
			// Request channel is closed, we're done.
			return
		}
		id := t.id
		n := t.n
		amin := t.amin
		fb := t.fb
		seed := t.seed

		rnd := rand.New(rand.NewSource(seed))

		// Pick an a.  Use a random product of factor base primes
		// that multiply to at least amin.
		af := map[uint]uint{}
		a := big.One
		for a.Cmp(amin) < 0 {
			f := uint(1 + rnd.Intn(len(fb)-1))
			af[f] += 1
			a = a.Mul64(fb[f])
		}

		// Pick b = sqrt(n) mod a
		var pp []math.PrimePower
		for i, k := range af {
			pp = append(pp, math.PrimePower{fb[i], k})
		}
		sort.Slice(pp, func(i, j int) bool {
			return fb[i] < fb[j]
		})
		b := math.BigSqrtModN(n.Mod(a), pp, rnd)

		// Set c = (b^2-n)/a
		c := b.Square().Sub(n).Div(a)

		// Find best point to sieve around
		x0 := n.SqrtCeil().Sub(b).Div(a)

		var rs []sieve.Result
		for _, r := range sieve.Smooth(a, b.Lsh(1), c, fb, x0.Sub64(sieverange/2), x0.Add64(sieverange/2), rnd) {
			r.X = a.Mul(r.X).Add(b)
			for f, k := range af {
				for i := uint(0); i < k; i++ {
					r.Factors = append(r.Factors, f)
				}
			}
			rs = append(rs, r)
		}
		results <- result{id: id, r: rs}
	}
}

// x^2 === prod(f) mod n
// f is a list of indexes into the factor base
type eqn struct {
	x big.Int
	f []uint
}

// pick some prime factors for our factor base.  If we happen
// upon a factor of n, return it instead.
func makeFactorBase(n big.Int) ([]int64, int64) {
	// upper limit on prime factors (TODO: dependent on n) that we sieve with
	const B = 10000
	fb := []int64{-1}
	s := &big.Scratch{}
	for i := 0; ; i++ {
		p := primes.Get(i)
		if p > B {
			return fb, 0
		}
		a := n.Mod64s(p, s)
		if a == 0 {
			return nil, p
		}
		if math.QuadraticResidue(a, p) {
			// if x^2 == n mod p has no solutions, then
			// x^2-n will never be divisible by p.  So it
			// doesn't need to be in the factor base.
			fb = append(fb, p)
		}
	}
}

// width of window to sieve at once.  TODO: make configurable?
const sieverange = 1 << 28
