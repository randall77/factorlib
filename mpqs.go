package factorlib

import (
	"github.com/randall77/factorlib/big"
	"github.com/randall77/factorlib/linear"
	"math/rand"
	"log"
)

func init() {
	factorizers["mpqs"] = mpqs
}

// mpqs = Multiple Polynomial Quadratic Sieve
//
// define f(x) = (ax+b)^2 - n
//
// f(x) = a^2x^2 + 2abx+b^2-n
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

// f(x0+m)/a = (a^2(x0+m)^2 + 2ab(x0+m) + b^2-n)/a
//           = a(x0+m)^2 + 2b(x0+m) + c
//           = ax0^2+2ax0m+am^2+2bx0+2bm+c
//           = ax0^2+2bx0+c + am^2+2ax0m+2bm
//           = 0 + (am+2ax0+2b)m
//           = (am + 2sqrt(n) - 2b + 2b) m
//           = (am + 2sqrt(n)) m
// choose a <= 2sqrt(n)/m, after that f(x0+m)/a starts growing linearly with a.

func mpqs(n big.Int, rnd *rand.Rand) []big.Int {
	// mpqs does not work for powers of a single prime.  Check that first.
	if f := primepower(n, rnd); f != nil {
		return f
	}

	// Pick a factor base
	fb, a := makeFactorBase(n)
	if a != 0 {
		return []big.Int{big.Int64(a), n.Div64(a)}
	}

	maxp := fb[len(fb)-1]

	// Figure out maximum possible a we want.
	amax := n.SqrtCeil().Lsh(2).Div64(sieverange)
	// Point to stop adding more factors to a.
	// It is ok if a is a bit small.
	amin := amax.Div64(maxp)

	// matrix is used to do gaussian elimination on mod 2 exponents.
	m := linear.NewMatrix(uint(len(fb)))
	
	// pair up large primes that we find using this table
	type largerecord struct {
		x big.Int
		f []uint
	}
	largeprimes := map[int64]largerecord{}
	
	for {
		// Pick an a.  Use a random product of factor base primes
		// that multiply to at least amin.
		af := map[uint]uint{}
		a := big.One
		for a.Cmp(amin) < 0 {
			f := uint(1+rand.Intn(len(fb)-1))
			af[f] += 1
			a = a.Mul64(fb[f])
		}

		// Pick b = sqrt(n) mod a
		var pp []primePower
		for i, k := range af {
			pp = append(pp, primePower{fb[i],k})
		}
		b := bigSqrtModN(n.Mod(a), pp, rnd)

		// Set c = (b^2-n)/a
		c := b.Square().Sub(n).Div(a)

		x0 := n.SqrtCeil().Sub(b).Div(a)
		for _, r := range sievesmooth(a, b.Lsh(1), c, fb, x0.Sub64(sieverange/2), rnd) {
			x := r.x
			factors := r.factors
			remainder := r.remainder
			/*
			fmt.Printf("%d*%d^2+%d*%d+%d=%d=", a, x, b, x, c, a.Mul(x).Add(b).Mul(x).Add(c))
			for i, f := range factors {
				if i != 0 {
					fmt.Printf("·")
				}
				fmt.Printf("%d", fb[f])
			}
			if remainder != 1 {
				fmt.Printf("·%d", remainder)
			}
			fmt.Println()
			*/
			for f, k := range af {
				for i := uint(0); i < k; i++ {
					factors = append(factors, f)
				}
			}
			x = a.Mul(x).Add(b)
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
					log.Printf("%d/%d falsepos=%d largeprimes=%d\n", m.Rows(), len(fb), falsepos, len(largeprimes))
					falsepos = 0
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
				log.Println("triv A")
				continue
			}
			if a.Add(b).Cmp(n) == 0 {
				// trivial equation, ignore it
				log.Println("triv B")
				continue
			}

			r := a.Add(b).GCD(n)
			return []big.Int{r, n.Div(r)}
		}
	}
}
