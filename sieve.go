package factorlib

import (
	"fmt"
	"github.com/randall77/factorlib/big"
	"github.com/randall77/factorlib/linear"
	"math/rand"
)

// sieve [-sieverange,sieverange) around mininum point.
const sieverange = 1 << 14

// use an array of this size to do the sieving
const window = 1 << 9

// check to see if anything returned from the sieve isn't actually usable
const checkFalsePositive = false

// check to see if anything usable isn't returned by the sieve
const checkFalseNegative = false

// Records f(x) == product(factors)*remainder
// The values in factors are indexes into the factor base
type sieveResult struct {
	x         big.Int
	factors   []uint
	remainder int64
}

// Find values of x for which f(x) = a x^2 + b x + c factors (within one bigprime) over the primes in fb.
// requires: a > 0
func sievesmooth(a, b, c big.Int, fb []int64, rnd *rand.Rand, fn func(big.Int, []uint, int64) bool) {
	maxp := fb[len(fb)-1]

	// find approximate zero crossings
	d := b.Square().Sub(a.Mul(c).Lsh(2))
	if d.Sign() < 0 {
		panic("polynomial has no roots")
		// TODO: choose min instead?  Then x = -b/2a
	}
	x := b.Neg().Add(d.SqrtFloor()).Div(a).Rsh(1)
	//x2 := b.Neg().Sub(d).Div(a).Rsh(1)
	// TODO: sieve around x2 also? (if d != 0)

	// starting point
	x0 := x.Sub64(sieverange)

	// results buffer
	var factors []uint

	// find starting points
	si := makeSieveInfo2(a, b, c, x0, fb, rnd)

	// pick threshold
	threshold := byte(a.Mul(x0).Add(b).Mul(x0).Add(c).BitLen()) - 2*log2(maxp) // TODO: subtract more?

	// sieve to find any potential smooth f(x)
	sieve := make([]byte, window) // TODO: cache this?
	res := sieveinner(sieve, si, threshold)

	s := &big.Scratch{}

	// check potential results using trial factorization
	for _, i := range res {
		// compute y=f(x)
		x := x0.Add64(int64(i))
		y := a.Mul(x).Add(b).Mul(x).Add(c)

		// trial divide y by the factor base
		// accumulate factor base indexes of factors
		factors = factors[:0]
		for k, p := range fb {
			if p == -1 {
				if y.Sign() < 0 {
					y = y.Neg()
					factors = append(factors, uint(k))
				}
				continue
			}
			for y.Mod64s(p, s) == 0 {
				y = y.Div64(p)
				factors = append(factors, uint(k))
			}
		}

		// if remainder > B^2, it's too big, might not be prime.
		if y.Cmp64(maxp*maxp) > 0 {
			if checkFalsePositive {
				var f []int64
				for _, j := range factors {
					f = append(f, fb[j])
				}
				fmt.Printf("  false positive x=%d f(x)=%d f=%v·%d\n", x, a.Mul(x).Add(b).Mul(x).Add(c), f, y)
			}
			continue
		}

		if fn(x, dup(factors), y.Int64()) {
			return
		}
	}
	if checkFalseNegative {
		// This is expensive, we have to trial divide all the numbers in the
		// sieve range.  Oh well, testing.
	checkloop:
		for i := 0; i < 2*sieverange; i++ {
			for _, r := range res {
				if r == i {
					continue checkloop
				}
			}
			x := x0.Add64(int64(i))
			y := a.Mul(x).Add(b).Mul(x).Add(c)
			var f []int64
			for _, p := range fb {
				if p == -1 {
					if y.Sign() < 0 {
						y = y.Neg()
						f = append(f, p)
					}
					continue
				}
				for y.Mod64s(p, s) == 0 {
					y = y.Div64(p)
					f = append(f, p)
				}
			}
			
			// if remainder <= B^2, leftover is a prime
			if y.Cmp64(maxp*maxp) <= 0 {
				fmt.Printf("  false negative x=%d f(x)=%d f=%v·%d\n", x, a.Mul(x).Add(b).Mul(x).Add(c), f, y)
			}
		}
	}
}

// TODO: write this in assembly?
func sieveinner(sieve []byte, si []sieveinfo2, threshold byte) []int {
	var r []int
	for i := 0; i < 2*sieverange; i += window {
		// clear sieve
		for j := 0; j < window; j++ {
			sieve[j] = 0
		}
		// increment sieve entries for f(x) that are divisible
		// by each factor base prime.
		for j := range si {
			f := &si[j]
			pk := int(f.pk)
			lg_p := f.lg_p
			j := int(f.off)
			for ; j < window; j += pk {
				sieve[j] += lg_p
			}
			f.off = int32(j - window) // for next time
		}
		for j := 0; j < window; j++ {
			if sieve[j] >= threshold {
				r = append(r, i+j)
			}
		}
	}
	return r
}

type sieveinfo2 struct {
	pk   int32 // p^k for this factor base entry
	lg_p uint8 // ~log_2(p)
	off  int32 // working offset in sieve array
}

func makeSieveInfo2(a, b, c big.Int, start big.Int, fb []int64, rnd *rand.Rand) []sieveinfo2 {
	var si []sieveinfo2
	s := &big.Scratch{}

	for _, p := range fb[1:] {
		pk := p
		for k := uint(1); ; k++ {
			if pk > fb[len(fb)-1] {
				// Kind of arbitrary, but use powers of p as long as p^k is
				// smaller than than the maximum factor base prime.
				break
			}
			st := start.Mod64s(pk, s)
			for _, r := range quadraticModPK(a.Mod64s(pk, s), b.Mod64s(pk, s), c.Mod64s(pk, s), p, k, pk, rnd) {
				// find first pk*i+r which is >= start
				off := (r - st + pk) % pk
				si = append(si, sieveinfo2{int32(pk), log2(p), int32(off)})
			}
			pk *= p
		}
	}
	return si
}

func init() {
	factorizers["qs2"] = qs2
}

func qs2(n big.Int, rnd *rand.Rand) []big.Int {
	// qs does not work for powers of a single prime.  Check that first.
	if f := primepower(n, rnd); f != nil {
		return f
	}

	// first, pick a factor base
	fb, a := makeFactorBase(n)
	if a != 0 {
		return []big.Int{big.Int64(a), n.Div64(a)}
	}

	// matrix is used to do gaussian elimination on mod 2 exponents.
	m := linear.NewMatrix(uint(len(fb)))

	// pair up large primes that we find
	type largerecord struct {
		x big.Int
		f []uint
	}
	largeprimes := map[int64]largerecord{}

	var result []big.Int

	// function to process sieve results
	fn := func(x big.Int, factors []uint, remainder int64) bool {
		fmt.Printf("%d^2-%d=%d=", x, n, x.Mul(x).Sub(n))
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

		if remainder != 1 {
			// try to find another record with the same largeprime
			lr, ok := largeprimes[remainder]
			if !ok {
				// haven't seen this large prime yet.  Save record for later
				largeprimes[remainder] = largerecord{x, factors}
				return false
			}
			// combine current equation with other largeprime equation
			// x1^2 === prod(f1) * largeprime
			// x2^2 === prod(f2) * largeprime
			fmt.Printf("  largeprime %d match\n", remainder)
			x = x.Mul(lr.x).Mod(n).Mul(big.Int64(remainder).ModInv(n)).Mod(n) // TODO: could remainder divide n?
			factors = append(factors, lr.f...)
		}

		idlist := m.AddRow(factors, eqn{x, factors})
		if idlist == nil {
			fmt.Println(m.Rows())
			return false
		}

		// we found a set of equations with all even powers
		// compute a and b where a^2 === b^2 mod n
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
			fmt.Println("triv A")
			return false
		}
		if a.Add(b).Cmp(n) == 0 {
			// trivial equation, ignore it
			fmt.Println("triv B")
			return false
		}

		r := a.Add(b).GCD(n)
		result = []big.Int{r, n.Div(r)}
		return true
	}

	sievesmooth(big.Int64(1), big.Int64(0), n.Neg(), fb, rnd, fn)

	// TODO: after sieving done, restart with larger interval?

	return result
}
