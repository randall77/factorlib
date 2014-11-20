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
	amax := n.SqrtCeil().Lsh(1).Div64(sieverange)
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
		// that multiply to at least amax.
		af := map[int64]uint{}
		a := big.One
		for a.Cmp(amin) < 0 {
			f := int64(1+rand.Intn(len(fb)-1))
			af[f] += 1
			a = a.Mul64(fb[f])
		}
		log.Printf("a=%d af=%v amin=%d amax=%d\n", a, af, amin, amax)

		// Pick b = sqrt(n) mod a
		var pp []primePower
		for i, k := range af {
			pp = append(pp, primePower{fb[i],k})
		}
		b := bigSqrtModN(n.Mod(a), pp, rnd)

		// Set c = (b^2-n)/a
		c := b.Square().Sub(n).Div(a)

		// function to process sieve results
		fn := func(x big.Int, factors []uint, remainder int64) []big.Int {
			_ = m
			_ = largeprimes
			return nil
		}

		x0 := n.SqrtCeil().Sub(b).Div(a)
		r := sievesmooth2(a, b.Lsh(1), c, fb, rnd, x0, fn)
		if r != nil {
			return r
		}
	}
	return nil
}

type byInt64Value []int64

func (a byInt64Value) Len() int           { return len(a) }
func (a byInt64Value) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a byInt64Value) Less(i, j int) bool { return a[i] < a[j] }
