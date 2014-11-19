package factorlib

import (
	"github.com/randall77/factorlib/big"
	"github.com/randall77/factorlib/primes"
	"math/rand"
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
// Where we choose b such that b^2==n mod a so the division works.
//
// Then we sieve ax^2+2bx+c which gives us a congruence of squares.
// Select a = prod primes in factor base

func mpqs(n big.Int, rnd *rand.Rand) []big.Int {
	// mpqs does not work for powers of a single prime.  Check that first.
	if f := primepower(n, rnd); f != nil {
		return f
	}

	// first, pick a factor base
	fb, a := makeFactorBase(n)
	if a != 0 {
		return []big.Int{big.Int64(a), n.Div64(a)}
	}

	maxp := fb[len(fb)-1]

	for {
		// Pick an a
		i := len(fb) / 2
		a := primes.Get(i)

		// Pick b = sqrt(n) mod a
		b := sqrtModP(n.Mod64(a), a, rnd)
		_ = b

	}
	_ = fb
	_ = maxp
	return nil
}
