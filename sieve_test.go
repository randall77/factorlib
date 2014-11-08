package factorlib

import (
	"github.com/randall77/factorlib/big"
	"github.com/randall77/factorlib/primes"
	"math/rand"
	"testing"
)

func TestSieve(t *testing.T) {
	a := big.Int64(23)
	b := big.Int64(-9813)
	c := big.Int64(1011)

	rnd := rand.New(rand.NewSource(123))

	// Pick a small factor base.  Include primes only if p divides ax^2+bx+c for some x.
	// Equivalently, only if ax^2+bx+c===0 mod p has a solution.
	// That happens only if b^2-4ac is a quadratic residue.
	d := a.Mul(a).Sub(b.Mul(b).Sub(a.Mul(c).Lsh(2)))
	fb := []int64{-1}
	for i := 0; i < 100; i++ {
		p := primes.Get(i)
		if quadraticResidue(d.Mod64(p), p) {
			fb = append(fb, p)
		}
	}

	cnt := 0
	
	fn := func(x big.Int, factors []uint, remainder int64) bool {
		y := a.Mul(x).Add(b).Mul(x).Add(c)
		z := big.One
		for _, f := range factors {
			z = z.Mul64(fb[f])
		}
		z = z.Mul64(remainder)
		if z.Cmp(y) != 0 {
			t.Errorf("bad sieve result z=%d, want %d", z, y)
		}
		cnt++
		return false
	}
	
	sievesmooth(a, b, c, fb, rnd, fn)

	if cnt == 0 {
		t.Errorf("no results for a=%d b=%d c=%d", a, b, c)
	}
}
