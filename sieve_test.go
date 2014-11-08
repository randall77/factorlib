package factorlib

import (
	"github.com/randall77/factorlib/big"
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
		p := getPrime(i)
		if quadraticResidue(d.Mod64(p), p) {
			fb = append(fb, p)
		}
	}

	results := sievesmooth(a, b, c, fb, rnd)
	
	if len(results) == 0 {
		t.Errorf("no results for a=%d b=%d c=%d", a, b, c)
	}
	
	for _, r := range results {
		x := r.x
		y := a.Mul(x).Add(b).Mul(x).Add(c)
		z := big.One
		for _, f := range r.factors {
			z = z.Mul64(fb[f])
		}
		z = z.Mul64(r.remainder)
		if z.Cmp(y) != 0 {
			t.Errorf("bad sieve result %v z=%d, want %d", r, z, y)
		}
	}
}
