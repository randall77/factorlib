package sieve

import (
	"math/rand"
	"testing"

	"github.com/randall77/factorlib/big"
	"github.com/randall77/factorlib/math"
	"github.com/randall77/factorlib/primes"
)

func TestSieve(t *testing.T) {
	a := big.Int64(23)
	b := big.Int64(-9813)
	c := big.Int64(1011)
	x0 := big.Int64(426) // ~root of ax^2+bx+c
	x1 := x0.Add64(1 << 24)

	rnd := rand.New(rand.NewSource(123))

	// Pick a small factor base.  Include primes only if p divides ax^2+bx+c for some x.
	// Equivalently, only if ax^2+bx+c===0 mod p has a solution.
	// That happens only if b^2-4ac is a quadratic residue.
	d := b.Mul(b).Sub(a.Mul(c).Lsh(2))
	fb := []int64{-1}
	for i := 0; i < 100; i++ {
		p := primes.Get(i)
		if math.QuadraticResidue(d.Mod64(p), p) {
			fb = append(fb, p)
		}
	}

	results := Smooth(a, b, c, fb, x0, x1, rnd)

	if len(results) == 0 {
		t.Errorf("no results for a=%d b=%d c=%d", a, b, c)
	}

	for _, r := range results {
		y := a.Mul(r.X).Add(b).Mul(r.X).Add(c)
		z := big.One
		for _, f := range r.Factors {
			z = z.Mul64(fb[f])
		}
		if r.Remainder != 1 && !big.Int64(r.Remainder).ProbablyPrime(1000) {
			t.Errorf("bad remainder %d: %v", r.Remainder, r)
		}
		z = z.Mul64(r.Remainder)
		if !z.Equals(y) {
			t.Errorf("bad sieve result z=%d, want %d", z, y)
		}
	}
}
