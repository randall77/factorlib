package primepower

import (
	"fmt"
	"log"
	"math/rand"

	"github.com/randall77/factorlib/big"
	"github.com/randall77/factorlib/primes"
)

// If n is a prime power, factor n.  Otherwise, return nil and an error.
func Factor(n big.Int, rnd *rand.Rand, logger *log.Logger) ([]big.Int, error) {
	// there is probably a faster way, but this is fast enough.
	for i := 0; ; i++ {
		p := primes.Get(i)
		x := root(n, p)
		if x.Cmp(big.One) <= 0 {
			return nil, fmt.Errorf("%d not a prime power", n)
		}
		if x.Exp(p).Cmp(n) == 0 {
			var a []big.Int
			for j := int64(0); j < p; j++ {
				a = append(a, x)
			}
			return a, nil
		}
	}
}

// computes floor(n^(1/k))
func root(n big.Int, k int64) big.Int {
	lo := big.One.Lsh(uint((int64(n.BitLen()) - 1) / k))
	hi := lo.Lsh(1)
	if lo.Exp(k).Cmp(n) > 0 {
		panic("low too high")
	}
	if hi.Exp(k).Cmp(n) <= 0 {
		panic("high too low")
	}
	// binary search
	for {
		m := lo.Add(hi).Rsh(1)
		if m.Cmp(lo) == 0 {
			// test result
			if lo.Exp(k).Cmp(n) > 0 {
				panic("root too big")
			}
			if lo.Add(big.One).Exp(k).Cmp(n) <= 0 {
				panic("root too small")
			}
			return lo
		}
		if m.Exp(k).Cmp(n) <= 0 {
			lo = m
		} else {
			hi = m
		}
	}
}
