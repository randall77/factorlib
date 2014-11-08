package factorlib

import (
	"math/rand"
	"github.com/randall77/factorlib/big"
	"github.com/randall77/factorlib/primes"
)

func init() {
	// Note: this isn't a complete algorithm, so it will
	// cause the driver to fail if n is not a prime power.
	factorizers["primepower"] = primepower
}

// If n is a prime power, factor n.  Otherwise, return nil.
func primepower(n big.Int, rnd *rand.Rand) []big.Int {
	// there is probably a faster way, but this is fast enough.
	for i := 0; ; i++ {
		p := primes.Get(i)
		x := root(n, p)
		if x.Cmp(big.One) <= 0 {
			return nil
		}
		if x.Exp(p).Cmp(n) == 0 {
			var a []big.Int
			for j := int64(0); j < p; j++ {
				a = append(a, x)
			}
			return a
		}
	}
}

// computes floor(n^(1/k))
func root(n big.Int, k int64) big.Int {
	lo := big.One.Lsh(uint((int64(n.BitLen())-1)/k))
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
