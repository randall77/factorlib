package factorlib

import (
	"math/rand"
)

func init() {
	// Note: this isn't a complete algorithm, so it will
	// cause the driver to fail if n is not a prime power.
	factorizers["primepower"] = primepower
}

// If n is a prime power, factor n.  Otherwise, return nil.
func primepower(n bigint, rnd *rand.Rand) []bigint {
	// there is probably a faster way, but this is fast enough.
	for i := 0; ; i++ {
		p := getPrime(i)
		x := root(n, p)
		if x.Cmp(one) <= 0 {
			return nil
		}
		if x.Exp(p).Cmp(n) == 0 {
			var a []bigint
			for j := int64(0); j < p; j++ {
				a = append(a, x)
			}
			return a
		}
	}
}

// computes floor(n^(1/k))
func root(n bigint, k int64) bigint {
	lo := one.Lsh(uint((int64(n.BitLen())-1)/k))
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
			if lo.Add(one).Exp(k).Cmp(n) <= 0 {
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
