package factorlib

import (
	"math/big"
	"math/rand"
)

func init() {
	// Note: this isn't a complete algorithm, so it will
	// cause the driver to fail if n is not a prime power.
	factorizers["primepower"] = primepower
}

// If n is a prime power, factor n.  Otherwise, return nil.
func primepower(n *big.Int, rnd *rand.Rand) []*big.Int {
	// there is probably a faster way, but this is fast enough.
	e := big.NewInt(0)
	k := big.NewInt(0)
	for i := 0; ; i++ {
		k.SetInt64(getPrime(i))
		x := root(n, k)
		if x.Cmp(one) <= 0 {
			return nil
		}
		if e.Exp(x, k, nil).Cmp(n) == 0 {
			var a []*big.Int
			for j := int64(0); j < getPrime(i); j++ {
				a = append(a, x)
			}
			return a
		}
	}
}

// computes floor(n^(1/k))
func root(n *big.Int, k *big.Int) *big.Int {
	b := int64(n.BitLen())
	i := k.Int64()
	lo := big.NewInt(1)
	lo.Lsh(lo, uint((b-1)/i))
	hi := big.NewInt(0)
	hi.Lsh(lo, 1)
	if big.NewInt(0).Exp(lo, k, nil).Cmp(n) > 0 {
		panic("low too high")
	}
	if big.NewInt(0).Exp(hi, k, nil).Cmp(n) <= 0 {
		panic("high too low")
	}
	// binary search
	m := big.NewInt(0)
	e := big.NewInt(0)
	for {
		m.Add(lo, hi)
		m.Rsh(m, 1)
		if m.Cmp(lo) == 0 {
			// test result
			if e.Exp(m, k, nil).Cmp(n) > 0 {
				panic("root too big")
			}
			m.Add(m, one)
			if e.Exp(m, k, nil).Cmp(n) <= 0 {
				panic("root too small")
			}
			return lo
		}
		if e.Exp(m, k, nil).Cmp(n) <= 0 {
			lo.Set(m)
		} else {
			hi.Set(m)
		}
	}
}
