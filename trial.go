package factorlib

import (
	"math/big"
	"math/rand"
)

func init() {
	factorizers["trial"] = trial
}

// trial tries dividing by 2,3,5,7,11,... until a factor is found.
func trial(n big.Int, rnd *rand.Rand) []big.Int {
	var p, r big.Int
	for i := 0; ; i++ {
		p.SetInt64(getPrime(i))
		r.Mod(&n, &p)
		if r.Sign() == 0 {
			r.Div(&n, &p)
			v := [2]big.Int{p, r}
			return v[:]
		}
	}
}
