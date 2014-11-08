package factorlib

import (
	"math/rand"
	"github.com/randall77/factorlib/big"
)

func init() {
	factorizers["trial"] = trial
}

// trial tries dividing by 2,3,5,7,11,... until a factor is found.
func trial(n big.Int, rnd *rand.Rand) []big.Int {
	for i := 0; ; i++ {
		p := getPrime(i)
		if n.Mod64(p) == 0 {
			v := [2]big.Int{big.Int64(p), n.Div64(p)}
			return v[:]
		}
	}
}
