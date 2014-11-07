package factorlib

import (
	"math/rand"
)

func init() {
	factorizers["trial"] = trial
}

// trial tries dividing by 2,3,5,7,11,... until a factor is found.
func trial(n bigint, rnd *rand.Rand) []bigint {
	for i := 0; ; i++ {
		p := getPrime(i)
		if n.Mod64(p) == 0 {
			v := [2]bigint{NewBig(p), n.Div64(p)}
			return v[:]
		}
	}
}
