package factorlib

import (
	"log"
	"math/rand"

	"github.com/randall77/factorlib/big"
	"github.com/randall77/factorlib/primes"
)

func init() {
	factorizers["trial"] = trial
}

// trial tries dividing by 2,3,5,7,11,... until a factor is found.
func trial(n big.Int, rnd *rand.Rand, logger *log.Logger) ([]big.Int, error) {
	s := &big.Scratch{}
	for i := 0; ; i++ {
		p := primes.Get(i)
		if n.Mod64s(p, s) == 0 {
			v := [2]big.Int{big.Int64(p), n.Div64(p)}
			return v[:], nil
		}
	}
}
