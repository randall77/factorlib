package factorlib

import (
	"fmt"
	"github.com/randall77/factorlib/big"
	"github.com/randall77/factorlib/primes"
	"math/rand"
)

func init() {
	factorizers["qs"] = qs
}

func mpqs(n big.Int, rnd *rand.Rand) []big.Int {
	// first, pick a factor base
	fb, a := makeFactorBase(n)
	if a != 0 {
		return []big.Int{big.Int64(a), n.Div64(a)}
	}

	maxp := fb[len(fb)-1]
	fmt.Printf("maxp=%d len(fb)=%d\n", maxp, len(fb))

	for i := 0; ; i++ {
		// Pick a
		a := primes.Get(i)

		// Pick b = sqrt(n) mod a
		b := sqrtModP(n.Mod64(a), a, rnd)
		_ = b

		
	}
}
