package factorlib

import (
	"fmt"
	"math/rand"
)

func init() {
	factorizers["qs"] = qs
}

func mpqs(n bigint, rnd *rand.Rand) []bigint {
	// first, pick a factor base
	fb, a := makeFactorBase(n)
	if a != 0 {
		return []bigint{NewBig(a), n.Div64(a)}
	}

	maxp := fb[len(fb)-1]
	fmt.Printf("maxp=%d len(fb)=%d\n", maxp, len(fb))

	for i := 0; ; i++ {
		// Pick a
		a := getPrime(i)

		// Pick b = sqrt(n) mod a
		b := sqrtModP(n.Mod64(a), a, rnd)
		_ = b

		
	}
}
