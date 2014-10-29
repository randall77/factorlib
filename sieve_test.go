package factorlib

import (
	"math/big"
	"math/rand"
	"testing"
	"fmt"
)

func TestSieve(t *testing.T) {
	rnd := rand.New(rand.NewSource(123))

	fb := []int64{2,3,5,7,11,13,17,19}

	var a, b, c big.Int
	a.SetInt64(23)
	b.SetInt64(83)
	c.SetInt64(1011)

	for _, r := range sievesmooth(a, b, c, fb, rnd) {
		fmt.Printf("f(%d)=%d = %#v %d\n", &r.x, 23*r.x.Int64()*r.x.Int64()+83*r.x.Int64()+1011,r.factors, r.remainder)
	}
}
