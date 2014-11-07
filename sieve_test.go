package factorlib

import (
	"math/rand"
	"testing"
	"fmt"
)

func TestSieve(t *testing.T) {
	rnd := rand.New(rand.NewSource(123))

	fb := []int64{2,3,5,7,11,13,17,19}

	for _, r := range sievesmooth(NewBig(23), NewBig(83), NewBig(1011), fb, rnd) {
		fmt.Printf("f(%d)=%d = %#v %d\n", &r.x, 23*r.x.Int64()*r.x.Int64()+83*r.x.Int64()+1011,r.factors, r.remainder)
	}
}
