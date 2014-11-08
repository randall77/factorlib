package big

import (
	"math/big"
	"testing"
)

func TestCmp(t *testing.T) {
	// possible Ints
	b := []string{
		"-4", "-3", "-2", "-1", "0", "1", "2", "3", "4",
		"9223372036854775807",  // 2^63-1
		"9223372036854775808",  // 2^63
		"9223372036854775809",  // 2^63+1
		"-9223372036854775807", // -2^63+1
		"-9223372036854775808", // -2^63
		"-9223372036854775809", // -2^63-1
	}
	// posible int64s
	i := []int64{
		-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5,
		1<<63 - 2,
		1<<63 - 1,
		-1 << 63,
		-1<<63 + 1,
		-1<<63 + 2,
	}
	for _, s := range b {
		x, ok := ParseInt(s)
		if !ok {
			t.Errorf("bad format %s", s)
		}
		for _, y := range i {
			c1 := x.Cmp64(y)
			c2 := x.v.Cmp(big.NewInt(y))
			if c1 != c2 {
				t.Errorf("bad Cmp for %s and %d (Int=%d big.Int=%d)", s, y, c1, c2)
			}
		}
	}
}

func TestSqrtFloor(t *testing.T) {
	for i := int64(0); i < 1000; i++ {
		b := Int64(i)
		c := b.SqrtFloor()
		if c.Mul(c).Cmp(b) > 0 {
			t.Errorf("%d.SqrtFloor() = %d, too large", b, c)
		}
		d := c.Add64(1)
		if d.Mul(d).Cmp(b) <= 0 {
			t.Errorf("%d.SqrtFloor() = %d, too small", b, c)
		}
	}
}

func TestSqrtCeil(t *testing.T) {
	for i := int64(1); i < 1000; i++ {
		b := Int64(i)
		c := b.SqrtCeil()
		if c.Mul(c).Cmp(b) < 0 {
			t.Errorf("sqrtCeil(%d) = %d, too small", b, c)
		}
		d := c.Sub64(1)
		if d.Mul(d).Cmp(b) >= 0 {
			t.Errorf("sqrtCeil(%d) = %d, too large", b, c)
		}
	}
	if !Zero.SqrtCeil().IsZero() {
		t.Errorf("sqrtCeil(0) != 0")
	}
}
