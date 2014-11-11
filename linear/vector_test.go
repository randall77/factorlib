package linear

import (
	"testing"
)

func TestVector(t *testing.T) {
	v := newVector(1000)
	if !v.empty() {
		t.Errorf("empty failed")
	}
	v.setBit(900)
	if v.getBit(100) {
		t.Errorf("bit 100 shouldn't be set")
	}
	if !v.getBit(900) {
		t.Errorf("bit 900 should be set")
	}
	if v.firstBit() != 900 {
		t.Errorf("first bit didn't return 900")
	}
	if v.empty() {
		t.Errorf("empty failed")
	}
	v.setBit(700)
	v.setBit(800)
	if v.firstBit() != 700 {
		t.Errorf("bad first bit")
	}
	v.toggleBit(700)
	if v.getBit(700) {
		t.Errorf("bad toggle")
	}
}
func TestEmptyFirstBit(t *testing.T) {
	v := newVector(500)
	defer func() {
		a := recover()
		if a == nil {
			panic("firstBit of empty vector did not panic")
		}
	}()
	_ = v.firstBit()
}

func TestXor(t *testing.T) {
	for _, size := range []uint{0, 1, 10, 63, 64, 65, 999} {
		x := newVector(size)
		y := newVector(size)

		for i := uint(0); i < size; i += 3 {
			x.setBit(i)
		}
		for i := uint(0); i < size; i += 5 {
			y.setBit(i)
		}
		x.xor(y)
		for i := uint(0); i < size; i++ {
			want := false
			if i%3 == 0 || i%5 == 0 {
				want = true
				if i%15 == 0 {
					want = false
				}
			}
			if x.getBit(i) != want {
				t.Errorf("bad xor size=%d bit=%d, want %t", size, i, want)
			}
		}
	}
}
