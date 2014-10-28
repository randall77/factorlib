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
	v.toggleBit(800)
	v.toggleBit(900)
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