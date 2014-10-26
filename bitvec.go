package factorlib

type bitVec []uint64

func newBitVec(n int) bitVec {
	return make([]uint64, (n+63)>>6)
}

func (b bitVec) getBit(i int) bool {
	return b[i>>6]>>uint(i&63) & 1 != 0
}

func (b bitVec) setBit(i int) {
	b[i>>6] |= uint64(1) << uint(i&63)
}

// inverts the setting of bit i.  Returns the old setting.
func (b bitVec) toggleBit(i int) bool {
	r := b[i>>6] >> uint(i&63) & 1 != 0
	b[i>>6] ^= uint64(1) << uint(i&63)
	return r
}

func (b bitVec) empty() bool {
	for _, x := range b {
		if x != 0 {
			return false
		}
	}
	return true
}

func (b bitVec) firstBit() int {
	for i, x := range b {
		if x != 0 {
			for j := 0; j < 64; j++ {
				if x & 1 != 0 {
					return i<<6 + j
				}
				x >>= 1
			}
		}
	}
	panic("firstBit not defined on empty bitVec")
}

// b ^= c
func (b bitVec) xor(c bitVec) {
	for i, x := range c {
		b[i] ^= x
	}
}

func (b bitVec) erase() {
	for i := range b {
		b[i] = 0
	}
}

func (b bitVec) copy() bitVec {
	a := make(bitVec, len(b))
	copy(a, b)
	return a
}
