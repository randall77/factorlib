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

func (b bitVec) toggleBit(i int) {
	b[i>>6] ^= uint64(1) << uint(i&63)
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

// b ^= c.  b and c must be constructed using the same length.
func (b bitVec) xor(c bitVec) {
	for i, x := range c {
		b[i] ^= x
	}
}
