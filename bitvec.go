package factorlib

type bitVec []uint

// wordsize == {32,64} when uint == uint{32,64}
const wordsize = 32 << (^uint(0) >> 63)

func newBitVec(n uint) bitVec {
	return make([]uint, (n+wordsize-1)/wordsize)
}

func (b bitVec) getBit(i uint) bool {
	return b[i/wordsize]>>(i%wordsize) & 1 != 0
}

func (b bitVec) setBit(i uint) {
	b[i/wordsize] |= uint(1) << (i%wordsize)
}

func (b bitVec) toggleBit(i uint) {
	b[i/wordsize] ^= uint(1) << (i%wordsize)
}

func (b bitVec) empty() bool {
	for _, x := range b {
		if x != 0 {
			return false
		}
	}
	return true
}

func (b bitVec) firstBit() uint {
	for i, x := range b {
		if x != 0 {
			for j := uint(0); j < wordsize; j++ {
				if x & 1 != 0 {
					return uint(i)*wordsize + j
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
