package linear

type vector []uint

// wordsize == {32,64} when uint == uint{32,64}
const wordsize = 32 << (^uint(0) >> 63)

func newVector(n uint) vector {
	return make([]uint, (n+wordsize-1)/wordsize)
}

func (b vector) getBit(i uint) bool {
	return b[i/wordsize]>>(i%wordsize) & 1 != 0
}

func (b vector) setBit(i uint) {
	b[i/wordsize] |= uint(1) << (i%wordsize)
}

func (b vector) toggleBit(i uint) {
	b[i/wordsize] ^= uint(1) << (i%wordsize)
}

func (b vector) empty() bool {
	for _, x := range b {
		if x != 0 {
			return false
		}
	}
	return true
}

func (b vector) firstBit() uint {
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
	panic("firstBit not defined on empty vector")
}

// b ^= c.  b and c must be constructed using the same length.
func (b vector) xor(c vector) {
	for i, x := range c {
		b[i] ^= x
	}
}
