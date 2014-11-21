package factorlib

import (
	"fmt"
	"math"
	"math/rand"

	"github.com/randall77/factorlib/big"
)

// width of window to sieve at once.  TODO: make configurable?
const sieverange = 1 << 24

// use an array of this size to do the sieving
const window = 1 << 14

// check to see if anything returned from the sieve isn't actually usable
const checkFalsePositive = false

// check to see if anything usable isn't returned by the sieve
const checkFalseNegative = false

// TODO: do this more systematically
var falsepos int

// Records f(x) == product(factors)*remainder
// The values in the factors slice are indexes into the factor base
type sieveResult struct {
	x         big.Int
	factors   []uint
	remainder int64
}

// Find values of x for which f(x) = a x^2 + b x + c factors (within one bigprime) over the primes in fb.
// Returns each x found, togegther with the factorization of f(x) into factor base primes and a
// possible bigprime remainder.
// Searches in the window [x0, x0+sieverange).
// requires: a > 0
func sievesmooth(a, b, c big.Int, fb []int64, x0 big.Int, rnd *rand.Rand) []sieveResult {
	var result []sieveResult

	maxp := fb[len(fb)-1]

	// Compute scaling factor for logarithms
	// Note: we assume that the maximum f(x) is achieved at one of the endpoints of the sieve range.
	y0 := a.Mul(x0).Add(b).Mul(x0).Add(c).Abs()
	x1 := x0.Add64(sieverange-1)
	y1 := a.Mul(x1).Add(b).Mul(x1).Add(c).Abs()
	maxf := y0
	if y1.Cmp(maxf) > 0 {
		maxf = y1
	}
	scale := 255/maxf.Log()
	
	// find starting points
	si := makeSieveInfo(a, b, c, fb, x0, scale, rnd)

	// compute thresholds for each window
	var thresholds [sieverange / window]byte
	for i := 0; i < sieverange; i += window {
		x := x0.Add64(int64(i))
		y0 = a.Mul(x).Add(b).Mul(x).Add(c).Abs()
		x = x.Add64(window)
		y1 = a.Mul(x).Add(b).Mul(x).Add(c).Abs()
		m := y0
		if y1.Cmp(m) > 0 {
			m = y1
		}
		thresholds[i/window] = byte(scale * (m.Log() - 2 * math.Log(float64(maxp))))
	}

	// sieve to find any potential smooth f(x)
	res := sieveinner(si, thresholds[:])

	// results buffer
	var factors []uint
	s := &big.Scratch{}

	// check potential results using trial factorization
	for _, i := range res {
		// compute y=f(x)
		x := x0.Add64(int64(i))
		y := a.Mul(x).Add(b).Mul(x).Add(c)

		// trial divide y by the factor base
		// accumulate factor base indexes of factors
		factors = factors[:0]
		for k, p := range fb {
			if p == -1 {
				if y.Sign() < 0 {
					y = y.Neg()
					factors = append(factors, uint(k))
				}
				continue
			}
			for y.Mod64s(p, s) == 0 {
				y = y.Div64(p)
				factors = append(factors, uint(k))
			}
		}

		// if remainder > B^2, it's too big, might not be prime.
		if y.Cmp64(maxp*maxp) > 0 {
			if checkFalsePositive {
				var f []int64
				for _, j := range factors {
					f = append(f, fb[j])
				}
				fmt.Printf("  false positive x=%d f(x)=%d f=%v·%d\n", x, a.Mul(x).Add(b).Mul(x).Add(c), f, y)
			}
			falsepos++
			continue
		}

		result = append(result, sieveResult{x, dup(factors), y.Int64()})
		if len(result) > 2*len(fb) {
			// early out when we're factoring small n
			// TODO: include in spec for sievesmooth function?
			return result
		}
	}
	if checkFalseNegative {
		// This is expensive, we have to trial divide all the numbers in the
		// sieve range.  Oh well, testing.
	checkloop:
		for i := 0; i < sieverange; i++ {
			for _, r := range res {
				if r == i {
					continue checkloop
				}
			}
			x := x0.Add64(int64(i))
			y := a.Mul(x).Add(b).Mul(x).Add(c)
			var f []int64
			for _, p := range fb {
				if p == -1 {
					if y.Sign() < 0 {
						y = y.Neg()
						f = append(f, p)
					}
					continue
				}
				for y.Mod64s(p, s) == 0 {
					y = y.Div64(p)
					f = append(f, p)
				}
			}

			// if remainder <= B^2, leftover is a prime
			if y.Cmp64(maxp*maxp) <= 0 {
				fmt.Printf("  false negative x=%d f(x)=%d f=%v·%d\n", x, a.Mul(x).Add(b).Mul(x).Add(c), f, y)
			}
		}
	}
	return result
}

// sieveinfo
type sieveinfo_inner struct {
	pk   int32  // p^k, the amount we step by through the sieve array
	off  uint16 // offset within the window which holds the next muliple of p^k
	lg_p uint8  // log_2(p), the amount we add to each sieve array slot
}

// a block of sieveinfo_inners
type block struct {
	next *block
	num  int
	data [64 - 2]sieveinfo_inner
}

func sieveinner(si []sieveinfo, thresholds []byte) []int {
	if window >= 1<<16 {
		panic("window too big")
	}

	// compute stats about our sieve work
	maxpk := int32(0)
	nsmall := 0
	for _, s := range si {
		pk := s.pk
		if pk < window {
			nsmall++
		}
		if pk > maxpk {
			maxpk = pk
		}
	}
	nlarge := len(si) - nsmall

	// Slice to store sieve work with pk < window
	small := make([]sieveinfo_inner, 0, nsmall)

	// For pk >= window, arrange sieve work in buckets of size w by offset.
	// The offsets stored in this array are relative to the start
	// of the bucket.  So all entries s in large[3] have an effective
	// offset of 3*window+s.off
	nw := 2 + int(maxpk/window)
	large := make([]*block, nw)

	// Allocate enough free blocks to hold everything.
	// There will be at most one partially empty block for each window.
	// Compute as the maximum # of full blocks + maximum # of nonfull blocks.
	n := nlarge/len(block{}.data) + nw
	blockstore := make([]block, n)
	// start each window with one empty block
	for i := 0; i < nw; i++ {
		large[i] = &blockstore[i]
	}
	// put the rest in a free list
	var freeblocks *block
	if nw < n {
		for i := nw; i < n-1; i++ {
			blockstore[i].next = &blockstore[i+1]
		}
		freeblocks = &blockstore[nw]
	}

	// put each sieve entry into the appropriate bucket
	for _, s := range si {
		pk := s.pk
		off := s.off
		lg_p := s.lg_p
		if pk < window {
			small = append(small, sieveinfo_inner{pk, uint16(off), lg_p})
			continue
		}
		i := off / window
		b := large[i]
		if b.num == len(b.data) {
			c := freeblocks
			freeblocks = c.next
			c.next = b
			large[i] = c
			b = c
		}
		b.data[b.num] = sieveinfo_inner{pk, uint16(off % window), lg_p}
		b.num++
	}
	/*
		log.Printf("len(si)=%d, len(small)=%d", len(si), len(small))
		for i := 0; i < nw; i++ {
			n := 0
			for b := large[i]; b != nil; b = b.next {
				n += b.num
			}
			log.Printf("len(large[%d])=%d", i, n)
		}
	*/

	var r []int
	for i := 0; i < sieverange; i += window {
		// start with a clear sieve
		var sieve [window]byte

		// increment sieve entries for f(x) that are divisible
		// by each factor base prime.

		// first, the small primes
		for j, s := range small {
			pk := int(s.pk)
			off := int(s.off)
			lg_p := s.lg_p
			for ; off < window; off += pk {
				sieve[off] += lg_p
			}
			small[j] = sieveinfo_inner{int32(pk), uint16(off - window), lg_p}
		}

		// second, the large primes
		for b := large[0]; b != nil; {
			for _, s := range b.data[:b.num] {
				pk := int(s.pk)
				off := int(s.off)
				lg_p := s.lg_p
				sieve[off] += lg_p
				off += pk
				j := off / window
				c := large[j]
				if c.num == len(c.data) {
					d := freeblocks
					freeblocks = d.next
					d.next = c
					large[j] = d
					c = d
				}
				c.data[c.num] = sieveinfo_inner{int32(pk), uint16(off % window), lg_p}
				c.num++
			}
			c := b
			b = b.next
			c.num = 0
			c.next = freeblocks
			freeblocks = c
		}

		// shift large buckets down by one
		copy(large, large[1:])
		// allocate a new empty bucket for the last one
		b := freeblocks
		freeblocks = b.next
		b.next = nil
		large[len(large)-1] = b

		// check for smooth numbers
		threshold := thresholds[i/window]
		for j := 0; j < window; j++ {
			if sieve[j] >= threshold {
				r = append(r, i+j)
			}
		}
	}
	return r
}

type sieveinfo struct {
	pk   int32 // p^k for this factor base entry
	lg_p uint8 // ~log_2(p)
	off  int32 // working offset in sieve array
}

func makeSieveInfo(a, b, c big.Int, fb []int64, x0 big.Int, scale float64, rnd *rand.Rand) []sieveinfo {
	var si []sieveinfo
	s := &big.Scratch{}
	maxp := fb[len(fb)-1]

	for _, p := range fb[1:] {
		if a.Mod64(p) == 0 {
			// This can happen in mpqs if p is one of the factors we used to construct a.
			// Ignore this prime - we might miss a few smooth f(x), but such is life.
			continue
		}
		lg_p := byte(scale * math.Log(float64(p))) // note: round down to avoid overflow in sieve buckets
		pk := p
		for k := uint(1); ; k++ {
			if pk > maxp {
				// Kind of arbitrary, but use powers of p as long as p^k is
				// smaller than than the maximum factor base prime.
				break
			}
			st := x0.Mod64s(pk, s)
			for _, r := range quadraticModPK(a.Mod64s(pk, s), b.Mod64s(pk, s), c.Mod64s(pk, s), p, k, pk, rnd) {
				// find first pk*i+r which is >= x0
				off := (r - st + pk) % pk
				si = append(si, sieveinfo{int32(pk), lg_p, int32(off)})
			}
			pk *= p
		}
	}
	return si
}
