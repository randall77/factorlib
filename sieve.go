package factorlib

import (
	"fmt"
	"math/rand"

	"github.com/randall77/factorlib/big"
)

// sieve [-sieverange,sieverange) around mininum point.
const sieverange = 1 << 20

// use an array of this size to do the sieving
const window = 1 << 14

// check to see if anything returned from the sieve isn't actually usable
const checkFalsePositive = false

// check to see if anything usable isn't returned by the sieve
const checkFalseNegative = false

// TODO: do this more systematically
var falsepos int

// Records f(x) == product(factors)*remainder
// The values in factors are indexes into the factor base
type sieveResult struct {
	x         big.Int
	factors   []uint
	remainder int64
}

func sievesmooth(a, b, c big.Int, fb []int64, rnd *rand.Rand, fn func(big.Int, []uint, int64) []big.Int) []big.Int {
	// find approximate zero crossings
	d := b.Square().Sub(a.Mul(c).Lsh(2))
	if d.Sign() < 0 {
		panic("polynomial has no roots")
		// TODO: choose min instead?  Then x = -b/2a
	}
	x := b.Neg().Add(d.SqrtFloor()).Div(a).Rsh(1)
	//x2 := b.Neg().Sub(d).Div(a).Rsh(1)
	// TODO: sieve around x2 also? (if d != 0)

	return sievesmooth2(a, b, c, fb, rnd, x.Sub64(sieverange/2), fn)
}

// Find values of x for which f(x) = a x^2 + b x + c factors (within one bigprime) over the primes in fb.
// requires: a > 0
func sievesmooth2(a, b, c big.Int, fb []int64, rnd *rand.Rand, start big.Int, fn func(big.Int, []uint, int64) []big.Int) []big.Int {
	maxp := fb[len(fb)-1]

	// results buffer
	var factors []uint

	// find starting points
	si := makeSieveInfo(a, b, c, start, fb, rnd)

	// pick threshold
	// TODO: sample a few locations in the range.  Just first and last?
	var thresholds [sieverange / window]byte
	for i := 0; i < sieverange; i += window {
		x := start.Add64(int64(i))
		y := a.Mul(x).Add(b).Mul(x).Add(c)
		thresholds[i/window] = byte(y.BitLen()) - 2*log2(maxp) // TODO: subtract more?
	}

	// sieve to find any potential smooth f(x)
	res := sieveinner(si, thresholds[:])

	s := &big.Scratch{}

	// check potential results using trial factorization
	for _, i := range res {
		// compute y=f(x)
		x := start.Add64(int64(i))
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

		if r := fn(x, dup(factors), y.Int64()); r != nil {
			return r
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
			x := start.Add64(int64(i))
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
	return nil
}

// sieveinfo
type sieveinfo3 struct {
	pk   int32  // p^k, the amount we step by through the sieve array
	off  uint16 // offset within the window which holds the next muliple of p^k
	lg_p uint8  // log_2(p), the amount we add to each sieve array slot
}

// a block of sieveinfos
type block struct {
	next *block
	num  int
	data [64 - 2]sieveinfo3
}

func sieveinner(si []sieveinfo2, thresholds []byte) []int {
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
	small := make([]sieveinfo3, 0, nsmall)

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
			small = append(small, sieveinfo3{pk, uint16(off), lg_p})
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
		b.data[b.num] = sieveinfo3{pk, uint16(off % window), lg_p}
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
			small[j] = sieveinfo3{int32(pk), uint16(off - window), lg_p}
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
				c.data[c.num] = sieveinfo3{int32(pk), uint16(off % window), lg_p}
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

func sieveinner_old(si []sieveinfo2, thresholds []byte) []int {
	var r []int
	for i := 0; i < sieverange; i += window {
		// clear sieve
		var sieve [window]byte

		threshold := thresholds[i/window]
		// increment sieve entries for f(x) that are divisible
		// by each factor base prime.
		for j := range si {
			f := &si[j]
			pk := int(f.pk)
			lg_p := f.lg_p
			j := int(f.off)
			for ; j < window; j += pk {
				sieve[j] += lg_p
			}
			f.off = int32(j - window) // for next time
		}
		for j := 0; j < window; j++ {
			if sieve[j] >= threshold {
				r = append(r, i+j)
			}
		}
	}
	return r
}

type sieveinfo2 struct {
	pk   int32 // p^k for this factor base entry
	lg_p uint8 // ~log_2(p)
	off  int32 // working offset in sieve array
}

func makeSieveInfo(a, b, c big.Int, start big.Int, fb []int64, rnd *rand.Rand) []sieveinfo2 {
	var si []sieveinfo2
	s := &big.Scratch{}

	for _, p := range fb[1:] {
		if a.Mod64(p) == 0 {
			// can happen in mpqs if p is one of the factors we used to construct a
			continue
		}
		pk := p
		for k := uint(1); ; k++ {
			if pk > fb[len(fb)-1] {
				// Kind of arbitrary, but use powers of p as long as p^k is
				// smaller than than the maximum factor base prime.
				break
			}
			st := start.Mod64s(pk, s)
			for _, r := range quadraticModPK(a.Mod64s(pk, s), b.Mod64s(pk, s), c.Mod64s(pk, s), p, k, pk, rnd) {
				// find first pk*i+r which is >= start
				off := (r - st + pk) % pk
				si = append(si, sieveinfo2{int32(pk), log2(p), int32(off)})
			}
			pk *= p
		}
	}
	return si
}
