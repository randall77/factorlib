package factorlib

import (
	"fmt"
	"math/rand"
)

// sieve [-sieverange,sieverange) around mininum point.
const sieverange = 1 << 10

// use an array of this size to do the sieving
const window = 1 << 8

// Records f(x) == product(factors)*remainder
// The values in factors are indexes into the factor base
type sieveResult struct {
	x         bigint
	factors   []uint
	remainder int64
}

func sieveinner(sieve []byte, si []sieveinfo2, threshold byte) []int {
	var r []int
	for i := 0; i < 2*sieverange; i += window {
		// clear sieve
		for j := 0; j < window; j++ {
			sieve[j] = 0
		}
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

// Find values of x for which f(x) = a x^2 + b x + c factors (within one bigprime) over fb.
// requires: a > 0
func sievesmooth(a, b, c bigint, fb []int64, rnd *rand.Rand) []sieveResult {
	var result []sieveResult

	maxp := fb[len(fb)-1]
	bigmaxp2 := NewBig(maxp).Square()
	fmt.Printf("maxp=%d maxp2=%d len(fb)=%d\n", maxp, bigmaxp2, len(fb))

	// find approximate zero crossings
	d := b.Square().Sub(a.Mul(c).Lsh(2))
	if d.Sign() < 0 {
		panic("polynomial has no roots")
		// TODO: choose min instead?  Then x = -b/2a
	}
	d = d.SqrtFloor()
	x := b.Neg().Add(d).Div(a).Rsh(1)
	//x2 := b.Neg().Sub(d).Div(a).Rsh(1)
	// TODO: sieve around x2 also? (if d != 0)

	// starting point
	x0 := x.Sub64(sieverange)

	sieve := make([]byte, window) // TODO: cache this?

	var factors []uint

	// find starting points
	si := makeSieveInfo2(a, b, c, x0, fb, rnd)

	// pick threshold
	startf := a.Mul(x0).Add(b).Mul(x0).Add(c)
	fmt.Printf("startf: %d\n", startf)
	threshold := byte(startf.BitLen()) - 2*log2(maxp) // TODO: subtract more?
	fmt.Printf("threshold: %d\n", threshold)
	
	// sieve to find any potential smooth f(x)
	res := sieveinner(sieve, si, threshold)
	
	// check potential results using trial factorization
	for _, i := range res {
		// compute f(x)
		x := x0.Add64(int64(i))
		f := a.Mul(x).Add(b).Mul(x).Add(c)
		
		// trial divide f by the factor base
		// accumulate factor base indexes of factors
		factors = factors[:0]
		for k, p := range fb {
			for f.Mod64(p) == 0 {
				f = f.Div64(p)
				factors = append(factors, uint(k))
			}
		}
		
		// if remainder > B^2, it's too big, might not be prime.
		if f.Cmp(bigmaxp2) > 0 {
			//fmt.Printf("  false positive y=%d z=%d threshold=%d sieve[i]=%d log2(y)=%d log2(y/z)=%d\n", y, bigz, threshold, sieve[i], y.BitLen(), x.Div(y, bigz).BitLen())
			continue
		}
		
		result = append(result, sieveResult{x, dup(factors), f.Int64()})
	}
	return result
}

type sieveinfo2 struct {
	pk   int32 // p^k for this factor base entry
	lg_p uint8 // ~log_2(p)
	off  int32 // working offset in sieve array
}

func makeSieveInfo2(a, b, c bigint, start bigint, fb []int64, rnd *rand.Rand) []sieveinfo2 {
	var si []sieveinfo2

	for _, p := range fb {
		pk := p
		for k := uint(1); k < 2; k++ { // TODO: quadraticModPK
			if pk > fb[len(fb)-1] {
				// Kind of arbitrary, but use powers of p as long as p^k is
				// smaller than than the maximum factor base prime.
				break
			}
			am := a.Mod64(pk)
			bm := b.Mod64(pk)
			cm := c.Mod64(pk)

			for _, r := range quadraticModP(am, bm, cm, pk, rnd) {
				// find first pk*i+r which is >= start
				s := start.Mod64(pk)
				off := (r + pk - s) % pk
				si = append(si, sieveinfo2{int32(pk), log2(p), int32(off)})
				fmt.Printf("%#v\n", si[len(si)-1])
			}
			pk *= p
		}
	}
	return si
}
