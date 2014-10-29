package factorlib

import (
	"fmt"
	"math/big"
	"math/rand"
)

// sieve [-sieverange,sieverange) around mininum point.
const sieverange = 1 << 10

// use an array of this size to do the sieving
const window = 1 << 8

// Records f(x) == product(factors)*remainder
// The values in factors are indexes into the factor base
type sieveResult struct {
	x         big.Int
	factors   []uint
	remainder int64
}

// Find values of x for which f(x) = a x^2 + b x + c factors (or almost factors) over fb.
func sievesmooth(a, b, c big.Int, fb []int64, rnd *rand.Rand) []sieveResult {
	var result []sieveResult

	maxp := fb[len(fb)-1]
	var bigmaxp2 big.Int
	bigmaxp2.SetInt64(maxp * maxp)
	fmt.Printf("maxp=%d maxp2=%d len(fb)=%d\n", maxp, &bigmaxp2, len(fb))

	// find x0 which is the minarg of ax^2+bx+c:
	//   d/dx(ax^2+bx+c) = 0
	//   2ax + b = 0
	//   x = -b/2a
	// so choose x0=-b/2a as the midpoint of our sieve

	var x0 big.Int
	x0.Div(&b, &a)
	x0.Rsh(&x0, 1)
	x0.Neg(&x0)
	fmt.Printf("x0:%d\n", &x0)

	// starting point
	x0.Sub(&x0, big.NewInt(sieverange))

	sieve := make([]byte, window)

	// temporaries
	var f, g, x, r, bigp, bigwindow big.Int
	bigwindow.SetInt64(window)
	var factors []uint

	// find starting points
	si := makeSieveInfo2(a, b, c, x0, fb, rnd)

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

		// check sieve entries for big numbers, indicating smooth f(x).
		f.Mul(&a, &x0)
		f.Add(&f, &b)
		f.Mul(&f, &x0)
		f.Add(&f, &c)
		fmt.Printf("x:%d f:%d\n", &x0, &f)

		x0.Add(&x0, &bigwindow)

		g.Mul(&a, &x0)
		g.Add(&g, &b)
		g.Mul(&g, &x0)
		g.Add(&g, &c)

		min := &f
		if g.Cmp(min) < 0 {
			min = &g
		}
		fmt.Printf("min: %d\n", min)
		threshold := byte(min.BitLen()) - 2*log2(maxp)
		fmt.Printf("threshold: %d\n", threshold)

		for j := 0; j < window; j++ {
			if sieve[j] < threshold {
				// TODO: in testing mode, check if f(x) is smooth.  If so, report a false negative.
				continue
			}

			// compute f(x)
			x.SetInt64(int64(j))
			x.Add(&x, &x0)
			f.Mul(&a, &x)
			f.Add(&f, &b)
			f.Mul(&f, &x)
			f.Add(&f, &c)

			// trial divide f by the factor base
			// accumulate factor base indexes of factors
			factors = factors[:0]
			for k, p := range fb {
				bigp.SetInt64(p)
				for r.Mod(&f, &bigp).Sign() == 0 {
					f.Div(&f, &bigp)
					factors = append(factors, uint(k))
				}
			}

			// if remainder > B^2, it's too big, might not be prime.
			if f.Cmp(&bigmaxp2) > 0 {
				//fmt.Printf("  false positive y=%d z=%d threshold=%d sieve[i]=%d log2(y)=%d log2(y/z)=%d\n", y, bigz, threshold, sieve[i], y.BitLen(), x.Div(y, bigz).BitLen())
				continue
			}

			var sr sieveResult
			sr.x.Set(&x)
			sr.factors = dup(factors)
			sr.remainder = f.Int64()
			result = append(result, sr)
		}
	}
	return result
}

type sieveinfo2 struct {
	pk   int32 // p^k for this factor base entry
	lg_p uint8 // ~log_2(p)
	off  int32 // working offset in sieve array
}

func makeSieveInfo2(a, b, c big.Int, start big.Int, fb []int64, rnd *rand.Rand) []sieveinfo2 {
	var si []sieveinfo2
	var biga, bigb, bigc, bigpk, m big.Int

	for _, p := range fb {
		pk := p
		for k := uint(1); k < 2; k++ { // TODO: quadraticModPK
			if pk > fb[len(fb)-1] {
				// Kind of arbitrary, but use powers of p as long as p^k is
				// smaller than than the maximum factor base prime.
				break
			}
			bigpk.SetInt64(pk)
			biga.Mod(&a, &bigpk)
			bigb.Mod(&b, &bigpk)
			bigc.Mod(&c, &bigpk)

			for _, r := range quadraticModP(biga.Int64(), bigb.Int64(), bigc.Int64(), pk, rnd) {
				// find first pk*i+r which is >= start
				s := m.Mod(&start, &bigpk).Int64()
				off := (r + pk - s) % pk
				si = append(si, sieveinfo2{int32(pk), log2(p), int32(off)})
				fmt.Printf("%#v\n", si[len(si)-1])
			}
			pk *= p
		}
	}
	return si
}
