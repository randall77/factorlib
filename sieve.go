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

// Find values of x for which f(x) = a x^2 + b x + c factors (within one bigprime) over the primes in fb.
// requires: a > 0
func sievesmooth(a, b, c bigint, fb []int64, rnd *rand.Rand) []sieveResult {
	var result []sieveResult

	maxp := fb[len(fb)-1]

	// find approximate zero crossings
	d := b.Square().Sub(a.Mul(c).Lsh(2))
	if d.Sign() < 0 {
		panic("polynomial has no roots")
		// TODO: choose min instead?  Then x = -b/2a
	}
	x := b.Neg().Add(d.SqrtFloor()).Div(a).Rsh(1)
	//x2 := b.Neg().Sub(d).Div(a).Rsh(1)
	// TODO: sieve around x2 also? (if d != 0)

	// starting point
	x0 := x.Sub64(sieverange)

	// results buffer
	var factors []uint

	// find starting points
	si := makeSieveInfo2(a, b, c, x0, fb, rnd)

	// pick threshold
	threshold := byte(a.Mul(x0).Add(b).Mul(x0).Add(c).BitLen()) - 2*log2(maxp) // TODO: subtract more?
	
	// sieve to find any potential smooth f(x)
	sieve := make([]byte, window) // TODO: cache this?
	res := sieveinner(sieve, si, threshold)
	
	// check potential results using trial factorization
	for _, i := range res {
		// compute f(x)
		x := x0.Add64(int64(i))
		y := a.Mul(x).Add(b).Mul(x).Add(c)
		
		// trial divide f by the factor base
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
			for y.Mod64(p) == 0 {
				y = y.Div64(p)
				factors = append(factors, uint(k))
			}
		}
		
		// if remainder > B^2, it's too big, might not be prime.
		if y.Cmp64(maxp*maxp) > 0 {
			//fmt.Printf("  false positive y=%d z=%d threshold=%d sieve[i]=%d log2(y)=%d log2(y/z)=%d\n", y, bigz, threshold, sieve[i], y.BitLen(), x.Div(y, bigz).BitLen())
			continue
		}
		
		result = append(result, sieveResult{x, dup(factors), y.Int64()})
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

	for _, p := range fb[1:] {
		pk := p
		for k := uint(1); k < 2; k++ { // TODO: quadraticModPK
			if pk > fb[len(fb)-1] {
				// Kind of arbitrary, but use powers of p as long as p^k is
				// smaller than than the maximum factor base prime.
				break
			}
			s := start.Mod64(pk)
			for _, r := range quadraticModP(a.Mod64(pk), b.Mod64(pk), c.Mod64(pk), pk, rnd) {
				// find first pk*i+r which is >= start
				off := (r - s + pk) % pk
				si = append(si, sieveinfo2{int32(pk), log2(p), int32(off)})
			}
			pk *= p
		}
	}
	return si
}

func init() {
	factorizers["qs2"] = qs2
}

func qs2(n bigint, rnd *rand.Rand) []bigint {
	// qs does not work for powers of a single prime.  Check that first.
	if f := primepower(n, rnd); f != nil {
		return f
	}

	// first, pick a factor base
	fb, a := makeFactorBase(n)
	if a != 0 {
		return []bigint{Big(a), n.Div64(a)}
	}

	for _, r := range sievesmooth(Big(1), Big(0), n.Neg(), fb, rnd) {
		fmt.Printf("f(%d)= prod %v * %d\n", r.x, r.factors, r.remainder)
	}
	return nil
}
