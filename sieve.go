package factorlib

import (
	"fmt"
	"github.com/randall77/factorlib/big"
	"math/rand"
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
	sieve := make([]byte, window) // TODO: cache this?
	res := sieveinner(sieve, si, thresholds[:])

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

// TODO: write this in assembly?
func sieveinner(sieve []byte, si []sieveinfo2, thresholds []byte) []int {
	var r []int
	for i := 0; i < sieverange; i += window {
		// clear sieve
		for j := 0; j < window; j++ {
			sieve[j] = 0
		}
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
