package sieve

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"unsafe"

	"github.com/randall77/factorlib/big"
	fmath "github.com/randall77/factorlib/math"
)

// use an array of this size to do the sieving
// 16KB fits in L1 cache (L1 cache is typically 32KB)
const windowBits = 14
const window = 1 << windowBits

// check to see if anything usable isn't returned by the sieve
const checkFalseNegative = false

// Records f(x) == product(factors)*remainder
// The values in the factors slice are indexes into the factor base
type Result struct {
	X         big.Int
	Factors   []uint
	Remainder int64
}

// Smooth finds values of x for which f(x) = a x^2 + b x + c factors (within one bigprime) over the primes in fb.
// Returns each x found, togegther with the factorization of f(x) into factor base primes and a
// possible bigprime (< max(fb)^2) remainder.
// Searches in the window [x0, x1).
// Also returns the number of false positives (for logging only)
func Smooth(a, b, c big.Int, fb []int64, x0, x1 big.Int, rnd *rand.Rand) ([]Result, int64) {
	if fb[0] != -1 {
		panic("first factor base must be -1")
	}
	// Compute result.
	result, falsePos := smooth(a, b, c, fb, x0, x1, rnd)
	// Check result.
	for _, r := range result {
		f := a.Mul(r.X).Add(b).Mul(r.X).Add(c)
		g := big.Int64(r.Remainder)
		for _, i := range r.Factors {
			g = g.Mul64(fb[i])
		}
		if !f.Equals(g) {
			panic(fmt.Sprintf("Smooth result not correct a=%d b=%d c=%d f(%d)=%d != %d %v %d", a, b, c, r.X, f, g, r.Factors, r.Remainder))
		}
	}
	// Return result.
	return result, falsePos
}

func smooth(a, b, c big.Int, fb []int64, x0, x1 big.Int, rnd *rand.Rand) ([]Result, int64) {
	if x0.Equals(x1) {
		return nil, 0
	}
	if x0.Add64(1).Equals(x1) {
		// Testing just x0.
		var scratch big.Scratch
		var factors []uint
		f := a.Mul(x0).Add(b).Mul(x0).Add(c)
		if f.IsZero() {
			return nil, 0
		}
		for i, p := range fb {
			if p == -1 {
				if f.Sign() < 0 {
					f = f.Neg()
					factors = append(factors, uint(i))
				}
				continue
			}
			if f.Mod64s(p, &scratch) == 0 {
				f = f.Div64(p)
				factors = append(factors, uint(i))
			}
		}
		maxp := fb[len(fb)-1] // TODO: assumes last entry in factor base is the largest
		if f.Cmp64(maxp*maxp) < 0 {
			return []Result{Result{X: x0, Factors: factors, Remainder: f.Int64()}}, 0
		}
		return nil, 0
	}

	// Find min and max of f(x)=ax^2+bx+c over the range [x0,x1)
	// The min and max must be either one of the endpoints, or
	// the value at the extremum x=-b/2a, if that is in range.
	y0 := a.Mul(x0).Add(b).Mul(x0).Add(c)
	min := y0
	max := y0
	x1m1 := x1.Sub64(1)
	y1 := a.Mul(x1m1).Add(b).Mul(x1m1).Add(c)
	if max.Cmp(y1) < 0 {
		max = y1
	}
	if min.Cmp(y1) > 0 {
		min = y1
	}
	x := b.Div(a.Lsh(1)).Neg() // TODO: check both round up and round down?
	if x.Cmp(x0) > 0 && x.Cmp(x1m1) < 0 {
		y := a.Mul(x).Add(b).Mul(x).Add(c)
		if max.Cmp(y) < 0 {
			max = y
		}
		if min.Cmp(y) > 0 {
			min = y
		}
	}

	if max.Sign() != min.Sign() {
		// split if we change sign
		// TODO: just binary search for the midpoint? Or is this enough?
		xmid := x0.Add(x1).Rsh(1)
		r0, f0 := smooth(a, b, c, fb, x0, xmid, rnd)
		r1, f1 := smooth(a, b, c, fb, xmid, x1, rnd)
		return append(r0, r1...), f0 + f1
	}

	// sieve to find any potential smooth f(x)
	res := sieveinner(a, b, c, fb, x0, x1, min, max, rnd)

	// results buffer
	var result []Result
	var factors []uint
	s := &big.Scratch{}

	// check potential results using trial factorization
	var falsePos int64
	for _, x := range res {
		// compute y=f(x)
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
		maxp := fb[len(fb)-1]
		if y.Cmp64(maxp*maxp) > 0 {
			falsePos++
			continue
		}
		if y.Cmp(big.Int64(y.Int64())) != 0 {
			continue
		}

		result = append(result, Result{X: x, Factors: dup(factors), Remainder: y.Int64()})
		if len(result) > 2*len(fb) {
			// early out when we're factoring small n
			// TODO: include in spec for Smooth function?
			return result, falsePos
		}
	}
	if checkFalseNegative {
		// This is expensive, we have to trial divide all the numbers in the
		// sieve range.  Oh well, testing.
	checkloop:
		for x := x0; x.Cmp(x1) < 0; x = x.Add64(1) {
			for _, r := range res {
				if r.Equals(x) {
					continue checkloop
				}
			}
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
			maxp := fb[len(fb)-1]
			if y.Cmp64(maxp*maxp) <= 0 {
				log.Printf("  false negative x=%d f(x)=%d f=%v·%d\n", x, a.Mul(x).Add(b).Mul(x).Add(c), f, y)
			}
		}
	}
	return result, falsePos
}

// sieveEntry is the data type in which we accumulate the sum of
// the scaled logarithms of the factors we have found so far.
type sieveEntry uint8

// sieveInfo contains data about a power of a factor base prime.
type sieveInfo uint64

const (
	pkBits   = 32
	offBits  = windowBits
	logpBits = 8 * unsafe.Sizeof(sieveEntry(0))
)

func makeSieveInfo(pk uint, off uint, logp sieveEntry) sieveInfo {
	return sieveInfo(pk) + sieveInfo(off)<<pkBits + sieveInfo(logp)<<(pkBits+offBits)
}

// p^k, the amount we step by through the sieve array
func (si sieveInfo) pk() uint {
	return uint(si & (1<<pkBits - 1))
}

// offset within the window which holds the next muliple of p^k
func (si sieveInfo) off() uint {
	return uint(si >> pkBits & (1<<offBits - 1))
}

// scaled value of log(p), the amount we add to each sieve array slot
func (si sieveInfo) logp() sieveEntry {
	return sieveEntry(si >> (pkBits + offBits))
}

func init() {
	if pkBits+offBits+logpBits > 64 {
		panic("sieveInfo needs too many bits")
	}
	if logpBits != 8 {
		panic("bad")
	}
}

func sieveinner(a, b, c big.Int, fb []int64, x0, x1, min, max big.Int, rnd *rand.Rand) []big.Int {
	// Compute scaling factor for logarithms
	// This guarantees that the sieve values will not overflow
	var scale float64
	if max.Sign() > 0 {
		scale = 255.99 / max.Log()
	} else {
		scale = 255.99 / min.Neg().Log()
	}

	s := &big.Scratch{}
	maxp := fb[len(fb)-1]
	logmaxp := big.Int64(maxp).Log()

	// Our sieve will be an array of bytes of size "window".
	//
	// We divide the factor base into "small" and "large" factors.
	// Small factors are <window, so for them we need to mark at least
	// one entry in every window.Large factors are >window, where we don't
	// need to update the window every time for that factor.
	//
	// Small factors are kept in a slice and processed every window.
	// Large factors are kept in a slice of slices. The outer slice
	// is indexed by how many windows away the next occurence of the
	// prime powers in that slice are.

	// List of sieveInfo records for small primes
	var small []sieveInfo
	// List of lists of sieveInfo records for large primes. Indexed by the number
	// of windows away the next occurrence of that prime power is.
	large := make([][]sieveInfo, (maxp+window-1)/window+1)
	for _, p := range fb[1:] {
		if a.Mod64s(p, s) == 0 {
			// This can happen in mpqs if p is one of the factors we used to construct a.
			// Ignore this prime - we might miss a few smooth f(x), but such is life.
			continue
		}
		logp := sieveEntry(scale * math.Log(float64(p))) // note: rounds down to avoid overflow in sieveEntry type
		if logp == 0 {
			// Not the end of the world, but it means we're not doing
			// anything useful for this factor base prime.
			panic("logp too small")
		}
		// The upper bound is kind of arbitrary, but choose to use powers of p
		// as long as p^k is smaller than than the maximum factor base prime.
		for k, pk := uint(1), p; pk <= maxp; k, pk = k+1, pk*p {
			st := x0.Mod64s(pk, s)
			for _, r := range fmath.QuadraticModPK(a.Mod64s(pk, s), b.Mod64s(pk, s), c.Mod64s(pk, s), p, k, pk, rnd) {
				// find first pk*i+r which is >= x0
				off := (r - st + pk) % pk
				if pk < window {
					small = append(small, makeSieveInfo(uint(pk), uint(off), logp))
				} else {
					i := off / window
					large[i] = append(large[i], makeSieveInfo(uint(pk), uint(off%window), logp))
				}
			}
		}
	}

	// A smooth number is one for which we've found a lot of small
	// factors. This will be identifiable in the sieve array by a
	// large sieve value, representing the sum of the logarithms of the
	// factor base primes that divide f(x0 + sieve index).
	//
	// The threshold is the value the sieve array entry needs to reach
	// before we become interested in it. The max entry is 255, and if we
	// reach that value, we know we have fully factored f(x) over the
	// factor base. In practice, we don't need to reach all the way to
	// 255 for the entry to be useful.
	// The 7/8 is a fudge factor to account for rounding and other
	// second order effects.
	// Using min.Log accounts for non-uniform size of f(x).
	// 2*log(maxp) allows us to find a single prime remainder
	// less than the square of the max factor base prime.

	// This is the array of sums of logarithms of factors found.
	// Factors of f(x0+i) are added to sieve[i%window]
	// (We do a range of window numbers at a time.)
	var sieve [window]sieveEntry

	// Results we've found so far (the x values).
	var r []big.Int

	for x0.Cmp(x1) < 0 {
		// start with a clear sieve
		for j := range sieve {
			sieve[j] = 0
		}

		// increment sieve entries for f(x) that are divisible
		// by each factor base prime.

		// first, the small primes
		for i, s := range small {
			pk := s.pk()
			off := s.off()
			logp := s.logp()
			for ; off < window; off += pk {
				sieve[off] += logp
			}
			small[i] = makeSieveInfo(pk, off-window, logp)
		}

		// second, the large primes
		l := large[0]
		copy(large, large[1:])
		large[len(large)-1] = l[:0] // reuses slice backing array
		for _, s := range l {
			pk := s.pk()
			off := s.off()
			logp := s.logp()
			sieve[off] += logp
			off += pk
			j := off / window
			large[j] = append(large[j], makeSieveInfo(pk, off%window, logp))
		}

		// Check for smooth numbers.
		// Use the threshold corresponding to the smaller endpoint.
		y0 := a.Mul(x0).Add(b).Mul(x0).Add(c).Abs()
		z := x0.Add64(window)
		y1 := a.Mul(z).Add(b).Mul(z).Add(c).Abs()
		if y0.Cmp(y1) > 0 {
			y0 = y1
		}
		threshold := sieveEntry(scale * (y0.Abs().Log() - 2*logmaxp))
		for i := 0; i < window; i++ {
			if sieve[i] >= threshold {
				x := x0.Add64(int64(i))
				if x.Cmp(x1) < 0 {
					r = append(r, x)
				}
			}
		}

		// Next window.
		x0 = x0.Add64(window)
	}
	return r
}

func dup(x []uint) []uint {
	y := make([]uint, len(x))
	copy(y, x)
	return y
}
