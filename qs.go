package factorlib

import (
	"fmt"
	"github.com/randall77/factorlib/big"
	"github.com/randall77/factorlib/linear"
	"github.com/randall77/factorlib/primes"
	"math/rand"
)

func init() {
	factorizers["qs"] = qs
}

func qs(n big.Int, rnd *rand.Rand) []big.Int {
	// qs does not work for powers of a single prime.  Check that first.
	if f := primepower(n, rnd); f != nil {
		return f
	}

	// Quadratic sieve:
	// Define f(x) = x^2 - n.  We wish to find x such that f(x) is smooth.
	// We will take x starting at ceil(sqrt(n)).

	// is f(x) divisible by p^k?
	// calculating mod p^k, f(x + p^k) = (x+p^k)^2 - n = x^2 - n = f(x)
	// so if f(x) is divisible by p^k, so is f(x + i * p^k) for all i.
	//
	// if x_1 and x_2 are the two solutions to f(x) == 0 mod p^k, then all x which
	// have f(x) divisible by p^k are x_1 + i * p^k and x_2 + i * p^k.
	//
	// f(x) == 0 mod p^k is x^2 == n mod p^k.

	// So the factor base contains small primes and powers of primes.
	// For each factor base element p^k, we store two offsets satisfying
	// f(start + off1) == 0 mod p^k and f(start + off2) == 0 mod p^k.
	// Then to sieve, we just add log_p to sieve[off1], sieve[off1+p^k], sieve[off1+2*p^k], ...
	// and the same for off2.  When we reach the end of the sieve, we update off1 and off2
	// for the next sieve.

	// the sieve returns x such that f(x) = x^2 - n = prod_i factors[i]
	// all factors <= B
	// factors are sorted in increasing order.

	// first, pick a factor base
	fb, a := makeFactorBase(n)
	if a != 0 {
		return []big.Int{big.Int64(a), n.Div64(a)}
	}
	fb = fb[1:]

	maxp := fb[len(fb)-1]
	fmt.Printf("maxp=%d len(fb)=%d\n", maxp, len(fb))

	// The x we start with
	start := n.SqrtCeil()

	si := makeSieveInfo(n, start, fb, rnd)

	sievelen := int(maxp)
	if sievelen < 10000 {
		sievelen = 10000
	}

	sieve := make([]byte, sievelen)

	// factors we find
	var factors []uint

	// matrix is used to do gaussian elimination on mod 2 exponents.
	m := linear.NewMatrix(uint(len(fb)))

	// largeprimes records instances of the equation
	//   x^2 = prod(f) * p mod n
	// where f are small primes (indexes into factor base) and
	// p is a large prime.
	type largerecord struct {
		x big.Int
		f []uint
	}
	largeprimes := map[int64]largerecord{}

	s := &big.Scratch{}

	for {
		// clear sieve
		for i := 0; i < sievelen; i++ {
			sieve[i] = 0
		}

		// increment sieve entries for f(x) that are divisible
		// by each factor base prime.
		for i := range si {
			f := &si[i]
			pk := int(f.pk)
			lg_p := f.lg_p
			j := int(f.off1)
			for ; j < sievelen; j += pk {
				sieve[j] += lg_p
			}
			f.off1 = int32(j - sievelen) // for next time

			// same for off2
			j = int(f.off2)
			if j < 0 {
				continue // special case needed for p==2,k==1
			}
			for ; j < sievelen; j += pk {
				sieve[j] += lg_p
			}
			f.off2 = int32(j - sievelen) // for next time
		}

		// check sieve entries for big numbers, indicating smooth f(x).
		threshold := byte(start.Square().Sub(n).BitLen()) - 2*log2(maxp)
		for i := 0; i < sievelen; i++ {
			if sieve[i] < threshold {
				// TODO: in testing mode, check if x^2-n is smooth.
				continue
			}
			//fmt.Printf("%d+%d %d %d\n", start, i, sieve[i], threshold)

			// compute x, y=f(x)
			x := start.Add64(int64(i))
			y := x.Square().Sub(n)

			// trial divide y by the factor base
			// accumulate factor base indexes of factors
			factors = factors[:0]
			for i, p := range fb {
				for y.Mod64s(p, s) == 0 {
					y = y.Div64(p)
					factors = append(factors, uint(i))
				}
			}

			// if remainder > B^2, it's too big
			if y.Cmp64(maxp*maxp) > 0 {
				//fmt.Printf("  false positive y=%d z=%d threshold=%d sieve[i]=%d log2(y)=%d log2(y/z)=%d\n", y, bigz, threshold, sieve[i], y.BitLen(), x.Div(y, bigz).BitLen())
				continue
			}

			z := y.Int64()
			if z > 1 {
				// try to find another record with the same largeprime
				lr, ok := largeprimes[z]
				if !ok {
					var lr largerecord
					lr.x = x
					lr.f = dup(factors)
					// haven't seen this large prime yet.  Save record for later
					largeprimes[z] = lr
					//fmt.Printf("  savelarge %d %v\n", z, factors)
					continue
				}
				// combine current equation with other largeprime equation
				// x1^2 === prod(f1) * largeprime
				// x2^2 === prod(f2) * largeprime
				//fmt.Printf("  largeprime %d match\n", bigz)
				x = x.Mul(lr.x).Mod(n).Mul(y.ModInv(n)).Mod(n) // TODO: could y divide n?
				factors = append(factors, lr.f...)
			}
			fmt.Printf("eqn%d/%d %d^2 === ", m.Rows(), len(fb), x)
			for j, i := range factors {
				if j > 0 {
					fmt.Printf("·")
				}
				fmt.Printf("%d", fb[i])
			}
			fmt.Printf("\n")
			var e eqn
			e.x = x
			e.f = dup(factors)
			idlist := m.AddRow(factors, e)
			if idlist == nil {
				continue
			}

			// we found a set of equations with all even powers
			// compute a and b where a^2 === b^2 mod n
			a := big.One
			b := big.One
			odd := make([]bool, len(fb))
			for _, id := range idlist {
				e := id.(eqn)
				a = a.Mul(e.x).Mod(n)
				for _, i := range e.f {
					if !odd[i] {
						// first occurrence of this factor
						odd[i] = true
						continue
					}
					// second occurrence of this factor
					b = b.Mul64(fb[i]).Mod(n)
					odd[i] = false
				}
			}

			for i, p := range fb {
				if odd[i] {
					fmt.Printf("prime i=%d p=%d\n", i, p)
					panic("gauss elim failed")
				}
			}

			if a.Cmp(b) == 0 {
				// trivial equation, ignore it
				fmt.Println("triv A")
				continue
			}
			if a.Add(b).Cmp(n) == 0 {
				// trivial equation, ignore it
				fmt.Println("triv B")
				continue
			}

			r := a.Add(b).GCD(n)
			return []big.Int{r, n.Div(r)}
		}

		start = start.Add64(int64(sievelen))
	}
}

// x^2 === prod(f) mod n
// f is a list of indexes into the factor base
type eqn struct {
	x big.Int
	f []uint
}

// pick some prime factors for our factor base.  If we happen
// upon a factor of n, return it instead.
func makeFactorBase(n big.Int) ([]int64, int64) {
	// upper limit on prime factors (TODO: dependent on n) that we sieve with
	const B = 50000
	fb := []int64{-1}
	s := &big.Scratch{}
	for i := 0; ; i++ {
		p := primes.Get(i)
		if p > B {
			return fb, 0
		}
		a := n.Mod64s(p, s)
		if a == 0 {
			return nil, p
		}
		if quadraticResidue(a, p) {
			// if x^2 == n mod p has no solutions, then
			// x^2-n will never be divisible by p.  So it
			// doesn't need to be in the factor base.
			fb = append(fb, p)
		}
	}
}

type sieveinfo struct {
	pk         int32 // p^k for this factor base entry
	lg_p       uint8 // ~log_2(p)
	off1, off2 int32 // starting offsets in sieve array
}

func makeSieveInfo(n big.Int, start big.Int, fb []int64, rnd *rand.Rand) []sieveinfo {
	var si []sieveinfo
	s := &big.Scratch{}

	for _, p := range fb {
		pk := p
		for k := uint(1); ; k++ {
			if pk > fb[len(fb)-1] {
				break
			}
			z1 := sqrtModPK(n.Mod64s(pk, s), p, k, rnd) // find solution to x^2 == a mod p^k
			if z1 == 0 {
				panic("factor base divides n")
			}
			z2 := pk - z1

			// adjust to start
			s := start.Mod64s(pk, s)
			off1 := z1 - s
			if off1 < 0 {
				off1 += pk
			}
			off2 := z2 - s
			if off2 < 0 {
				off2 += pk
			}
			if off1 == off2 {
				// should only happen for pk==2
				if pk != 2 {
					panic("duplicate offsets")
				}
				off2 = -1
			}

			si = append(si, sieveinfo{int32(pk), log2(p), int32(off1), int32(off2)})
			if p == 2 && k == 1 {
				// TODO: handle powers of 2 correctly.  Hansel's lemma
				// doesn't work for them.
				break
			}
			pk *= p
		}
	}
	return si
}

func dup(x []uint) []uint {
	y := make([]uint, len(x))
	copy(y, x)
	return y
}
