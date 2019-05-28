package factorlib

import (
	"log"
	"math/rand"

	"github.com/randall77/factorlib/big"
	"github.com/randall77/factorlib/linear"
	"github.com/randall77/factorlib/primes"
)

func init() {
	factorizers["qs"] = qs
}

// The quadratic sieve is a general purpose factoring algorithm.
// It works by trying to find two numbers a, b satisfying:
//
//   a^2 == b^2 mod n with 0 < a, b < n   (1)
//
// Once we have found such a, b, we can rearrange equation (1) to get:
//   (a-b)(a+b) == 0 mod n
// A "good" pair of a, b have a != b and a+b != n.
// From any good pair satisfying (1) we can get a nontrivial factor of n
// by computing gcd(a+b, n). Why? If we have a good pair, then
//   a+b is not a multiple of n (a+b>0, a+b!=n, and a+b<2n).
//   a-b is not a multiple of n (a-b>-n, a-b!=0, and a-b<n).
// So the prime factors of n are nontrivially split between a+b and a-b.
//
// The chances of a random a, b satisfying (1) being "good" is at least ~1/2.
// So we only need to produce a constant number of a, b pairs before we expect
// to find a factor of n. (Not that we'll be producing a, b randomly, but
// as a practical matter we see the same behavior as if they were.)
//
// For instance, for n=527, we might get a=3 and b=65, because
// 3^2 == 65^2 == 9 mod 527.
// Then a factor of 527 is gcd(68, 527) = 17.
//
// So how do we find a, b satisfying (1)?
// First, we look for many x where x^2 - n factors into small primes.
//   x^2 - n = p1 · p2 · ... · pk
// Or equivalently,
//   x^2 == p1 · p2 · ... · pk mod n
// If we have the above equality for two different x, then we can
// generate another equality by multiplying each side of the two
// equalities (multiplying the x's on the left-hand side and adding
// the powers of matching primes on the right-hand side).
//
// If we can find a set of these equations where, when combined as
// above, all the powers of the primes on the right-hand side are
// even, then we can define a to be the product of all of the x's mod
// n, and b to be the square root of the product of all the primes on
// the right-hand side (the square root is easy to compute by just
// dividing all the prime powers by 2).

// For our 527 example, we'll get something like:
//   23^2 - 527 = 2
//   25^2 - 527 = 98 = 2 · 7 · 7
// combining those two equations,
//   23^2 · 25^2 == 2 · 2 · 7 · 7 mod 527
//   (23·25)^2 == (2·7)^2 mod 527
//   48^2 == 14^2 mod 527
// and sure enough, gcd(48+14, 527) == 31, which is a factor of 527.
//
// To find a subset of the x's whose right-hand sides combine to have
// all even prime powers, we use Gaussian elimination on a matrix over
// GF(2).
//
// It may seem somewhat recursive that in order to factor n, we need to
// first factor lots of x^2-n. Fortunately, even if factoring n is hard,
// factoring x^2-n need not be. Also, we don't have to factor all such x,
// only some of them. We can always give up on a specific x.
//
// We pick a set of small primes, called the factor base, over which
// we look for x^2-n that factor completely. You could imagine just
// using trial factorization with all the primes in the factor base on
// each trial x in turn, starting with x0 = ceil(sqrt(n)) and
// increasing by 1 each time.  That would work, but we can find x's
// where x^2-n factors over the factor base much faster using a sieve.
//
// A larger factor base means more of the x^2-n factor completely over
// the factor base. But the larger the factor base is, the more x we
// need to find before we have a good chance of finding a subset that
// has only even prime powers on the right-hand side.
//
// If, when we divide out all the primes in the factor base, we are left
// with a number smaller than the square of the maximum factor base prime,
// then that value is also a prime factor of x^2-n. So we can get effectively
// square the maximum factor base prime using this technique (or equivalently,
// sieve with only primes up to the square root of the max factor base prime).
func qs(n big.Int, rnd *rand.Rand, logger *log.Logger) ([]big.Int, error) {
	// qs does not work for powers of a single prime.  Check that first.
	if f, err := primepower(n, rnd, logger); err == nil {
		return f, nil
	}

	// first, pick a factor base
	fb, a := makeFactorBase(n)
	if a != 0 {
		return []big.Int{big.Int64(a), n.Div64(a)}, nil
	}

	// matrix is used to do gaussian elimination on mod 2 exponents.
	m := linear.NewMatrix(uint(len(fb)))

	// pair up large primes that we find using this table
	type largerecord struct {
		x big.Int
		f []uint
	}
	largeprimes := map[int64]largerecord{}

	x0 := n.SqrtCeil()
	for {
		//logger.Printf("sieving at %d\n", x0)
		for _, r := range sievesmooth(big.Int64(1), big.Int64(0), n.Neg(), fb, x0, rnd) {
			x := r.x
			factors := r.factors
			remainder := r.remainder
			/*
				fmt.Printf("%d^2-%d=%d=", x, n, x.Mul(x).Sub(n))
				for i, f := range factors {
					if i != 0 {
						fmt.Printf("·")
					}
					fmt.Printf("%d", fb[f])
				}
				if remainder != 1 {
					fmt.Printf("·%d", remainder)
				}
				fmt.Println()
			*/
			if remainder != 1 {
				// try to find another record with the same largeprime
				lr, ok := largeprimes[remainder]
				if !ok {
					// haven't seen this large prime yet.  Save record for later
					largeprimes[remainder] = largerecord{x, factors}
					continue
				}
				// combine current equation with other largeprime equation
				// x1^2 === prod(f1) * largeprime
				// x2^2 === prod(f2) * largeprime
				//fmt.Printf("  largeprime %d match\n", remainder)
				x = x.Mul(lr.x).Mod(n).Mul(big.Int64(remainder).ModInv(n)).Mod(n) // TODO: could remainder divide n?
				factors = append(factors, lr.f...)
			}

			// Add equation to the matrix
			idlist := m.AddRow(factors, eqn{x, factors})
			if idlist == nil {
				if m.Rows()%100 == 0 {
					logger.Printf("%d/%d falsepos=%d largeprimes=%d\n", m.Rows(), len(fb), falsepos, len(largeprimes))
					falsepos = 0
				}
				continue
			}

			// We found a set of equations with all even powers.
			// Compute a and b where a^2 === b^2 mod n
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
			for _, o := range odd {
				if o {
					panic("gauss elim failed")
				}
			}

			if a.Cmp(b) == 0 {
				// trivial equation, ignore it
				logger.Println("triv A")
				continue
			}
			if a.Add(b).Cmp(n) == 0 {
				// trivial equation, ignore it
				logger.Println("triv B")
				continue
			}

			r := a.Add(b).GCD(n)
			return []big.Int{r, n.Div(r)}, nil
		}
		x0 = x0.Add64(sieverange)
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

func dup(x []uint) []uint {
	y := make([]uint, len(x))
	copy(y, x)
	return y
}
