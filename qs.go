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

	// first, pick a factor base
	fb, a := makeFactorBase(n)
	if a != 0 {
		return []big.Int{big.Int64(a), n.Div64(a)}
	}

	// matrix is used to do gaussian elimination on mod 2 exponents.
	m := linear.NewMatrix(uint(len(fb)))

	// pair up large primes that we find using this table
	type largerecord struct {
		x big.Int
		f []uint
	}
	largeprimes := map[int64]largerecord{}

	// function to process sieve results
	fn := func(x big.Int, factors []uint, remainder int64) []big.Int {
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
				return nil
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
			if m.Rows() % 100 == 0 {
				fmt.Printf("%d/%d\n", m.Rows(), len(fb))
			}
			return nil
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
			fmt.Println("triv A")
			return nil
		}
		if a.Add(b).Cmp(n) == 0 {
			// trivial equation, ignore it
			fmt.Println("triv B")
			return nil
		}

		r := a.Add(b).GCD(n)
		return []big.Int{r, n.Div(r)}
	}

	x0 := n.SqrtCeil()
	for {
		fmt.Printf("sieving at %d\n", x0)
		r := sievesmooth2(big.Int64(1), big.Int64(0), n.Neg(), fb, rnd, x0, fn)
		if r != nil {
			return r
		}
		x0 = x0.Add64(2*sieverange)
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
