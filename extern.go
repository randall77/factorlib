package factorlib

import (
	"fmt"
	"math/big"
	"math/rand"
	"sort"
)

// Set of factoring algorithms to choose from.
// Algorithms can add themselves here in an initializer.
// Eventually, we should choose one automagically.
var factorizers = map[string]func(bigint, *rand.Rand) []bigint{}

// Factor returns the prime factorization of n.
// alg is a hint about which factoring algorithm to choose.
// rnd is a random source for use by the factoring algorithm.
func Factor(n_ *big.Int, alg string, rnd *rand.Rand) []*big.Int {
	n := NewBigFromBig(n_)
	// figure out the algorithm to use
	split := factorizers[alg]
	if split == nil {
		// TODO: better error
		var algs []string
		for k, _ := range factorizers {
			algs = append(algs, k)
		}
		sort.Strings(algs)
		fmt.Printf("unknown algorithm: %s\n", alg)
		fmt.Println("possible algorithms:")
		for _, s := range algs {
			fmt.Printf("  %s\n", s)
		}
		return nil
	}

	// Main loop.  Keep splitting factors until they are all prime.
	factors := []bigint{n}
	for i := 0; i < len(factors); {
		f := factors[i]

		// If the current factor is prime, leave it in the list.
		if f.ProbablyPrime(1000) {
			i++
			continue
		}

		// Otherwise, find a nontrivial factorization of it.
		// This is the main call from the driver into
		// the specific factoring algorithm chosen.
		d := split(f, rnd)

		// check answer
		m := one
		for _, x := range d {
			if x.Cmp(one) == 0 || x.Cmp(f) == 0 {
				panic("trivial factor")
			}
			m = m.Mul(x)
		}
		if m.Cmp(f) != 0 {
			panic("splitter failed")
		}

		// update list of factors
		factors[i] = d[0]
		factors = append(factors, d[1:]...)
	}

	// Sort for nice display.
	sort.Sort(byValue(factors))

	// Check final result.
	m := one
	for _, f := range factors {
		m = m.Mul(f)
	}
	if m.Cmp(n) != 0 {
		panic("factorization failed")
	}
	var factors_ []*big.Int
	for _, f := range factors {
		factors_ = append(factors_, f.GetBig())
	}
	return factors_
}

// sort *big.Int in nondecreasing order
type byValue []bigint

func (a byValue) Len() int           { return len(a) }
func (a byValue) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a byValue) Less(i, j int) bool { return a[i].Cmp(a[j]) < 0 }
