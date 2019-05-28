package factorlib

import (
	"fmt"
	"math/rand"
	"sort"
	"strings"

	"github.com/randall77/factorlib/big"
)

// Set of factoring algorithms to choose from.
// Algorithms can add themselves here in an initializer.
// Eventually, we should choose one automagically.
var factorizers = map[string]func(big.Int, *rand.Rand) ([]big.Int, error){}

// Factor returns the prime factorization of n.
// alg is a hint about which factoring algorithm to choose.
// rnd is a random source for use by the factoring algorithm.
func Factor(n big.Int, alg string, rnd *rand.Rand) ([]big.Int, error) {
	// figure out the algorithm to use
	split := factorizers[alg]
	if split == nil {
		var algs []string
		for k, _ := range factorizers {
			algs = append(algs, k)
		}
		sort.Strings(algs)
		return nil, fmt.Errorf("unknown algorithm: %s (possible algorithms: %s)", alg, strings.Join(algs, ", "))
	}

	// Main loop.  Keep splitting factors until they are all prime.
	factors := []big.Int{n}
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
		d, err := split(f, rnd)
		if err != nil {
			return nil, err
		}

		// check answer
		m := big.One
		for _, x := range d {
			if x.Cmp(big.One) == 0 || x.Cmp(f) == 0 {
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
	m := big.One
	for _, f := range factors {
		m = m.Mul(f)
	}
	if m.Cmp(n) != 0 {
		panic("factorization failed")
	}
	return factors, nil
}

// sort *big.Int in nondecreasing order
type byValue []big.Int

func (a byValue) Len() int           { return len(a) }
func (a byValue) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a byValue) Less(i, j int) bool { return a[i].Cmp(a[j]) < 0 }
