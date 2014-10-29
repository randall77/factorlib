package factorlib

import (
	"fmt"
	"math/big"
	"math/rand"
	"sort"
)

// helpful constants
var zero = *big.NewInt(0)
var one = *big.NewInt(1)
var two = *big.NewInt(2)

// Set of factoring algorithms to choose from.
// Algorithms can add themselves here in an initializer.
// Eventually, we should choose one automagically.
var factorizers = map[string]func(big.Int, *rand.Rand) []big.Int{}

// Factor returns the prime factorization of n.
// alg is a hint about which factoring algorithm to choose.
// rnd is a random source for use by the factoring algorithm.
func Factor(n big.Int, alg string, rnd *rand.Rand) []big.Int {
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
	var m big.Int
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
		d := split(f, rnd)

		// check answer
		m.SetInt64(1)
		for _, x := range d {
			if x.Cmp(&one) == 0 || x.Cmp(&f) == 0 {
				panic("trivial factor")
			}
			m.Mul(&m, &x)
		}
		if m.Cmp(&f) != 0 {
			panic("splitter failed")
		}

		// update list of factors
		factors[i] = d[0]
		factors = append(factors, d[1:]...)
	}

	// Sort for nice display.
	sort.Sort(byValue(factors))

	// Check final result.
	m.SetInt64(1)
	for _, f := range factors {
		m.Mul(&m, &f)
	}
	if m.Cmp(&n) != 0 {
		panic("factorization failed")
	}
	return factors
}

// sort *big.Int in nondecreasing order
type byValue []big.Int

func (a byValue) Len() int           { return len(a) }
func (a byValue) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a byValue) Less(i, j int) bool { return a[i].Cmp(&a[j]) < 0 }
