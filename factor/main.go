package main

import (
	"flag"
	"fmt"
	"github.com/randall77/factorlib"
	"github.com/randall77/factorlib/big"
	"math/rand"
	"os"
	"strconv"
)

var seed = flag.Int64("seed", 0, "seed for RNG")
var alg = flag.String("alg", "trial", "factoring algorithm to use")

func main() {
	flag.Parse()

	rnd := rand.New(rand.NewSource(*seed))

	// Figure out the number to factor
	args := flag.Args()
	if len(args) == 0 {
		fmt.Println("no number to factor")
		os.Exit(1)
	}
	if len(args) > 1 {
		fmt.Println("can't factor multiple numbers")
		os.Exit(1)
	}
	var n big.Int
	nstr := args[0]
	switch nstr[0] {
	case 'r':
		// random d-digit number
		d, err := strconv.Atoi(nstr[1:])
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
		k := big.Ten.Exp(int64(d) - 1)
		n = k.Mul64(9).Rand(rnd).Add(k)
	case 's':
		// random d-digit semiprime
		d, err := strconv.Atoi(nstr[1:])
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
		if d%2 != 0 {
			fmt.Println("semiprime must have an even number of digits")
			os.Exit(1)
		}
		min := big.Ten.Exp(int64(d) - 1)
		max := min.Mul64(10)
		for {
			x := randomPrime(d/2, rnd)
			y := randomPrime(d/2, rnd)
			n = x.Mul(y)
			if n.Cmp(min) >= 0 && n.Cmp(max) < 0 {
				break
			}
		}
	default:
		var ok bool
		n, ok = big.ParseInt(nstr)
		if !ok {
			fmt.Printf("parsing \"%s\": invalid number\n", nstr)
			os.Exit(1)
		}
		if n.Cmp(big.One) <= 0 {
			fmt.Printf("invalid n: %s\n", nstr)
			os.Exit(1)
		}
	}

	// Call into main library to do factoring
	fmt.Printf("factoring %d using algorithm %s\n", n, *alg)
	factors := factorlib.Factor(n, *alg, rnd)

	// Print result
	fmt.Printf("%d = ", n)
	for i, f := range factors {
		if i > 0 {
			fmt.Print("·")
		}
		fmt.Printf("%d", f)
	}
	fmt.Println()
}

// make a random prime with the given number of digits
func randomPrime(digits int, rnd *rand.Rand) big.Int {
	min := big.Ten.Exp(int64(digits - 1))
	w := min.Mul64(9)
	for {
		n := min.Add(w.Rand(rnd))
		if n.ProbablyPrime(1000) {
			return n
		}
	}
}
