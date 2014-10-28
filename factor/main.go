package main

import (
	"github.com/randall77/factorlib"
	"flag"
	"fmt"
	"math/big"
	"math/rand"
	"os"
	"strconv"
)

var seed = flag.Int64("seed", 0, "seed for RNG")
var alg = flag.String("alg", "trial", "factoring algorithm to use")

var one = big.NewInt(1)
var nine = big.NewInt(9)
var ten = big.NewInt(10)

func main() {
	flag.Parse()

	rnd := rand.New(rand.NewSource(*seed))

	// Figure out the number to factor
	n := big.NewInt(0)
	args := flag.Args()
	if len(args) == 0 {
		fmt.Println("no number to factor")
		os.Exit(1)
	}
	if len(args) > 1 {
		fmt.Println("can't factor multiple numbers")
		os.Exit(1)
	}
	nstr := args[0]
	switch nstr[0] {
	case 'r':
		// random d-digit number
		d, err := strconv.Atoi(nstr[1:])
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
		k := big.NewInt(0).Exp(ten, big.NewInt(int64(d)-1), nil)
		n.Mul(k, nine)
		n.Rand(rnd, n)
		n.Add(n, k)
	case 's':
		// random d-digit semiprime
		d, err := strconv.Atoi(nstr[1:])
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
		digits := int64(d)
		if digits%2 != 0 {
			fmt.Println("semiprime must have an even number of digits")
			os.Exit(1)
		}
		min := big.NewInt(0).Exp(ten, big.NewInt(digits-1), nil)
		max := big.NewInt(0).Exp(ten, big.NewInt(digits), nil)
		for {
			x := randomPrime(digits / 2, rnd)
			y := randomPrime(digits / 2, rnd)
			n.Mul(x, y)
			if n.Cmp(min) >= 0 && n.Cmp(max) < 0 {
				break
			}
		}
	default:
		_, ok := n.SetString(nstr, 10)
		if !ok {
			fmt.Printf("parsing \"%s\": invalid number\n", nstr)
			os.Exit(1)
		}
		if n.Cmp(one) <= 0 {
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
func randomPrime(digits int64, rnd *rand.Rand) *big.Int {
	k := big.NewInt(0).Exp(ten, big.NewInt(digits-1), nil)
	x := big.NewInt(0)
	x.Mul(k, nine)
	n := big.NewInt(0)
	for {
		n.Rand(rnd, x)
		n.Add(n, k)
		if n.ProbablyPrime(1000) {
			return n
		}
	}
}