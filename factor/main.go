package main

import (
	"flag"
	"fmt"
	"log"
	"math/rand"
	"os"
	"runtime/pprof"
	"strconv"

	"github.com/randall77/factorlib"
	"github.com/randall77/factorlib/big"
)

var seed = flag.Int64("seed", 0, "seed for RNG")
var alg = flag.String("alg", "trial", "factoring algorithm to use")
var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")

func main() {
	flag.Parse()

	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	rnd := rand.New(rand.NewSource(*seed))

	// Figure out the number to factor
	args := flag.Args()
	if len(args) == 0 {
		log.Fatal("no number to factor")
	}
	if len(args) > 1 {
		log.Fatal("can't factor multiple numbers")
	}
	var n big.Int
	nstr := args[0]
	switch nstr[0] {
	case 'r':
		// random d-digit number
		d, err := strconv.Atoi(nstr[1:])
		if err != nil {
			log.Fatal(err)
		}
		k := big.Ten.Exp(int64(d) - 1)
		n = k.Mul64(9).Rand(rnd).Add(k)
	case 's':
		// random d-digit semiprime
		d, err := strconv.Atoi(nstr[1:])
		if err != nil {
			log.Fatal(err)
		}
		if d%2 != 0 {
			log.Fatal("semiprime must have an even number of digits")
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
			log.Fatalf("parsing \"%s\": invalid number\n", nstr)
		}
		if n.Cmp(big.One) <= 0 {
			log.Fatalf("invalid n: %s\n", nstr)
		}
	}

	// Call into main library to do factoring
	var algstr string
	if *alg != "" {
		algstr = fmt.Sprintf(" using algorithm %s", *alg)
	}
	log.Printf("factoring %d%s\n", n, algstr)
	factors, err := factorlib.Factor(n, *alg, rnd)

	// Print result
	if err != nil {
		fmt.Printf("factorization failed: %s\n", err)
		return
	}
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
