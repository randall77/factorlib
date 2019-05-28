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
var logfile = flag.String("logfile", "", "logging file (default=stderr)")

func main() {
	flag.Parse()

	// Set up output streams and logger.
	var logger *log.Logger
	var fatal func(msg string, args ...interface{})
	var output func(msg string, args ...interface{})
	if *logfile == "" {
		logger = log.New(os.Stderr, "", log.LstdFlags)
		fatal = func(msg string, args ...interface{}) {
			fmt.Fprintf(os.Stderr, msg+"\n", args...)
			os.Exit(1)
		}
		output = func(msg string, args ...interface{}) {
			fmt.Fprintf(os.Stdout, msg+"\n", args...)
		}
	} else {
		w, err := os.Create(*logfile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "can't open log file %s\n", *logfile)
			os.Exit(1)
		}
		defer w.Close()
		logger = log.New(w, "", log.LstdFlags)
		fatal = func(msg string, args ...interface{}) {
			logger.Printf(msg, args...)
			fmt.Fprintf(os.Stderr, msg+"\n", args...)
			w.Close()
			os.Exit(1)
		}
		output = func(msg string, args ...interface{}) {
			logger.Printf(msg, args...)
			fmt.Fprintf(os.Stdout, msg+"\n", args...)
		}
	}

	// Set up profile output, if requested.
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			fatal("can't create profile file: %v", err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
		logger.Printf("cpu profile output to %s", *cpuprofile)
	}

	// Initialize random seed.
	rnd := rand.New(rand.NewSource(*seed))
	logger.Printf("seed: %d", *seed)

	// Figure out the number to factor
	args := flag.Args()
	logger.Printf("args: %s", args)
	if len(args) == 0 {
		fatal("no number to factor")
	}
	if len(args) > 1 {
		fatal("can't factor multiple numbers")
	}
	var n big.Int
	nstr := args[0]
	switch nstr[0] {
	case 'r':
		// random d-digit number
		d, err := strconv.Atoi(nstr[1:])
		if err != nil {
			fatal("%v", err)
		}
		k := big.Ten.Exp(int64(d) - 1)
		n = k.Mul64(9).Rand(rnd).Add(k)
	case 's':
		// random d-digit semiprime
		d, err := strconv.Atoi(nstr[1:])
		if err != nil {
			fatal("%v", err)
		}
		if d%2 != 0 {
			fatal("semiprime %s must request an even number of digits", nstr)
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
			fatal("parsing \"%s\": invalid number", nstr)
		}
	}
	logger.Printf("factoring %d", n)

	// Call into main library to do factoring
	if *alg != "" {
		logger.Printf("using algorithm %s", *alg)
	}
	factors, err := factorlib.Factor(n, *alg, rnd)

	// Print result
	if err != nil {
		output("factorization failed: %s", err)
		return
	}
	s := fmt.Sprintf("%d = ", n)
	for i, f := range factors {
		if i > 0 {
			s += "·"
		}
		s += fmt.Sprintf("%d", f)
	}
	output(s)
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
