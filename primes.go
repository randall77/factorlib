package factorlib

import (
	"fmt"
)

// storage for all the primes we've found so far.
var primes = []int64{2}

// getPrime(0) == 2, getPrime(1) == 3, ...
func getPrime(i int) int64 {
	genPrimes(i)
	return primes[i]
}

// pi(x) = number of primes <= x.  Used for double-checking our work.
var pi = []struct {
	x   int64
	pix int64
}{
	{10, 4},
	{100, 25},
	{1000, 168},
	{10000, 1229},
	{100000, 9592},
	{1000000, 78498},
	{10000000, 664579},
	{100000000, 5761455},
	{1000000000, 50847534},
	{10000000000, 455052511},
	{100000000000, 4118054813},
	{1000000000000, 37607912018},
	{10000000000000, 346065536839},
}

func addPrime(p int64) {
	// check against pi(x) table
	if len(pi) > 0 && p > pi[0].x {
		if int64(len(primes)) != pi[0].pix {
			panic(fmt.Sprintf("pi(%d)=%d, we got %d", pi[0].x, pi[0].pix, len(primes)))
		}
		pi = pi[1:]
	}
	primes = append(primes, p)
}

// must be a multiple of 30
const sievewidth = 30 * 20000

var sieved int64

// generate up to the ith prime.
// TODO: make this multithreaded-safe
func genPrimes(i int) {
again:
	if i < len(primes) {
		return
	}
	if sieved == 0 {
		// simple sieve for primes up to 60000
		// index i represents number 2*i+3
		const n = (60000 - 3 + 1) / 2
		if n*n <= sievewidth {
			panic("need to raise initial sieve size")
		}
		s := make([]bool, n)
		for i := int64(0); i < 128; i++ { // 2*128+3 > sqrt(60000)
			if s[i] {
				continue
			}
			addPrime(2*i + 3)
			// mark as composite values starting at p^2, stepping by 2p
			for j := 2*i*(i+3) + 3; j < n; j += 2*i + 3 {
				s[j] = true
			}
		}
		for i := int64(128); i < n; i++ {
			if !s[i] {
				addPrime(2*i + 3)
			}
		}
		sieved = 60000
		goto again
	}

	// Sieve up some more primes.
	// Each byte in the sieve covers 30 integers.  Each bit in the byte represents
	// one of the eight numbers in those 30 that are nonzero mod 2,3, and 5.
	start := sieved
	end := start + sievewidth

	// initialize sieve
	for i := 0; i < sievewidth/30; i++ {
		sieve[i] = 0
	}

	// mark all composite numbers in the sieve
	for i := 3; ; i++ { // "3" means start at p=7
		p := primes[i]
		if p*p >= end {
			// We don't need to sieve with primes >= sqrt(end).
			// If they have such a factor, they will also have a
			// factor < sqrt(end).
			break
		}

		// compute starting offset.  It is the offset from sieve start
		// of the first multiple of p which is >= p^2.
		x := p*p - start
		if x < 0 {
			x += (-x + p - 1) / p * p
		}
		for j := 0; j < 30; j++ {
			// mark all integers starting at x+j*p with stride 30*p.
			y := x + int64(j)*p
			b := deltaIdx[y%30]
			if b == 255 {
				// all values are multiples of 2,3, or 5
				continue
			}
			block := y / 30
			mask := byte(1) << b
			// This is the inner loop.  We sweep through the sieve
			// 8 times for each prime < sqrt(end).
			for block < sievewidth/30 {
				sieve[block] |= mask
				block += p
			}
		}
	}
	// pick primes out of the sieve
	for i := 0; i < sievewidth/30; i++ {
		s := sieve[i]
		for j := 0; j < 8; j++ {
			if s&1 == 0 {
				addPrime(start + int64(i)*30 + int64(deltas[j]))
			}
			s >>= 1
		}
	}
	sieved = end
	goto again
}

var sieve [sievewidth / 30]byte

// these are the 8 indexes mod 30 which are not 0 mod 2,3,5
var deltas = [8]byte{1, 7, 11, 13, 17, 19, 23, 29}
var deltaIdx [30]byte

func init() {
	for i := 0; i < 30; i++ {
		deltaIdx[i] = 255
	}
	for i := 0; i < 8; i++ {
		deltaIdx[deltas[i]] = byte(i)
	}
}
