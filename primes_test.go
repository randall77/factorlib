package factorlib

import (
	"bufio"
	"os"
	"strconv"
	"testing"
)

func TestPrimes(t *testing.T) {
	f, err := os.Open("primes.txt")
	if err != nil {
		t.Fatal(err)
	}
	scanner := bufio.NewScanner(f)
	scanner.Split(bufio.ScanWords)
	var a []int64
	for scanner.Scan() {
		p, err := strconv.Atoi(scanner.Text())
		if err != nil {
			t.Error(err)
		}
		a = append(a, int64(p))
	}
	if len(a) != 10000 {
		t.Errorf("want 10000 primes, got %d", len(a))
	}
	for i, p := range a {
		if getPrime(i) != p {
			t.Errorf("prime %d: want %d, got %d", i, p, getPrime(i))
		}
	}
}
