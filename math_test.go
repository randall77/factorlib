package factorlib

import (
	"math/rand"
	"testing"
)

func TestLog(t *testing.T) {
	for n := int64(0); n <= 65536; n++ {
		s := uint8(Big(n).BitLen())
		r := log2(n)
		if r != s {
			t.Errorf("n=%d want %d, got %d", n, s, r)
		}
	}
}

func TestGCD(t *testing.T) {
	const n = 200
	for x := int64(0); x < n; x++ {
		for y := int64(0); y < n; y++ {
			if x == 0 && y == 0 {
				// Result is undefined - don't bother testing
				continue
			}
			var g int64
			for z := int64(1); z < n; z++ {
				if x%z == 0 && y%z == 0 {
					g = z
				}
			}
			if gcd(x, y) != g {
				t.Errorf("gcd(%d,%d)=%d, want %d", x, y, gcd(x, y), g)
			}
		}
	}
}

func TestModInv(t *testing.T) {
	for n := int64(2); n < 1000; n++ {
		for x := int64(1); x < n; x++ {
			if gcd(n, x) != 1 {
				continue
			}
			y := modInv(x, n)
			if x*y%n != 1 {
				t.Errorf("n=%d x=%d y=%d xy=%d", n, x, y, x*y%n)
			}
		}
	}
}

func TestQR(t *testing.T) {
	for i := 0; i < 1000; i++ {
		p := getPrime(i)

		// find all quadratic residues
		m := map[int64]struct{}{}
		for a := int64(0); a < p; a++ {
			m[a*a%p] = struct{}{}
		}

		for a := int64(0); a < p; a++ {
			_, isQR := m[a]
			r := quadraticResidue(a, p)
			if r != isQR {
				t.Errorf("p=%d a=%d want %t, got %t", p, a, isQR, r)
			}
		}
	}
}

func TestSqrtModP(t *testing.T) {
	rnd := rand.New(rand.NewSource(123))
	for i := 0; i < 1000; i++ {
		p := getPrime(i)

		// compute roots mod p
		m := map[int64]int64{}
		for a := int64(0); a < p; a++ {
			m[a*a%p] = a
		}

		for a := int64(0); a < p; a++ {
			s, ok := m[a]
			if !ok {
				// a is a quadratic nonresidue
				continue
			}
			r := sqrtModP(a, p, rnd)
			if r != s && r != p-s {
				t.Errorf("p=%d a=%d want %d or %d, got %d", p, a, s, p-s, r)
			}
		}
	}
}

func TestSqrtModPK(t *testing.T) {
	rnd := rand.New(rand.NewSource(123))
	for i := 0; i < 1000; i++ {
		p := getPrime(i)
		for k := uint(1); ; k++ {
			pk := exp(p, k)
			if pk > 10000 {
				break
			}

			// compute roots mod p^k
			m := map[int64][]int64{}
			for a := int64(0); a < pk; a++ {
				m[a*a%pk] = append(m[a*a%pk], a)
			}

			for a := int64(0); a < pk; a++ {
				if a != 0 && gcd(a, pk) != 1 {
					// a is not relatively prime to p^k
					continue
				}
				s, ok := m[a]
				if !ok {
					// a is a quadratic nonresidue
					continue
				}
				r := sqrtModPK(a, p, k, rnd)
				ok = false
				for _, x := range s {
					if x == r {
						ok = true
						break
					}
				}
				if !ok {
					t.Errorf("pk=%d a=%d want element of %#v, got %d", pk, a, s, r)
				}
			}
		}
	}
}

func TestSqrtModN(t *testing.T) {
	rnd := rand.New(rand.NewSource(123))
	// test square roots mod 5^i 7^j 11^k for all quadratic residues mod those numbers.
	var primepowers = [3]primePower{{5, 0}, {7, 0}, {11, 0}}
	for i := uint(0); i < 4; i++ {
		primepowers[0].k = i
		for j := uint(0); j < 4; j++ {
			primepowers[1].k = j
			for k := uint(0); k < 3; k++ {
				primepowers[2].k = k

				n := int64(1)
				for _, pp := range primepowers {
					for z := uint(0); z < pp.k; z++ {
						n *= pp.p
					}
				}
				for a := int64(0); a < n; a++ {
					if gcd(a, n) != 1 {
						continue
					}
					if !quadraticResidue(a%5, 5) {
						continue
					}
					if !quadraticResidue(a%7, 7) {
						continue
					}
					if !quadraticResidue(a%11, 11) {
						continue
					}
					x := sqrtModN(a, primepowers[:], rnd)
					if x*x%n != a {
						t.Errorf("bad sqrtModN a=%d n=%d x=%d\n", a, n, x)
					}
				}
			}
		}
	}
}

func TestQuadraticModP(t *testing.T) {
	rnd := rand.New(rand.NewSource(123))
	for i := 0; i < 25; i++ {
		p := getPrime(i)
		for a := int64(1); a < p; a++ {
			for b := int64(0); b < p; b++ {
				for c := int64(0); c < p; c++ {
					s := quadraticModP(a, b, c, p, rnd)
					if len(s) > 2 {
						t.Errorf("too many quadratic solutions")
					}
					if len(s) == 2 && s[0] == s[1] {
						t.Errorf("returned same root twice")
					}
					for _, x := range s {
						if (a*x*x+b*x+c)%p != 0 {
							t.Errorf("p=%d a=%d b=%d c=%d x=%d (ax^2+bx+c)%%p=%d", p, a, b, c, x, ((a*x*x+b*x+c)%p+p)%p)
						}
					}
					cnt := 0
					for x := int64(0); x < p; x++ {
						if (a*x*x+b*x+c)%p == 0 {
							cnt++
						}
					}
					if cnt != len(s) {
						t.Errorf("expected %d results, got %d\n", cnt, len(s))
					}
				}
			}
		}
	}
}
