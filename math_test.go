package factorlib

import (
	"math/big"
	"math/rand"
	"testing"
)

func TestLog(t *testing.T) {
	b := big.NewInt(0)
	for n := int64(0); n <= 65536; n++ {
		s := uint8(b.SetInt64(n).BitLen())
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

func TestSqrtFloor(t *testing.T) {
	b := big.NewInt(0)
	m := big.NewInt(0)
	for i := int64(0); i < 1000; i++ {
		b.SetInt64(i)
		c := sqrtFloor(b)
		m.Set(c)
		m.Mul(m, m)
		if m.Cmp(b) > 0 {
			t.Errorf("sqrtFloor(%d) = %d, too large", b, c)
		}
		m.Set(c)
		m.Add(c, one)
		m.Mul(m, m)
		if m.Cmp(b) <= 0 {
			t.Errorf("sqrtFloor(%d) = %d, too small", b, c)
		}
	}
}
func TestSqrtCeil(t *testing.T) {
	b := big.NewInt(0)
	m := big.NewInt(0)
	for i := int64(1); i < 1000; i++ {
		b.SetInt64(i)
		c := sqrtCeil(b)
		m.Set(c)
		m.Mul(m, m)
		if m.Cmp(b) < 0 {
			t.Errorf("sqrtCeil(%d) = %d, too small", b, c)
		}
		m.Set(c)
		m.Sub(c, one)
		m.Mul(m, m)
		if m.Cmp(b) >= 0 {
			t.Errorf("sqrtCeil(%d) = %d, too large", b, c)
		}
	}
	if sqrtCeil(zero).Sign() != 0 {
		t.Errorf("sqrtCeil(0) != 0")
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
