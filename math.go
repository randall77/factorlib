package factorlib

import (
	"fmt"
	"github.com/randall77/factorlib/big"
	"math/rand"
)

// generic math routines
//
// There are often two versions of each routine, a int64 one and a
// bigint one.  The int64 ones should only be used if arguments are
// known to be < 2^31. (That leaves room for a*x+b*y computations
// without overflowing.)

// return ceil(log_2(n))
func log2(n int64) uint8 {
	var r uint8
	for n != 0 {
		n >>= 1
		r++
	}
	return r
}

// exp returns x^e
func exp(x int64, e uint) int64 {
	r := int64(1)
	for e != 0 {
		if e&1 != 0 {
			r = r * x
		}
		x = x * x
		e >>= 1
	}
	return r
}

// expMod returns x^e mod p
func expMod(x, e, p int64) int64 {
	r := int64(1)
	for e != 0 {
		if e&1 != 0 {
			r = r * x % p
		}
		x = x * x % p
		e >>= 1
	}
	return r
}

// gcd returns the greatest common divisor of x and y.
//   x >= 0 && y >= 0
//   x != 0 || y != 0
func gcd(x, y int64) int64 {
	if x < y {
		x, y = y, x
	}
	for { // x >= y
		if y == 0 {
			return x
		}
		x, y = y, x%y
	}
}

// modInv returns y such that x*y % n == 1.
//   gcd(x,n) == 1
func modInv(x, n int64) int64 {
	// http://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Computing_multiplicative_inverses_in_modular_structures
	t := int64(0)
	newt := int64(1)
	r := n
	newr := x
	for newr != 0 {
		q := r / newr
		t, newt = newt, t-q*newt
		r, newr = newr, r-q*newr
	}
	if r > 1 {
		panic(fmt.Sprintf("%d is uninvertible mod %d", x, n))
	}
	if t < 0 {
		t += n
	}
	return t
}

// returns true iff there exists an x such that x^2 == n mod p.
//   p must be prime
//   0 <= n < p
func quadraticResidue(n, p int64) bool {
	if n < 2 {
		return true
	}
	// a is a quadratic residue (x^2 == a has a solution) iff
	// a^((p-1)/2) == 1 mod p.
	return expMod(n, p>>1, p) == 1
}

// sqrtModP finds an x such that x^2 == n mod p.
//   p must be prime
//   0 <= n < p
//   quadraticResidue(n, p) must be true
// The algorithm is randomized, uses rnd for random bits.
func sqrtModP(n, p int64, rnd *rand.Rand) int64 {
	if n < 2 {
		return n
	}
	if p%4 == 3 {
		return expMod(n, (p+1)>>2, p)
	}
	// Cipolla's algorithm (http://en.wikipedia.org/wiki/Cipolla's_algorithm)
	var a, d int64
	for {
		a = 1 + rnd.Int63n(p-1)
		d = (a*a + p - n) % p
		if !quadraticResidue(d, p) {
			break
		}
	}

	x0 := a
	x1 := int64(1)
	r0 := int64(1)
	r1 := int64(0)
	q := (p + 1) >> 1
	for q != 0 {
		if q&1 != 0 {
			// r *= x
			r0, r1 = (r0*x0+r1*x1%p*d)%p, (r0*x1+r1*x0)%p
		}
		// x *= x
		x0, x1 = (x0*x0+x1*x1%p*d)%p, (2*x0*x1)%p
		q >>= 1
	}
	if r1 != 0 {
		panic("bad f_p^2 element")
	}
	return r0
}

// sqrtModPK finds an x such that x^2 == n mod p^k.
//   p must be prime
//   0 <= n < p^k
//   if n != 0, p does not divide n
//   k >= 1
//   quadraticResidue(n, p) must be true
// The algorithm is randomized, uses rnd for random bits.
func sqrtModPK(n, p int64, k uint, rnd *rand.Rand) int64 {
	if n < 2 {
		return n
	}
	if p == 2 {
		// p == 2 is weird.  We'll solve one bit at a time.
		L := []int64{0}
		mask := int64(1)
		for b := uint(0); b < k; b++ {
			var L2 []int64
			for _, v := range L {
				if v*v&mask == n&mask {
					L2 = append(L2, v)
				}
				w := v + int64(1)<<b
				if w*w&mask == n&mask {
					L2 = append(L2, w)
				}
			}
			mask = 2*mask + 1
			L = L2
		}
		return L[0]
	}

	// first solve x^2 == n mod p
	r := sqrtModP(n%p, p, rnd)

	pi := p
	for i := uint(1); i < k; i++ {
		// r is a root of x^2 - n mod p^i.  Find a root of x^2 - n mod p^(i+1)
		// use Hensel's lemma: http://en.wikipedia.org/wiki/Hensel's_lemma
		// f(x) = x^2 - n
		// f'(x) = 2x != 0 mod p^k for p>2
		// t = (n-r^2)/p^i * (2r)^-1 mod p
		// s = r + t * p^i
		// TODO: lift by doubling i instead of incrementing i
		t := (n + (pi*p-r)*r) / pi % p
		t = t * modInv(2*r, p) % p
		r += t * pi
		pi *= p
	}
	return r
}

type primePower struct {
	p int64
	k uint
}

// returns a solution to x^2 == a mod n, where n is a product of the listed prime powers.
//   gcd(a,n) == 1
//   n[i].k >= 1
//   a must be a quadratic residue mod each prime
func sqrtModN(a int64, n []primePower, rnd *rand.Rand) int64 {
	if a <= 1 {
		return a
	}

	// compute N = product of all primes
	N := int64(1)
	for _, pp := range n {
		for i := uint(0); i < pp.k; i++ {
			N *= pp.p
		}
	}

	// use Chinese Remainder Theorem to compute result one p^k at a time.
	r := int64(0)
	for _, pp := range n {
		// compute p^k
		pk := int64(1)
		for i := uint(0); i < pp.k; i++ {
			pk *= pp.p
		}

		// find a solution to x^2 == a mod p^k
		x := sqrtModPK(a%pk, pp.p, pp.k, rnd)

		// add it in to total result
		M := N / pk
		r += M * x * modInv(M%pk, pk) // TODO: check for overflow
		r %= N
	}
	// check result
	if r*r%N != a {
		panic("bad sqrt")
	}
	return r
}

// returns a solution to x^2 == a mod n, where n is a product of the listed prime powers.
//   gcd(a,n) == 1
//   n[i].k >= 1
//   a must be a quadratic residue mod each prime
func bigSqrtModN(a big.Int, n []primePower, rnd *rand.Rand) big.Int {
	if a.Cmp(big.One) <= 0 {
		return a
	}

	// compute N = product of all primes
	N := big.Int64(1)
	for _, pp := range n {
		for i := uint(0); i < pp.k; i++ {
			N = N.Mul64(pp.p)
		}
	}

	// use Chinese Remainder Theorem to compute result one p^k at a time.
	r := big.Int64(0)
	for _, pp := range n {
		// compute p^k
		pk := int64(1)
		for i := uint(0); i < pp.k; i++ {
			pk *= pp.p
		}

		// find a solution to x^2 == a mod p^k
		x := sqrtModPK(a.Mod64(pk), pp.p, pp.k, rnd)

		// add it in to total result
		M := N.Div64(pk)
		r = r.Add(M.Mul64(x).Mul64(modInv(M.Mod64(pk), pk))).Mod(N)
	}
	// check result
	if !r.Square().Mod(N).Equals(a) {
		panic("bad sqrt")
	}
	return r
}

// solve ax^2+bx+c==0 mod p
//   0 <= a,b,c < p
//   a != 0
func quadraticModP(a, b, c, p int64, rnd *rand.Rand) []int64 {
	if p == 2 {
		// special case, easy to handle.
		// (2 is not a unit mod 2, so 1/2a doesn't work when p==2)
		if a^b == 0 {
			if c == 0 {
				return []int64{0, 1}
			}
			return nil
		}
		if c == 0 {
			return []int64{0}
		}
		return []int64{1}
	}

	d := (b*b + 4*(p-a)*c) % p
	if !quadraticResidue(d, p) {
		return nil
	}
	d = sqrtModP(d, p, rnd)
	i := modInv(2*a%p, p)
	r := []int64{(p - b + d) * i % p}
	if d != 0 {
		r = append(r, (2*p-b-d)*i%p)
	}
	return r
}

// solve ax^2+bx+c==0 mod p^k
//   0 <= a,b,c < p^k
//   gcd(a,p) == 1
//   pk == p^k
func quadraticModPK(a, b, c, p int64, k uint, pk int64, rnd *rand.Rand) []int64 {
	if p == 2 {
		if k > 1 {
			// TODO: k > 1
			return nil
		}
		// special case, easy to handle.
		// (2 is not a unit mod 2, so 1/2a doesn't work when p==2)
		if a^b == 0 {
			if c == 0 {
				return []int64{0, 1}
			}
			return nil
		}
		if c == 0 {
			return []int64{0}
		}
		return []int64{1}
	}

	d := (b*b + 4*(pk-a)*c) % pk
	e := d % p
	if e == 0 || !quadraticResidue(e, p) {
		// TODO: there might be solutions if e == 0.  Figure that out.
		return nil
	}
	d = sqrtModPK(d, p, k, rnd)
	i := modInv(2*a%pk, pk)
	r := []int64{(pk - b + d) * i % pk}
	if d != 0 {
		r = append(r, (2*pk-b-d)*i%pk)
	}
	return r
}
