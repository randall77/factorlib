package factorlib

import (
	"math/big"
	"math/rand"
)

// generic math routines
//
// There are often two versions of each routine, a int64 one and a
// big.Int one.  The int64 ones should only be used if arguments are
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
		panic("uninvertible")
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
	var a, d int64
	for {
		a = 1 + rnd.Int63n(p-1)
		d = (a*a - n) % p
		if d < 0 {
			d += p
		}
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

// returns true iff x^2 == a mod p has a solution for x.
// requires: 0 <= a < p, p prime
func quadraticResidueBig(a, p *big.Int) bool {
	if a.Cmp(two) < 0 {
		return true
	}
	// a is a quadratic residue (x^2 == a has a solution) iff
	// a^((p-1)/2) == 1 mod p.
	e := big.NewInt(0).Rsh(p, 1)
	return e.Exp(a, e, p).Cmp(one) == 0
}

// solves x^2 == a mod p
// returns nil if there is no solution.
func sqrtModPBig(n, p *big.Int, rnd *rand.Rand) *big.Int {
	if n.Cmp(two) < 0 {
		return n
	}
	if p.Bit(1) == 1 { // p == 3 mod 4
		// result = n^((p+1)/4) mod p
		x := big.NewInt(0).Add(p, one)
		x.Rsh(x, 2)
		x.Exp(n, x, p)
		return x
	}
	// p == 1 mod 4
	// Cipolla's algorithm (http://en.wikipedia.org/wiki/Cipolla%27s_algorithm)

	// pick a quadratic nonresidue for a
	a := big.NewInt(0)
	for {
		a.Rand(rnd, p)
		if !quadraticResidueBig(a, p) {
			break
		}
	}

	// compute a^2-n
	d := big.NewInt(0).Exp(a, two, p)
	d.Sub(d, n)
	if d.Sign() < 0 {
		d.Add(d, p)
	}

	// number in F_p^2, x0 + x1 * w, where w = sqrt(a^2-n)
	// compute (x0 + x1*w)^((p-1)/2)
	x0 := a
	x1 := big.NewInt(1)
	r0 := big.NewInt(1)
	r1 := big.NewInt(0)

	m1 := big.NewInt(0)
	m2 := big.NewInt(0)
	m3 := big.NewInt(0)
	m4 := big.NewInt(0)

	// repeated sqaring of x, multiplying into r
	q := big.NewInt(0).Add(p, one)
	q.Rsh(q, 1)
	for q.Sign() != 0 {
		if q.Bit(0) != 0 {
			// r *= x
			m1.Mul(r0, x0)
			m2.Mul(r0, x1)
			m3.Mul(r1, x0)
			m4.Mul(r1, x1)
			r0.Mul(m4, d)
			r0.Add(r0, m1)
			r1.Add(m2, m3)
			r0.Mod(r0, p)
			r1.Mod(r1, p)
		}
		// x *= x
		m1.Mul(x0, x0)
		m2.Mul(x0, x1)
		m3.Mul(x1, x0)
		m4.Mul(x1, x1)
		x0.Mul(m4, d)
		x0.Add(x0, m1)
		x1.Add(m2, m3)
		x0.Mod(x0, p)
		x1.Mod(x1, p)

		q.Rsh(q, 1)
	}
	if r1.Sign() != 0 {
		panic("still in Z_p^2 - Z_p")
	}
	// check our answer
	z := big.NewInt(0).Mul(r0, r0)
	z.Mod(z, p)
	if z.Cmp(n) != 0 {
		panic("bad sqrtMod")
	}
	return r0
}

// return floor(sqrt(x)).
func sqrtFloor(x *big.Int) *big.Int {
	if x.Sign() == 0 {
		return x
	}
	b := uint(x.BitLen())
	
	// invariant lo <= sqrt(x) < hi
	lo := big.NewInt(1)
	lo.Lsh(lo, (b-1)/2)
	hi := big.NewInt(1)
	hi.Lsh(hi, (b+1)/2)
	m := big.NewInt(0)
	m2 := big.NewInt(0)
	for {
		m.Add(lo, hi).Rsh(m, 1)
		if m.Cmp(lo) == 0 {
			return lo
		}
		m2.Mul(m, m)
		if m2.Cmp(x) <= 0 {
			lo, m = m, lo
		} else {
			hi, m = m, hi
		}
	}
}
func sqrtCeil(x *big.Int) *big.Int {
	y := sqrtFloor(x)
	y2 := big.NewInt(0)
	y2.Mul(y, y)
	if y2.Cmp(x) != 0 {
		y = y2.Add(y, one)
	}
	return y
}
