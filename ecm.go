package factorlib

import (
	"github.com/randall77/factorlib/big"
	"github.com/randall77/factorlib/primes"
	"math/rand"
)

// see http://en.wikipedia.org/wiki/Lenstra_elliptic_curve_factorization

// The elliptic curve factorization method (ECM) uses points on the curve
//    y^2 = x^3 + ax + 1 mod n
// where a is chosen at random.
//
// When n is prime, the points on the curve plus a "point at infinity"
// form a group.  Mod a non-prime, they don't quite.  ECM tries to
// find one of these non-group cases and can extract a factorization
// of n from it.
//
// The group + operation on points on the curve is tricky.
// http://en.wikipedia.org/wiki/Elliptic_curve#The_group_law
// The point at infinity is the 0 for the group (x+0=0+x=x).
//  r = p + q:
//    if px != qx
//      s = (qy - py) / (qx - px)
//      rx = s^2 - px - qx
//      ry = s * (px - rx) - py
//    else
//      if py == -qy
//        r = 0,0
//      else
//        must be the case that p == q
//        s = (3*px^2 + a) / (2*py)
//        rx = s^2 - px - qx
//        ry = s * (px - rx) - py
// Where all computations are done mod n.

// When n is not prime, one of the divide steps might fail.  That
// happens because gcd(n, divisor) > 1, and that gives us a factor of
// n.

// We choose a at random, making sure that the elliptic curve is
// nonsingular by checking that 4a^3+27 mod n is not zero.

// Starting at a point p = (0, 1) on the curve, we compute kp for k's
// with lots of small factors.  If we reach the 0 element, then choose
// a new a and try again.  If the divide fails, we've found a factor
// of n.

func init() {
	factorizers["ecm"] = ecm
}

// we use the elliptic curve y^2 = x^3 + ax + 1 for a random a in Z_n
type point struct {
	x, y big.Int
}

func (p point) Zero() bool {
	// since x==y==0 is never a solution to the elliptic curve, we
	// reserve that bit pattern for encoding infinity.
	return p.x.IsZero() && p.y.IsZero()
}

func (p point) Check(n, a big.Int) bool {
	if p.Zero() {
		return true
	}
	lhs := p.y.Square()
	rhs := p.x.Square().Add(a).Mul(p.x).Add(big.One)
	return lhs.Sub(rhs).Mod(n).IsZero()
}

func (p point) Equals(q point) bool {
	return p.x.Equals(q.x) && p.y.Equals(q.y)
}

func (p point) Add(q point, n, a big.Int) point {
	if p.Zero() {
		return q
	}
	if q.Zero() {
		return p
	}
	var num, denom big.Int
	if !p.x.Equals(q.x) {
		num = p.y.Sub(q.y)
		denom = p.x.Sub(q.x)
	} else if p.y.Equals(q.y) && !p.y.IsZero() {
		// double point
		num = p.x.Square().Mul(big.Three).Add(a)
		denom = p.y.Lsh(1)
	} else {
		return point{big.Zero, big.Zero}
	}
	denom = denom.Mod(n)
	f := denom.GCD(n)
	if !f.Equals(big.One) {
		panic(f)
	}
	s := num.Mul(denom.ModInv(n)).Mod(n)
	rx := s.Square().Sub(p.x).Sub(q.x).Mod(n)
	ry := s.Mul(p.x.Sub(rx)).Sub(p.y).Mod(n)
	return point{rx, ry}
}

func (p point) Mul(k int64, n, a big.Int) point {
	// compute q=kp by repeated doubling
	q := point{big.Zero, big.Zero}
	for ; k > 1; k >>= 1 {
		if k&1 != 0 {
			q = q.Add(p, n, a)
		}
		p = p.Add(p, n, a)
	}
	return q.Add(p, n, a)
}

func ecm(n big.Int, rnd *rand.Rand) (r []big.Int) {
	// ecm does not work for powers of a single prime.  Check that first.
	if f := primepower(n, rnd); f != nil {
		return f
	}

	// somewhere in the guts of the ecm routine we panic a factor of n and
	// catch it here.
	defer func() {
		f := recover().(big.Int)
		r = []big.Int{f, n.Div(f)}
	}()
	for {
		a := n.Rand(rnd)
		if a.Cube().Lsh(2).Add64(27).Mod(n).IsZero() {
			// n divides 4a^3+27 - curve has repeating factors, so skip it.
			continue
		}

		p := point{big.Zero, big.One}
		for i := 0; ; i++ {
			p = p.Mul(primes.Get(i), n, a)
			if p.Zero() {
				// this curve didn't work
				break
			}
		}
	}
}
