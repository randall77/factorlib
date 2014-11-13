package factorlib

import (
	"github.com/randall77/factorlib/big"
	"github.com/randall77/factorlib/primes"
	"math/rand"
)

func init() {
	factorizers["ecm"] = ecm
}

// we use the elliptic curve y^2 = x^3 + ax + 1 for a random a in Z_n
type point struct {
	x, y big.Int
}

func (p point) Inf() bool {
	// since x==y==0 is never a solution to the elliptic curve, we
	// reserve that bit pattern for encoding infinity.
	return p.x.IsZero() && p.y.IsZero()
}

func (p point) Check(n, a big.Int) bool {
	if p.Inf() {
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
	//log.Printf("add (%t %d %d) (%t %d %d)", p.inf, p.x, p.y, q.inf, q.x, q.y)
	if p.Inf() {
		return q
	}
	if q.Inf() {
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
	s := num.Mul(denom.Mod(n).ModInv(n)).Mod(n)
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
		d := a.Mul(a).Mul(a).Lsh(2).Add64(27).Mod(n)
		if d.IsZero() {
			// n divides 4a^3+27 - curve has repeating factors, so skip it.
			continue
		}

		p := point{big.Zero, big.One}
		for i := 0; ; i++ {
			p = p.Mul(primes.Get(i), n, a)
			if p.Inf() {
				// this curve didn't work
				break
			}
		}
	}
}
