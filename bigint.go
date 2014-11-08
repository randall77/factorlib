package factorlib

import (
	"fmt"
	"math"
	"math/big"
)

// A wrapper around math/big.Int which makes operations easier to express.
// The big difference is that a bigint is immutable.  It requires more allocations
// to make them immutable, but it's so much easier to write the code.

type bigint struct {
	v *big.Int
}

func Big(x int64) bigint {
	return bigint{big.NewInt(x)}
}
func BigFromBig(x *big.Int) bigint {
	return bigint{new(big.Int).Set(x)}
}

func (x bigint) Add(y bigint) bigint {
	return bigint{new(big.Int).Add(x.v, y.v)}
}

func (x bigint) Sub(y bigint) bigint {
	return bigint{new(big.Int).Sub(x.v, y.v)}
}

func (x bigint) Mul(y bigint) bigint {
	return bigint{new(big.Int).Mul(x.v, y.v)}
}

func (x bigint) Div(y bigint) bigint {
	return bigint{new(big.Int).Div(x.v, y.v)}
}
func (x bigint) Mod(y bigint) bigint {
	return bigint{new(big.Int).Mod(x.v, y.v)}
}

func (x bigint) Add64(y int64) bigint {
	z := big.NewInt(y)
	return bigint{z.Add(x.v, z)}
}
func (x bigint) Sub64(y int64) bigint {
	z := big.NewInt(y)
	return bigint{z.Sub(x.v, z)}
}
func (x bigint) Mul64(y int64) bigint {
	z := big.NewInt(y)
	return bigint{z.Mul(x.v, z)}
}
func (x bigint) Div64(y int64) bigint {
	z := big.NewInt(y)
	return bigint{z.Div(x.v, z)}
}
func (x bigint) Mod64(y int64) int64 {
	z := big.NewInt(y)
	return z.Mod(x.v, z).Int64()
}

func (x bigint) Square() bigint {
	return x.Mul(x)
}
func (x bigint) Neg() bigint {
	return bigint{new(big.Int).Neg(x.v)}
}
func (x bigint) IsZero() bool {
	return x.v.Sign() == 0
}

func (x bigint) Lsh(n uint) bigint {
	return bigint{new(big.Int).Lsh(x.v, n)}
}

func (x bigint) Rsh(n uint) bigint {
	return bigint{new(big.Int).Rsh(x.v, n)}
}

func (x bigint) Exp(k int64) bigint {
	b := big.NewInt(k)
	b.Exp(x.v, b, nil)
	return bigint{b}
}

func (x bigint) ExpMod(e, m bigint) bigint {
	b := new(big.Int)
	b.Exp(x.v, e.v, m.v)
	return bigint{b}
}

func (x bigint) ModInv(n bigint) bigint {
	// TODO: check gcd(x,n)==1?
	return bigint{new(big.Int).ModInverse(x.v, n.v)}
}

func (x bigint) GCD(y bigint) bigint {
	return bigint{new(big.Int).GCD(nil, nil, x.v, y.v)}
}

// return floor(sqrt(x)).
func (x bigint) SqrtFloor() bigint {
	if x.IsZero() {
		return x
	}
	b := uint(x.BitLen())

	// invariant lo <= sqrt(x) < hi
	lo := one.Lsh((b-1)/2)
	hi := lo.Lsh(1)
	for {
		m := lo.Add(hi).Rsh(1)
		if m.Cmp(lo) == 0 {
			return lo
		}
		if m.Square().Cmp(x) <= 0 {
			lo, m = m, lo
		} else {
			hi, m = m, hi
		}
	}
}
func (x bigint) SqrtCeil() bigint {
	y := x.SqrtFloor()
	if y.Square().Cmp(x) != 0 {
		y = y.Add(one)
	}
	return y
}

func (x bigint) Sign() int {
	return x.v.Sign()
}
func (x bigint) Cmp(y bigint) int {
	return x.v.Cmp(y.v)
}
var minInt64 = Big(math.MinInt64)
func (x bigint) Cmp64(y int64) int {
	if x.BitLen() >= 64 {
		if x.Sign() > 0 {
			return 1
		}
		if x.Cmp(minInt64) == 0 && y == math.MinInt64 {
			return 0
		}
		return -1
	}
	z := x.Int64()
	if z > y {
		return 1
	}
	if z < y {
		return -1
	}
	return 0
}
func (x bigint) BitLen() int {
	return x.v.BitLen()
}
func (x bigint) Bit(i int) uint {
	return x.v.Bit(i)
}
func (x bigint) ProbablyPrime(n int) bool {
	return x.v.ProbablyPrime(n)
}
func (x bigint) Format(s fmt.State, ch rune) {
	x.v.Format(s, ch)
}

func (x bigint) GetBig() *big.Int {
	return new(big.Int).Set(x.v)
}

func (x bigint) Int64() int64 {
	return x.v.Int64()
}

// helpful constants
var zero = Big(0)
var one = Big(1)
var two = Big(2)
