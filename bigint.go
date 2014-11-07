package factorlib

import (
	"math/big"
)

// A wrapper around math/big.Int which makes operations easier to express.
// The big difference is that a bigint is immutable.  It requires more allocations
// to make them immutable, but it's so much easier to write the code.

type bigint struct {
	v *big.Int
}

func NewBig(x int64) bigint {
	return bigint{big.NewInt(x)}
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

func (x bigint) Lsh(n uint) bigint {
	return bigint{new(big.Int).Lsh(x.v, n)}
}

func (x bigint) Rsh(n uint) bigint {
	return bigint{new(big.Int).Rsh(x.v, n)}
}

func (x bigint) SqrtFloor() bigint {
	return bigint{sqrtFloor(x.v)}
}

func (x bigint) Sign() int {
	return x.v.Sign()
}
func (x bigint) Cmp(y bigint) int {
	return x.v.Cmp(y.v)
}
func (x bigint) BitLen() int {
	return x.v.BitLen()
}
func (x bigint) String() string {
	return x.v.String()
}

func (x bigint) GetBig() *big.Int {
	return new(big.Int).Set(x.v)
}

func (x bigint) Int64() int64 {
	return x.v.Int64()
}
