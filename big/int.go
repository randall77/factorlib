package big

import (
	"fmt"
	"math"
	"math/big"
	"math/rand"
)

// A wrapper around math/big.Int which makes operations easier to express.
// The big difference is that a bigint is immutable.  It requires more allocations
// to make them immutable, but it's so much easier to write the code.

type Int struct {
	v *big.Int
}

func Int64(x int64) Int {
	return Int{big.NewInt(x)}
}

func (x Int) Add(y Int) Int {
	return Int{new(big.Int).Add(x.v, y.v)}
}

func (x Int) Sub(y Int) Int {
	return Int{new(big.Int).Sub(x.v, y.v)}
}

func (x Int) Mul(y Int) Int {
	return Int{new(big.Int).Mul(x.v, y.v)}
}

func (x Int) Div(y Int) Int {
	return Int{new(big.Int).Div(x.v, y.v)}
}
func (x Int) Mod(y Int) Int {
	return Int{new(big.Int).Mod(x.v, y.v)}
}

func (x Int) Add64(y int64) Int {
	z := big.NewInt(y)
	return Int{z.Add(x.v, z)}
}
func (x Int) Sub64(y int64) Int {
	z := big.NewInt(y)
	return Int{z.Sub(x.v, z)}
}
func (x Int) Mul64(y int64) Int {
	z := big.NewInt(y)
	return Int{z.Mul(x.v, z)}
}
func (x Int) Div64(y int64) Int {
	z := big.NewInt(y)
	return Int{z.Div(x.v, z)}
}
func (x Int) Mod64(y int64) int64 {
	z := big.NewInt(y)
	return z.Mod(x.v, z).Int64()
}

func (x Int) Square() Int {
	return x.Mul(x)
}
func (x Int) Neg() Int {
	return Int{new(big.Int).Neg(x.v)}
}
func (x Int) IsZero() bool {
	return x.v.Sign() == 0
}

func (x Int) Lsh(n uint) Int {
	return Int{new(big.Int).Lsh(x.v, n)}
}

func (x Int) Rsh(n uint) Int {
	return Int{new(big.Int).Rsh(x.v, n)}
}

func (x Int) Exp(k int64) Int {
	b := big.NewInt(k)
	b.Exp(x.v, b, nil)
	return Int{b}
}

func (x Int) ExpMod(e, m Int) Int {
	b := new(big.Int)
	b.Exp(x.v, e.v, m.v)
	return Int{b}
}

func (x Int) ModInv(n Int) Int {
	// TODO: check gcd(x,n)==1?
	return Int{new(big.Int).ModInverse(x.v, n.v)}
}

func (x Int) GCD(y Int) Int {
	return Int{new(big.Int).GCD(nil, nil, x.v, y.v)}
}

// return floor(sqrt(x)).
func (x Int) SqrtFloor() Int {
	if x.IsZero() {
		return x
	}
	b := uint(x.BitLen())

	// invariant lo <= sqrt(x) < hi
	lo := One.Lsh((b-1)/2)
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
func (x Int) SqrtCeil() Int {
	y := x.SqrtFloor()
	if y.Square().Cmp(x) != 0 {
		y = y.Add(One)
	}
	return y
}

func (x Int) Sign() int {
	return x.v.Sign()
}
func (x Int) Cmp(y Int) int {
	return x.v.Cmp(y.v)
}
var minInt64 = Int64(math.MinInt64)
func (x Int) Cmp64(y int64) int {
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
func (x Int) BitLen() int {
	return x.v.BitLen()
}
func (x Int) Bit(i int) uint {
	return x.v.Bit(i)
}
func (x Int) ProbablyPrime(n int) bool {
	return x.v.ProbablyPrime(n)
}
func (x Int) Format(s fmt.State, ch rune) {
	x.v.Format(s, ch)
}

func (x Int) Int64() int64 {
	return x.v.Int64()
}

// Rand returns a random number in [0,x)
func (x Int) Rand(rnd *rand.Rand) Int {
	return Int{new(big.Int).Rand(rnd, x.v)}
}

func ParseInt(s string) (Int, bool) {
	y := new(big.Int)
	_, ok := y.SetString(s, 10)
	return Int{y}, ok
}

// helpful constants
var Zero = Int64(0)
var One = Int64(1)
var Two = Int64(2)
var Ten = Int64(10)
