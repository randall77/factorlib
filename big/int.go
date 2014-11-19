package big

import (
	"fmt"
	"math"
	"math/big"
	"math/rand"
)

// A wrapper around math/big.Int which makes operations easier to express.
// The big difference is that this big.Int is immutable.  Operations on
// these big.Ints are easier to write code with, but require more
// allocations under the hood.  Totally worth it.

type Int struct {
	v *big.Int
}

// Constructors

func Int64(x int64) Int {
	return Int{big.NewInt(x)}
}

func ParseInt(s string) (Int, bool) {
	y, ok := new(big.Int).SetString(s, 10)
	return Int{y}, ok
}

// Arithmetic

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

func (x Int) Neg() Int {
	return Int{new(big.Int).Neg(x.v)}
}

func (x Int) Lsh(n uint) Int {
	return Int{new(big.Int).Lsh(x.v, n)}
}

func (x Int) Rsh(n uint) Int {
	return Int{new(big.Int).Rsh(x.v, n)}
}

// Info extraction

func (x Int) Int64() int64 {
	return x.v.Int64()
}

func (x Int) IsZero() bool {
	return x.v.Sign() == 0
}

func (x Int) Sign() int {
	return x.v.Sign()
}

func (x Int) Cmp(y Int) int {
	return x.v.Cmp(y.v)
}

func (x Int) Equals(y Int) bool {
	return x.v.Cmp(y.v) == 0
}

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

// Other math

func (x Int) Square() Int {
	return x.Mul(x)
}
func (x Int) Cube() Int {
	y := new(big.Int)
	y.Mul(x.v, x.v)
	y.Mul(y, x.v)
	return Int{y}
}

func (x Int) Exp(k int64) Int {
	b := big.NewInt(k)
	return Int{b.Exp(x.v, b, nil)}
}

func (x Int) SqrtFloor() Int {
	if x.IsZero() {
		return x
	}
	b := uint(x.BitLen())

	// invariant lo <= sqrt(x) < hi
	lo := One.Lsh((b - 1) / 2)
	hi := lo.Lsh(1)
	for {
		m := lo.Add(hi).Rsh(1)
		if m.Cmp(lo) == 0 {
			return lo
		}
		if m.Square().Cmp(x) <= 0 {
			lo = m
		} else {
			hi = m
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

// Discrete math stuff

func (x Int) ExpMod(k, m Int) Int {
	return Int{new(big.Int).Exp(x.v, k.v, m.v)}
}

func (x Int) ModInv(n Int) Int {
	// TODO: check gcd(x,n)==1?
	return Int{new(big.Int).ModInverse(x.v, n.v)}
}

func (x Int) GCD(y Int) Int {
	return Int{new(big.Int).GCD(nil, nil, x.v, y.v)}
}

// For printing

func (x Int) Format(s fmt.State, ch rune) {
	x.v.Format(s, ch)
}

// Rand returns a random number in [0,x)
func (x Int) Rand(rnd *rand.Rand) Int {
	return Int{new(big.Int).Rand(rnd, x.v)}
}

// Optimized routines

// Scratch space for use by Mod64s.  Mod64s is the same
// as Mod64 except it uses the scratch space to avoid allocation.
type Scratch [3]big.Int

func (x Int) Mod64s(y int64, s *Scratch) int64 {
	// Note: use DivMod here instead of Mod so we can reuse
	// storage for the dividend.  Mod allocates storage for
	// the (thrown away) dividend on each call.
	s[0].DivMod(x.v, s[1].SetInt64(y), &s[2])
	return s[2].Int64()
}

// helpful constants
var Zero = Int64(0)
var One = Int64(1)
var Two = Int64(2)
var Three = Int64(3)
var Ten = Int64(10)

var minInt64 = Int64(math.MinInt64)
