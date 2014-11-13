package factorlib

import (
	"github.com/randall77/factorlib/big"
	"testing"
)

func TestElliptic(t *testing.T) {
	// Test that for prime n, points on an elliptic curve form a group
	const n = 17
	const a = 5

	N := big.Int64(n)
	A := big.Int64(a)

	// find all the points on the curve
	onCurve := []point{{true, big.Int{}, big.Int{}}}
	for x := 0; x < n; x++ {
		for y := 0; y < n; y++ {
			p := point{false, big.Int64(int64(x)), big.Int64(int64(y))}
			if p.Check(N, A) {
				onCurve = append(onCurve, p)
			}
		}
	}

	// make sure addition is closed
	for _, p := range onCurve {
	testloop:
		for _, q := range onCurve {
			r := p.Add(q, N, A)
			for _, s := range onCurve {
				if r.Equals(s) {
					continue testloop
				}
			}
			t.Fatalf("(%t %d %d)+(%t %d %d)=(%t %d %d) not on curve", p.inf, p.x, p.y, q.inf, q.x, q.y, r.inf, r.x, r.y)
		}
	}
}
