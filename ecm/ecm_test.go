package ecm

import (
	"testing"

	"github.com/randall77/factorlib/big"
)

func TestElliptic(t *testing.T) {
	// Test that for prime n, points on an elliptic curve form a group
	const n = 17
	const a = 5

	N := big.Int64(n)
	A := big.Int64(a)

	// find all the points on the curve
	onCurve := []point{{big.Zero, big.Zero}}
	for x := 0; x < n; x++ {
		for y := 0; y < n; y++ {
			p := point{big.Int64(int64(x)), big.Int64(int64(y))}
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
			t.Fatalf("(%d %d)+(%d %d)=(%d %d) not on curve", p.x, p.y, q.x, q.y, r.x, r.y)
		}
	}
}
