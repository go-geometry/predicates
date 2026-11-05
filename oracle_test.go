package predicates

import (
	"math"
	"math/big"
	"math/rand/v2"
	"testing"
)

// The oracle: exact rational evaluation with math/big. It shares no code with
// the package (big.Rat, cofactor expansion; the package's fallback uses
// big.Int and Bareiss elimination).

func rat(x float64) *big.Rat { return new(big.Rat).SetFloat64(x) }

func ratSub(a, b float64) *big.Rat { return new(big.Rat).Sub(rat(a), rat(b)) }

// ratDet is the determinant by cofactor expansion along the first row.
func ratDet(m [][]*big.Rat) *big.Rat {
	n := len(m)
	if n == 1 {
		return new(big.Rat).Set(m[0][0])
	}
	sum := new(big.Rat)
	for j := 0; j < n; j++ {
		if m[0][j].Sign() == 0 {
			continue
		}
		sub := make([][]*big.Rat, n-1)
		for i := 1; i < n; i++ {
			row := make([]*big.Rat, 0, n-1)
			row = append(row, m[i][:j]...)
			row = append(row, m[i][j+1:]...)
			sub[i-1] = row
		}
		t := new(big.Rat).Mul(m[0][j], ratDet(sub))
		if j%2 == 1 {
			sum.Sub(sum, t)
		} else {
			sum.Add(sum, t)
		}
	}
	return sum
}

// ratLiftedRow returns [p-q..., |p-q|²].
func ratLiftedRow(p, q []float64) []*big.Rat {
	row := make([]*big.Rat, 0, len(p)+1)
	lift := new(big.Rat)
	for i := range p {
		d := ratSub(p[i], q[i])
		row = append(row, d)
		lift.Add(lift, new(big.Rat).Mul(d, d))
	}
	return append(row, lift)
}

// ratDiffRow returns [p-q...].
func ratDiffRow(p, q []float64) []*big.Rat {
	row := make([]*big.Rat, len(p))
	for i := range p {
		row[i] = ratSub(p[i], q[i])
	}
	return row
}

func oracleOrient2D(a, b, c [2]float64) int {
	return ratDet([][]*big.Rat{ratDiffRow(a[:], c[:]), ratDiffRow(b[:], c[:])}).Sign()
}

func oracleInCircle(a, b, c, d [2]float64) int {
	return ratDet([][]*big.Rat{ratLiftedRow(a[:], d[:]), ratLiftedRow(b[:], d[:]), ratLiftedRow(c[:], d[:])}).Sign()
}

func sign(x float64) int {
	switch {
	case x > 0:
		return 1
	case x < 0:
		return -1
	}
	return 0
}

// --- corpora -------------------------------------------------------------

// nudge moves x by k ulps (k may be negative).
func nudge(x float64, k int) float64 {
	for ; k > 0; k-- {
		x = math.Nextafter(x, math.Inf(1))
	}
	for ; k < 0; k++ {
		x = math.Nextafter(x, math.Inf(-1))
	}
	return x
}

func newRand(seed uint64) *rand.Rand { return rand.New(rand.NewPCG(seed, 0x9e3779b97f4a7c15)) }

func randPt2(r *rand.Rand, scale float64) [2]float64 {
	return [2]float64{(2*r.Float64() - 1) * scale, (2*r.Float64() - 1) * scale}
}

// randExp returns a random float64 with a random exponent in [lo, hi].
func randExp(r *rand.Rand, lo, hi int) float64 {
	x := math.Ldexp(1+r.Float64(), lo+r.IntN(hi-lo+1))
	if r.IntN(2) == 0 {
		x = -x
	}
	return x
}

// corpus2 yields triples of 2-D points of every kind the tests need.
type kind string

const (
	kRandom     kind = "random"
	kNear       kind = "near-degenerate"
	kExact      kind = "exactly-degenerate"
	kTiltedGrid kind = "tilted-grid"
	kHuge       kind = "huge/tiny"
)

type case2 struct {
	kind kind
	p    [4][2]float64 // up to four points
}

// collinear2 returns three points on a line, the third rounded onto it and
// then nudged by up to ±k ulps.
func collinear2(r *rand.Rand, scale float64, k int) [3][2]float64 {
	a, b := randPt2(r, scale), randPt2(r, scale)
	t := r.Float64()*3 - 1
	c := [2]float64{a[0] + t*(b[0]-a[0]), a[1] + t*(b[1]-a[1])}
	if k > 0 {
		c[0] = nudge(c[0], r.IntN(2*k+1)-k)
		c[1] = nudge(c[1], r.IntN(2*k+1)-k)
	}
	return [3][2]float64{a, b, c}
}

// cocircular2 returns four points on a circle (rounded), possibly nudged.
func cocircular2(r *rand.Rand, k int) [4][2]float64 {
	cx, cy, rad := r.Float64(), r.Float64(), 0.5+r.Float64()
	var p [4][2]float64
	// Counterclockwise order for the first three.
	base := r.Float64() * 2 * math.Pi
	angles := [4]float64{base, base + 2, base + 4, base + 1 + r.Float64()}
	for i, th := range angles {
		p[i] = [2]float64{cx + rad*math.Cos(th), cy + rad*math.Sin(th)}
	}
	if k > 0 {
		p[3][0] = nudge(p[3][0], r.IntN(2*k+1)-k)
		p[3][1] = nudge(p[3][1], r.IntN(2*k+1)-k)
	}
	return p
}

// exactCircle returns four points exactly on the circle x²+y² = 25, scaled by
// 2^s and translated by an integer offset (all exact), in counterclockwise
// order for the first three.
func exactCircle(r *rand.Rand) [4][2]float64 {
	pts := [][2]float64{{5, 0}, {4, 3}, {3, 4}, {0, 5}, {-3, 4}, {-4, 3}, {-5, 0}, {-4, -3}, {-3, -4}, {0, -5}, {3, -4}, {4, -3}}
	i := r.IntN(4)
	idx := [4]int{i, i + 4, i + 8, (i + 2 + r.IntN(10)) % 12}
	s := math.Ldexp(1, r.IntN(60)-30)
	ox, oy := float64(r.IntN(2001)-1000), float64(r.IntN(2001)-1000)
	var p [4][2]float64
	for k, j := range idx {
		p[k] = [2]float64{pts[j][0]*s + ox, pts[j][1]*s + oy}
	}
	return p
}

func corpus2(t testing.TB, n int) []case2 {
	t.Helper()
	r := newRand(1)
	var cs []case2
	add := func(k kind, p ...[2]float64) {
		var c case2
		c.kind = k
		copy(c.p[:], p)
		cs = append(cs, c)
	}
	for i := 0; i < n; i++ {
		add(kRandom, randPt2(r, 1), randPt2(r, 1), randPt2(r, 1), randPt2(r, 1))
		q := collinear2(r, 1, 0)
		add(kNear, q[0], q[1], q[2], randPt2(r, 1))
		q = collinear2(r, 1e3, 3)
		add(kNear, q[0], q[1], q[2], randPt2(r, 1))
		cc := cocircular2(r, i%4)
		add(kNear, cc[0], cc[1], cc[2], cc[3])
		ec := exactCircle(r)
		add(kExact, ec[0], ec[1], ec[2], ec[3])
		// Exactly collinear: integer multiples of a direction, translated.
		dx, dy := float64(r.IntN(1000)-500), float64(r.IntN(1000)-500)
		ox, oy := float64(r.IntN(1<<20)), float64(r.IntN(1<<20))
		add(kExact, [2]float64{ox, oy}, [2]float64{ox + dx, oy + dy}, [2]float64{ox + 3*dx, oy + 3*dy}, [2]float64{ox - dx, oy - dy})
		// A repeated point under inexact translations (mixed exponents): the
		// determinant is exactly zero, yet every rounded stage sees noise, so
		// the evaluation runs down to stage D.
		pq := [2]float64{randExp(r, -30, 60), randExp(r, -30, 60)}
		add(kExact, [2]float64{randExp(r, -30, 60), randExp(r, -30, 60)}, pq, pq,
			[2]float64{randExp(r, -30, 60), randExp(r, -30, 60)})
		// Huge, tiny and mixed magnitudes (the out-of-range path).
		add(kHuge, [2]float64{randExp(r, 900, 1000), randExp(r, 900, 1000)},
			[2]float64{randExp(r, 900, 1000), randExp(r, 900, 1000)},
			[2]float64{randExp(r, 900, 1000), randExp(r, 900, 1000)},
			[2]float64{randExp(r, 900, 1000), randExp(r, 900, 1000)})
		add(kHuge, [2]float64{randExp(r, -1074, -1000), randExp(r, -1074, -1000)},
			[2]float64{randExp(r, -1074, -1000), randExp(r, -1074, -1000)},
			[2]float64{randExp(r, -1074, -1000), randExp(r, -1074, -1000)},
			[2]float64{randExp(r, -1074, -1000), randExp(r, -1074, -1000)})
		add(kHuge, [2]float64{randExp(r, -600, 600), randExp(r, -600, 600)},
			[2]float64{randExp(r, -600, 600), randExp(r, -600, 600)},
			[2]float64{randExp(r, -600, 600), randExp(r, -600, 600)},
			[2]float64{randExp(r, -600, 600), randExp(r, -600, 600)})
		// Huge but exactly collinear: the fallback must answer zero.
		s := math.Ldexp(1, 900+r.IntN(100))
		add(kHuge, [2]float64{s, 2 * s}, [2]float64{2 * s, 4 * s}, [2]float64{-s, -2 * s}, [2]float64{0, 0})
	}
	// Shewchuk's tilted grid: points one ulp apart near (0.5, 0.5), against
	// the line through (12, 12) and (24, 24); Orient2DFast misclassifies many.
	for i := 0; i < 64; i++ {
		for j := 0; j < 64; j += 7 {
			p := [2]float64{nudge(0.5, i), nudge(0.5, j)}
			add(kTiltedGrid, p, [2]float64{12, 12}, [2]float64{24, 24}, [2]float64{nudge(0.5, j), nudge(0.5, i)})
		}
	}
	return cs
}
