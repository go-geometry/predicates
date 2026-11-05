package predicates

import (
	"math"
	"testing"
)

// tally counts outcomes per corpus kind so that a test can show it exercised
// every kind, and how often the non-robust evaluation was wrong.
type tally struct {
	cases, fastWrong map[kind]int
}

func newTally() *tally { return &tally{map[kind]int{}, map[kind]int{}} }

func (tl *tally) log(t *testing.T, name string) {
	t.Helper()
	for _, k := range []kind{kRandom, kNear, kExact, kTiltedGrid, kHuge} {
		t.Logf("%s %-20s cases %6d   Fast wrong %5d", name, k, tl.cases[k], tl.fastWrong[k])
	}
}

func TestOrient2DSignMatchesOracle(t *testing.T) {
	tl := newTally()
	for _, c := range corpus2(t, 600) {
		a, b, cc := c.p[0], c.p[1], c.p[2]
		want := oracleOrient2D(a, b, cc)
		if got := sign(Orient2D(a, b, cc)); got != want {
			t.Fatalf("Orient2D(%v, %v, %v) sign %d, exact %d (%s)", a, b, cc, got, want, c.kind)
		}
		if got := sign(Orient2DExact(a, b, cc)); got != want {
			t.Fatalf("Orient2DExact(%v, %v, %v) sign %d, exact %d (%s)", a, b, cc, got, want, c.kind)
		}
		tl.cases[c.kind]++
		if sign(Orient2DFast(a, b, cc)) != want {
			tl.fastWrong[c.kind]++
		}
	}
	tl.log(t, "Orient2D")
	// The corpus must be hard: the plain evaluation fails on it.
	if tl.fastWrong[kTiltedGrid] == 0 || tl.fastWrong[kNear] == 0 {
		t.Fatalf("near-degenerate corpus too easy: Fast never wrong (%v)", tl.fastWrong)
	}
}

func TestInCircleSignMatchesOracle(t *testing.T) {
	tl := newTally()
	for _, c := range corpus2(t, 300) {
		a, b, cc, d := c.p[0], c.p[1], c.p[2], c.p[3]
		want := oracleInCircle(a, b, cc, d)
		if got := sign(InCircle(a, b, cc, d)); got != want {
			t.Fatalf("InCircle(%v, %v, %v, %v) sign %d, exact %d (%s)", a, b, cc, d, got, want, c.kind)
		}
		if got := sign(InCircleExact(a, b, cc, d)); got != want {
			t.Fatalf("InCircleExact(%v, %v, %v, %v) sign %d, exact %d (%s)", a, b, cc, d, got, want, c.kind)
		}
		tl.cases[c.kind]++
		if sign(InCircleFast(a, b, cc, d)) != want {
			tl.fastWrong[c.kind]++
		}
	}
	tl.log(t, "InCircle")
	if tl.fastWrong[kNear] == 0 {
		t.Fatalf("near-degenerate corpus too easy: Fast never wrong (%v)", tl.fastWrong)
	}
}

func TestDocumentedSigns2D(t *testing.T) {
	a, b, c := [2]float64{0, 0}, [2]float64{1, 0}, [2]float64{0, 1}
	if Orient2D(a, b, c) <= 0 || Orient2D(a, c, b) >= 0 || Orient2D(a, b, [2]float64{2, 0}) != 0 {
		t.Fatal("Orient2D: counterclockwise must be positive, clockwise negative, collinear zero")
	}
	p, q, r := [2]float64{1, 0}, [2]float64{0, 1}, [2]float64{-1, 0}
	if InCircle(p, q, r, [2]float64{0, 0}) <= 0 || InCircle(p, q, r, [2]float64{2, 2}) >= 0 ||
		InCircle(p, q, r, [2]float64{0, -1}) != 0 {
		t.Fatal("InCircle: inside must be positive, outside negative, cocircular zero")
	}
}

func TestNonFiniteYieldsNaN2D(t *testing.T) {
	for _, bad := range []float64{math.NaN(), math.Inf(1), math.Inf(-1)} {
		p := [2]float64{bad, 0}
		o := [2]float64{0, 0}
		u := [2]float64{1, 0}
		v := [2]float64{0, 1}
		for _, got := range []float64{
			Orient2D(p, u, v), Orient2D(o, u, p), Orient2DExact(u, p, v),
			InCircle(p, u, v, o), InCircle(o, u, v, p), InCircleExact(o, p, u, v),
		} {
			if !math.IsNaN(got) {
				t.Fatalf("non-finite input %v: got %v, want NaN", bad, got)
			}
		}
	}
}

// The fallback keeps the sign when the determinant underflows float64 and
// saturates to ±Inf when it overflows.
func TestBigResultMagnitude(t *testing.T) {
	tiny := math.SmallestNonzeroFloat64
	if got := Orient2D([2]float64{0, 0}, [2]float64{tiny, 0}, [2]float64{0, tiny}); got != tiny {
		t.Fatalf("underflowing determinant: got %v, want the smallest subnormal", got)
	}
	if got := Orient2D([2]float64{0, 0}, [2]float64{0, tiny}, [2]float64{tiny, 0}); got != -tiny {
		t.Fatalf("underflowing determinant: got %v, want minus the smallest subnormal", got)
	}
	h := math.MaxFloat64
	if got := Orient2D([2]float64{0, 0}, [2]float64{h, 0}, [2]float64{0, h}); !math.IsInf(got, 1) {
		t.Fatalf("overflowing determinant: got %v, want +Inf", got)
	}
}

// bigDet's pivoting: a zero leading entry with a nonzero one below, and a
// zero column.
func TestBigDetPivoting(t *testing.T) {
	one := func(x float64) *bigInt { return toBig(x) }
	m := [][]*bigInt{{one(0), one(1)}, {one(1), one(0)}}
	if got := bigDet(m, 2); got != -1 {
		t.Fatalf("det [[0 1] [1 0]] = %v, want -1", got)
	}
	m = [][]*bigInt{{one(0), one(1)}, {one(0), one(2)}}
	if got := bigDet(m, 2); got != 0 {
		t.Fatalf("singular: got %v, want 0", got)
	}
}

// Planted breaks: with the error bound of a stage removed, the same corpus
// must produce wrong signs. If it does not, the corpus cannot see that stage.
func TestPlantedBreakFilterA(t *testing.T) {
	cs := corpus2(t, 300)
	defer func(a, b float64) { ccwErrBoundA, iccErrBoundA = a, b }(ccwErrBoundA, iccErrBoundA)
	ccwErrBoundA, iccErrBoundA = 0, 0
	wrongO, wrongI := 0, 0
	for _, c := range cs {
		if sign(Orient2D(c.p[0], c.p[1], c.p[2])) != oracleOrient2D(c.p[0], c.p[1], c.p[2]) {
			wrongO++
		}
		if sign(InCircle(c.p[0], c.p[1], c.p[2], c.p[3])) != oracleInCircle(c.p[0], c.p[1], c.p[2], c.p[3]) {
			wrongI++
		}
	}
	t.Logf("filter A removed: Orient2D wrong %d, InCircle wrong %d of %d", wrongO, wrongI, len(cs))
	if wrongO == 0 || wrongI == 0 {
		t.Fatal("planted break not detected: the corpus cannot see stage A's error bound")
	}
}

func TestPlantedBreakStageB(t *testing.T) {
	cs := corpus2(t, 300)
	defer func(a, b float64) { ccwErrBoundB, iccErrBoundB = a, b }(ccwErrBoundB, iccErrBoundB)
	ccwErrBoundB, iccErrBoundB = 0, 0
	wrongO, wrongI := 0, 0
	for _, c := range cs {
		if sign(Orient2D(c.p[0], c.p[1], c.p[2])) != oracleOrient2D(c.p[0], c.p[1], c.p[2]) {
			wrongO++
		}
		if sign(InCircle(c.p[0], c.p[1], c.p[2], c.p[3])) != oracleInCircle(c.p[0], c.p[1], c.p[2], c.p[3]) {
			wrongI++
		}
	}
	t.Logf("stage B bound removed: Orient2D wrong %d, InCircle wrong %d of %d", wrongO, wrongI, len(cs))
	if wrongO == 0 || wrongI == 0 {
		t.Fatal("planted break not detected: the corpus cannot see stage B's error bound")
	}
}
