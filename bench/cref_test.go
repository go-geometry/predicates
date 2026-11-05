package bench

import (
	"math"
	"math/rand/v2"
	"testing"

	"github.com/go-geometry/predicates"
)

const records = 256

func randExp(r *rand.Rand) float64 {
	x := math.Ldexp(1+r.Float64(), r.IntN(91)-30)
	if r.IntN(2) == 0 {
		x = -x
	}
	return x
}

// inputs returns records of `per` float64: "random" points in [-1, 1), which
// the floating-point filter decides, and "degenerate" records in which one
// point is repeated under inexact translations (mixed exponents), which run
// to the exact stage D.
func inputs(per int, degenerate bool, repeat [][2]int) []float64 {
	r := rand.New(rand.NewPCG(31, uint64(per)))
	pts := make([]float64, records*per)
	for i := range pts {
		if degenerate {
			pts[i] = randExp(r)
		} else {
			pts[i] = 2*r.Float64() - 1
		}
	}
	if degenerate {
		dim := 2
		if per == 15 {
			dim = 3
		}
		for k := 0; k < records; k++ {
			rec := pts[k*per : (k+1)*per]
			for _, rp := range repeat {
				copy(rec[rp[1]*dim:(rp[1]+1)*dim], rec[rp[0]*dim:(rp[0]+1)*dim])
			}
		}
	}
	return pts
}

type pred struct {
	name   string
	kind   Kind
	per    int
	repeat [][2]int // point copies that make the degenerate records exact zeros
	goF    func(p []float64) float64
}

func p2(p []float64, i int) [2]float64 { return [2]float64{p[2*i], p[2*i+1]} }
func p3(p []float64, i int) [3]float64 { return [3]float64{p[3*i], p[3*i+1], p[3*i+2]} }

var preds = []pred{
	{"Orient2D", Orient2D, 8, [][2]int{{1, 2}}, func(p []float64) float64 {
		return predicates.Orient2D(p2(p, 0), p2(p, 1), p2(p, 2))
	}},
	{"InCircle", InCircle, 8, [][2]int{{1, 2}}, func(p []float64) float64 {
		return predicates.InCircle(p2(p, 0), p2(p, 1), p2(p, 2), p2(p, 3))
	}},
	{"Orient3D", Orient3D, 15, [][2]int{{1, 2}}, func(p []float64) float64 {
		return predicates.Orient3D(p3(p, 0), p3(p, 1), p3(p, 2), p3(p, 3))
	}},
	{"InSphere", InSphere, 15, [][2]int{{1, 2}}, func(p []float64) float64 {
		return predicates.InSphere(p3(p, 0), p3(p, 1), p3(p, 2), p3(p, 3), p3(p, 4))
	}},
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

// The Go predicates and predicates.c agree in sign on every benchmark record.
func TestSignsAgreeWithPredicatesC(t *testing.T) {
	for _, pr := range preds {
		for _, deg := range []bool{false, true} {
			pts := inputs(pr.per, deg, pr.repeat)
			for k := 0; k < records; k++ {
				rec := pts[k*pr.per : (k+1)*pr.per]
				c := Loop(pr.kind, rec, 1, 1)
				g := pr.goF(rec)
				if sign(c) != sign(g) {
					t.Fatalf("%s record %d: C %v, Go %v", pr.name, k, c, g)
				}
			}
		}
	}
}

var sink float64

// goLoop mirrors loops.c: the records are converted to the Go argument types
// once, then each predicate is called directly (no function value), so both
// sides time the same work.
func goLoop(k Kind, pts []float64, iters int) float64 {
	s := 0.0
	switch k {
	case Orient2D, InCircle:
		recs := make([][4][2]float64, records)
		for i := range recs {
			p := pts[8*i:]
			recs[i] = [4][2]float64{p2(p, 0), p2(p, 1), p2(p, 2), p2(p, 3)}
		}
		if k == Orient2D {
			for i := 0; i < iters; i++ {
				r := &recs[i%records]
				s += predicates.Orient2D(r[0], r[1], r[2])
			}
			return s
		}
		for i := 0; i < iters; i++ {
			r := &recs[i%records]
			s += predicates.InCircle(r[0], r[1], r[2], r[3])
		}
		return s
	}
	recs := make([][5][3]float64, records)
	for i := range recs {
		p := pts[15*i:]
		recs[i] = [5][3]float64{p3(p, 0), p3(p, 1), p3(p, 2), p3(p, 3), p3(p, 4)}
	}
	if k == Orient3D {
		for i := 0; i < iters; i++ {
			r := &recs[i%records]
			s += predicates.Orient3D(r[0], r[1], r[2], r[3])
		}
		return s
	}
	for i := 0; i < iters; i++ {
		r := &recs[i%records]
		s += predicates.InSphere(r[0], r[1], r[2], r[3], r[4])
	}
	return s
}

func BenchmarkVsPredicatesC(b *testing.B) {
	for _, pr := range preds {
		for _, deg := range []bool{false, true} {
			kind := "random"
			if deg {
				kind = "degenerate"
			}
			pts := inputs(pr.per, deg, pr.repeat)
			b.Run(pr.name+"/"+kind+"/c", func(b *testing.B) {
				sink = Loop(pr.kind, pts, records, b.N)
			})
			b.Run(pr.name+"/"+kind+"/go", func(b *testing.B) {
				b.ReportAllocs()
				sink = goLoop(pr.kind, pts, b.N)
			})
		}
	}
}
