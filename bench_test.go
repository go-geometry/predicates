package predicates

import "testing"

// Benchmarks: the adaptive predicates on inputs the filter decides (random)
// and on inputs that run to stage D (exactly degenerate under inexact
// translations), against the Fast and Exact variants. bench/ compares them
// with predicates.c through cgo.

var sink float64

func benchInputs2() (rnd, deg [][4][2]float64) {
	r := newRand(21)
	for i := 0; i < 256; i++ {
		rnd = append(rnd, [4][2]float64{randPt2(r, 1), randPt2(r, 1), randPt2(r, 1), randPt2(r, 1)})
		q := [2]float64{randExp(r, -30, 60), randExp(r, -30, 60)}
		deg = append(deg, [4][2]float64{{randExp(r, -30, 60), randExp(r, -30, 60)}, q, q, {randExp(r, -30, 60), randExp(r, -30, 60)}})
	}
	return
}

func benchInputs3() (rnd, deg [][5][3]float64) {
	r := newRand(22)
	for i := 0; i < 256; i++ {
		var p [5][3]float64
		for k := range p {
			p[k] = [3]float64{randPt2(r, 1)[0], randPt2(r, 1)[0], randPt2(r, 1)[0]}
		}
		rnd = append(rnd, p)
		q := randExp3(r, -30, 60)
		deg = append(deg, [5][3]float64{randExp3(r, -30, 60), q, q, randExp3(r, -30, 60), q})
	}
	return
}

func BenchmarkOrient2D(b *testing.B) {
	rnd, deg := benchInputs2()
	for _, c := range []struct {
		name string
		in   [][4][2]float64
		f    func(a, b, c [2]float64) float64
	}{
		{"fast/random", rnd, Orient2DFast},
		{"adaptive/random", rnd, Orient2D},
		{"exact/random", rnd, Orient2DExact},
		{"adaptive/degenerate", deg, Orient2D},
		{"exact/degenerate", deg, Orient2DExact},
	} {
		b.Run(c.name, func(b *testing.B) {
			b.ReportAllocs()
			for i := 0; b.Loop(); i++ {
				p := &c.in[i&255]
				sink = c.f(p[0], p[1], p[2])
			}
		})
	}
}

func BenchmarkInCircle(b *testing.B) {
	rnd, deg := benchInputs2()
	for _, c := range []struct {
		name string
		in   [][4][2]float64
		f    func(a, b, c, d [2]float64) float64
	}{
		{"fast/random", rnd, InCircleFast},
		{"adaptive/random", rnd, InCircle},
		{"exact/random", rnd, InCircleExact},
		{"adaptive/degenerate", deg, InCircle},
		{"exact/degenerate", deg, InCircleExact},
	} {
		b.Run(c.name, func(b *testing.B) {
			b.ReportAllocs()
			for i := 0; b.Loop(); i++ {
				p := &c.in[i&255]
				sink = c.f(p[0], p[1], p[2], p[3])
			}
		})
	}
}

func BenchmarkOrient3D(b *testing.B) {
	rnd, deg := benchInputs3()
	for _, c := range []struct {
		name string
		in   [][5][3]float64
		f    func(a, b, c, d [3]float64) float64
	}{
		{"fast/random", rnd, Orient3DFast},
		{"adaptive/random", rnd, Orient3D},
		{"exact/random", rnd, Orient3DExact},
		{"adaptive/degenerate", deg, Orient3D},
		{"exact/degenerate", deg, Orient3DExact},
	} {
		b.Run(c.name, func(b *testing.B) {
			b.ReportAllocs()
			for i := 0; b.Loop(); i++ {
				p := &c.in[i&255]
				sink = c.f(p[0], p[1], p[2], p[3])
			}
		})
	}
}

func BenchmarkInSphere(b *testing.B) {
	rnd, deg := benchInputs3()
	for _, c := range []struct {
		name string
		in   [][5][3]float64
		f    func(a, b, c, d, e [3]float64) float64
	}{
		{"fast/random", rnd, InSphereFast},
		{"adaptive/random", rnd, InSphere},
		{"exact/random", rnd, InSphereExact},
		{"adaptive/degenerate", deg, InSphere},
		{"exact/degenerate", deg, InSphereExact},
	} {
		b.Run(c.name, func(b *testing.B) {
			b.ReportAllocs()
			for i := 0; b.Loop(); i++ {
				p := &c.in[i&255]
				sink = c.f(p[0], p[1], p[2], p[3], p[4])
			}
		})
	}
}

func BenchmarkOrientN(b *testing.B) {
	rnd, _ := benchInputs3()
	pts := make([][][]float64, len(rnd))
	for i, p := range rnd {
		pts[i] = [][]float64{p[0][:], p[1][:], p[2][:], p[3][:]}
	}
	b.ReportAllocs()
	for i := 0; b.Loop(); i++ {
		sink, _ = OrientN(pts[i&255])
	}
}
