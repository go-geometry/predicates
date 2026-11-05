package predicates

import (
	"math"
	"testing"
)

// The fuzz targets compare each predicate's sign with the math/big oracle on
// arbitrary float64 inputs (including NaN, ±Inf, subnormals and huge
// values). The seed corpus lives in testdata/fuzz.

func allFinite(xs ...float64) bool { return finite(xs...) }

func FuzzOrient2D(f *testing.F) {
	f.Add(0.5, 0.5, 12.0, 12.0, 24.0, 24.0)
	f.Add(0.0, 0.0, 1e300, 1e300, -1e300, -1e300)
	f.Add(math.SmallestNonzeroFloat64, 0.0, 0.0, math.SmallestNonzeroFloat64, 1.0, 1.0)
	f.Fuzz(func(t *testing.T, ax, ay, bx, by, cx, cy float64) {
		a, b, c := [2]float64{ax, ay}, [2]float64{bx, by}, [2]float64{cx, cy}
		got := Orient2D(a, b, c)
		if !allFinite(ax, ay, bx, by, cx, cy) {
			if !math.IsNaN(got) {
				t.Fatalf("non-finite input: got %v, want NaN", got)
			}
			return
		}
		if want := oracleOrient2D(a, b, c); sign(got) != want || sign(Orient2DExact(a, b, c)) != want {
			t.Fatalf("Orient2D(%v, %v, %v) = %v, exact sign %d", a, b, c, got, want)
		}
	})
}

func FuzzInCircle(f *testing.F) {
	f.Add(1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0)
	f.Add(5.0, 0.0, 4.0, 3.0, -3.0, 4.0, 0.0, -5.0)
	f.Add(1e200, 0.0, 0.0, 1e-200, -1e100, 3.0, 0.1, 0.2)
	f.Fuzz(func(t *testing.T, ax, ay, bx, by, cx, cy, dx, dy float64) {
		a, b, c, d := [2]float64{ax, ay}, [2]float64{bx, by}, [2]float64{cx, cy}, [2]float64{dx, dy}
		got := InCircle(a, b, c, d)
		if !allFinite(ax, ay, bx, by, cx, cy, dx, dy) {
			if !math.IsNaN(got) {
				t.Fatalf("non-finite input: got %v, want NaN", got)
			}
			return
		}
		if want := oracleInCircle(a, b, c, d); sign(got) != want || sign(InCircleExact(a, b, c, d)) != want {
			t.Fatalf("InCircle(%v, %v, %v, %v) = %v, exact sign %d", a, b, c, d, got, want)
		}
	})
}

func FuzzOrient3D(f *testing.F) {
	f.Add(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 5.0, 7.0, 0.0)
	f.Add(0.5, 0.5, 0.5, 12.0, 12.0, 12.0, 24.0, 24.0, 24.0, 0.0, 0.0, 1.0)
	f.Fuzz(func(t *testing.T, ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz float64) {
		a, b, c, d := [3]float64{ax, ay, az}, [3]float64{bx, by, bz}, [3]float64{cx, cy, cz}, [3]float64{dx, dy, dz}
		got := Orient3D(a, b, c, d)
		if !allFinite(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz) {
			if !math.IsNaN(got) {
				t.Fatalf("non-finite input: got %v, want NaN", got)
			}
			return
		}
		if want := oracleOrient3D(a, b, c, d); sign(got) != want || sign(Orient3DExact(a, b, c, d)) != want {
			t.Fatalf("Orient3D(%v, %v, %v, %v) = %v, exact sign %d", a, b, c, d, got, want)
		}
	})
}

// FuzzInSphere fixes nothing: 15 coordinates.
func FuzzInSphere(f *testing.F) {
	f.Add(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0)
	f.Add(3.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 2.0, 0.0, -3.0, 0.0, -2.0, -2.0, 1.0)
	f.Fuzz(func(t *testing.T, ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez float64) {
		a, b, c := [3]float64{ax, ay, az}, [3]float64{bx, by, bz}, [3]float64{cx, cy, cz}
		d, e := [3]float64{dx, dy, dz}, [3]float64{ex, ey, ez}
		got := InSphere(a, b, c, d, e)
		if !allFinite(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez) {
			if !math.IsNaN(got) {
				t.Fatalf("non-finite input: got %v, want NaN", got)
			}
			return
		}
		if want := oracleInSphere(a, b, c, d, e); sign(got) != want || sign(InSphereExact(a, b, c, d, e)) != want {
			t.Fatalf("InSphere(%v, %v, %v, %v, %v) = %v, exact sign %d", a, b, c, d, e, got, want)
		}
	})
}

// FuzzNDim decodes a dimension and a flat coordinate list from the input
// bytes, 8 bytes per coordinate, and checks OrientN and InSphereN.
func FuzzNDim(f *testing.F) {
	f.Add(uint8(2), []byte{})
	f.Add(uint8(4), make([]byte, 8*30))
	f.Fuzz(func(t *testing.T, dim uint8, raw []byte) {
		d := int(dim%4) + 1
		need := d * (d + 2)
		if len(raw) < 8*need {
			return
		}
		coord := func(i int) float64 {
			var u uint64
			for k := 0; k < 8; k++ {
				u = u<<8 | uint64(raw[8*i+k])
			}
			return math.Float64frombits(u)
		}
		pts := make([][]float64, d+2)
		ok := true
		for i := range pts {
			pts[i] = make([]float64, d)
			for j := range pts[i] {
				pts[i][j] = coord(i*d + j)
			}
			ok = ok && allFinite(pts[i]...)
		}
		o, errO := OrientN(pts[:d+1])
		s, errS := InSphereN(pts)
		if !ok {
			if errO == nil && !allFinite(flat(pts[:d+1])...) || errS == nil {
				t.Fatalf("non-finite input accepted: %v %v", errO, errS)
			}
			return
		}
		if errO != nil || errS != nil {
			t.Fatalf("finite input refused: %v %v", errO, errS)
		}
		if want := oracleOrientN(pts[:d+1]); sign(o) != want {
			t.Fatalf("OrientN(%v) = %v, exact sign %d", pts[:d+1], o, want)
		}
		if want := oracleInSphereN(pts); sign(s) != want {
			t.Fatalf("InSphereN(%v) = %v, exact sign %d", pts, s, want)
		}
	})
}

func flat(pts [][]float64) []float64 {
	var out []float64
	for _, p := range pts {
		out = append(out, p...)
	}
	return out
}
