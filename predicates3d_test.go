package predicates

import (
	"math"
	"math/big"
	"math/rand/v2"
	"testing"
)

func oracleOrient3D(a, b, c, d [3]float64) int {
	return ratDet([][]*big.Rat{ratDiffRow(a[:], d[:]), ratDiffRow(b[:], d[:]), ratDiffRow(c[:], d[:])}).Sign()
}

func oracleInSphere(a, b, c, d, e [3]float64) int {
	return ratDet([][]*big.Rat{ratLiftedRow(a[:], e[:]), ratLiftedRow(b[:], e[:]), ratLiftedRow(c[:], e[:]), ratLiftedRow(d[:], e[:])}).Sign()
}

type case3 struct {
	kind kind
	p    [5][3]float64
}

func randPt3(r *rand.Rand, scale float64) [3]float64 {
	return [3]float64{(2*r.Float64() - 1) * scale, (2*r.Float64() - 1) * scale, (2*r.Float64() - 1) * scale}
}

func nudge3(r *rand.Rand, p [3]float64, k int) [3]float64 {
	if k == 0 {
		return p
	}
	for i := range p {
		p[i] = nudge(p[i], r.IntN(2*k+1)-k)
	}
	return p
}

func randExp3(r *rand.Rand, lo, hi int) [3]float64 {
	return [3]float64{randExp(r, lo, hi), randExp(r, lo, hi), randExp(r, lo, hi)}
}

// exactSphere returns five points exactly on x²+y²+z² = 9, scaled by 2^s and
// translated by an integer offset, all exact.
func exactSphere(r *rand.Rand) [5][3]float64 {
	var on [][3]float64
	for _, v := range [][3]float64{{3, 0, 0}, {1, 2, 2}} {
		for _, perm := range [][3]int{{0, 1, 2}, {1, 2, 0}, {2, 0, 1}} {
			for sg := 0; sg < 8; sg++ {
				p := [3]float64{v[perm[0]], v[perm[1]], v[perm[2]]}
				for i := 0; i < 3; i++ {
					if sg&(1<<i) != 0 {
						p[i] = -p[i]
					}
				}
				on = append(on, p)
			}
		}
	}
	s := math.Ldexp(1, r.IntN(40)-20)
	off := [3]float64{float64(r.IntN(2001) - 1000), float64(r.IntN(2001) - 1000), float64(r.IntN(2001) - 1000)}
	var out [5][3]float64
	for k := range out {
		p := on[r.IntN(len(on))]
		for i := range p {
			out[k][i] = p[i]*s + off[i]
		}
	}
	return out
}

func corpus3(t testing.TB, n int) []case3 {
	t.Helper()
	r := newRand(3)
	var cs []case3
	add := func(k kind, p ...[3]float64) {
		var c case3
		c.kind = k
		copy(c.p[:], p)
		cs = append(cs, c)
	}
	for i := 0; i < n; i++ {
		add(kRandom, randPt3(r, 1), randPt3(r, 1), randPt3(r, 1), randPt3(r, 1), randPt3(r, 1))
		// Near-coplanar: d rounded onto the plane of a, b, c, then nudged.
		a, b, c := randPt3(r, 1), randPt3(r, 1), randPt3(r, 1)
		s, u := r.Float64()*2-0.5, r.Float64()*2-0.5
		var d [3]float64
		for k := range d {
			d[k] = a[k] + s*(b[k]-a[k]) + u*(c[k]-a[k])
		}
		add(kNear, a, b, c, nudge3(r, d, i%3), randPt3(r, 1))
		// Near-cospherical: five points rounded onto a sphere, the last nudged.
		ctr, rad := randPt3(r, 1), 0.5+r.Float64()
		var sp [5][3]float64
		for k := range sp {
			z, th := 2*r.Float64()-1, 2*math.Pi*r.Float64()
			w := math.Sqrt(1 - z*z)
			sp[k] = [3]float64{ctr[0] + rad*w*math.Cos(th), ctr[1] + rad*w*math.Sin(th), ctr[2] + rad*z}
		}
		sp[4] = nudge3(r, sp[4], i%4)
		add(kNear, sp[0], sp[1], sp[2], sp[3], sp[4])
		es := exactSphere(r)
		add(kExact, es[0], es[1], es[2], es[3], es[4])
		// Exactly coplanar lattice points (z = 0 plane after an integer shift).
		oz := float64(r.IntN(1 << 20))
		lat := func() [3]float64 {
			return [3]float64{float64(r.IntN(2001) - 1000), float64(r.IntN(2001) - 1000), oz}
		}
		add(kExact, lat(), lat(), lat(), lat(), lat())
		// A repeated point under inexact translations: exactly zero.
		q := randExp3(r, -20, 50)
		add(kExact, randExp3(r, -20, 50), q, q, randExp3(r, -20, 50), randExp3(r, -20, 50))
		add(kExact, randExp3(r, -20, 50), randExp3(r, -20, 50), q, randExp3(r, -20, 50), q)
		// Huge, tiny and mixed magnitudes.
		add(kHuge, randExp3(r, 900, 1000), randExp3(r, 900, 1000), randExp3(r, 900, 1000), randExp3(r, 900, 1000), randExp3(r, 900, 1000))
		add(kHuge, randExp3(r, -1074, -1000), randExp3(r, -1074, -1000), randExp3(r, -1074, -1000), randExp3(r, -1074, -1000), randExp3(r, -1074, -1000))
		add(kHuge, randExp3(r, -400, 400), randExp3(r, -400, 400), randExp3(r, -400, 400), randExp3(r, -400, 400), randExp3(r, -400, 400))
	}
	// A tilted grid in 3-D: points a few ulps apart near (0.5, 0.5, 0.5)
	// against the plane through three points of the line x = y = z.
	for i := 0; i < 16; i++ {
		for j := 0; j < 16; j += 3 {
			p := [3]float64{nudge(0.5, i), nudge(0.5, j), nudge(0.5, i-j)}
			add(kTiltedGrid, [3]float64{12, 12, 12}, [3]float64{24, 24, 24}, [3]float64{0, 0, 1}, p, [3]float64{nudge(0.5, j), 0.5, nudge(0.5, i)})
		}
	}
	return cs
}

func TestOrient3DSignMatchesOracle(t *testing.T) {
	tl := newTally()
	for _, c := range corpus3(t, 300) {
		a, b, cc, d := c.p[0], c.p[1], c.p[2], c.p[3]
		want := oracleOrient3D(a, b, cc, d)
		if got := sign(Orient3D(a, b, cc, d)); got != want {
			t.Fatalf("Orient3D(%v, %v, %v, %v) sign %d, exact %d (%s)", a, b, cc, d, got, want, c.kind)
		}
		if got := sign(Orient3DExact(a, b, cc, d)); got != want {
			t.Fatalf("Orient3DExact(%v, %v, %v, %v) sign %d, exact %d (%s)", a, b, cc, d, got, want, c.kind)
		}
		tl.cases[c.kind]++
		if sign(Orient3DFast(a, b, cc, d)) != want {
			tl.fastWrong[c.kind]++
		}
	}
	tl.log(t, "Orient3D")
	if tl.fastWrong[kNear] == 0 {
		t.Fatalf("near-degenerate corpus too easy: Fast never wrong (%v)", tl.fastWrong)
	}
}

func TestInSphereSignMatchesOracle(t *testing.T) {
	tl := newTally()
	for _, c := range corpus3(t, 200) {
		a, b, cc, d, e := c.p[0], c.p[1], c.p[2], c.p[3], c.p[4]
		want := oracleInSphere(a, b, cc, d, e)
		if got := sign(InSphere(a, b, cc, d, e)); got != want {
			t.Fatalf("InSphere(%v, %v, %v, %v, %v) sign %d, exact %d (%s)", a, b, cc, d, e, got, want, c.kind)
		}
		if got := sign(InSphereExact(a, b, cc, d, e)); got != want {
			t.Fatalf("InSphereExact(%v, %v, %v, %v, %v) sign %d, exact %d (%s)", a, b, cc, d, e, got, want, c.kind)
		}
		tl.cases[c.kind]++
		if sign(InSphereFast(a, b, cc, d, e)) != want {
			tl.fastWrong[c.kind]++
		}
	}
	tl.log(t, "InSphere")
	if tl.fastWrong[kNear] == 0 {
		t.Fatalf("near-degenerate corpus too easy: Fast never wrong (%v)", tl.fastWrong)
	}
}

func TestDocumentedSigns3D(t *testing.T) {
	a, b, c := [3]float64{0, 0, 0}, [3]float64{1, 0, 0}, [3]float64{0, 1, 0}
	below, above := [3]float64{0, 0, -1}, [3]float64{0, 0, 1}
	if Orient3D(a, b, c, below) <= 0 || Orient3D(a, b, c, above) >= 0 || Orient3D(a, b, c, [3]float64{5, 7, 0}) != 0 {
		t.Fatal("Orient3D: below must be positive, above negative, coplanar zero")
	}
	// A positively oriented tetrahedron on the unit sphere.
	p, q, r, s := [3]float64{1, 0, 0}, [3]float64{0, 1, 0}, [3]float64{-1, 0, 0}, [3]float64{0, 0, -1}
	if Orient3D(p, q, r, s) <= 0 {
		t.Fatal("test tetrahedron not positively oriented")
	}
	if InSphere(p, q, r, s, [3]float64{0, 0, 0}) <= 0 || InSphere(p, q, r, s, [3]float64{2, 2, 2}) >= 0 ||
		InSphere(p, q, r, s, [3]float64{0, 0, 1}) != 0 {
		t.Fatal("InSphere: inside must be positive, outside negative, cospherical zero")
	}
}

func TestNonFiniteYieldsNaN3D(t *testing.T) {
	o, x, y, z := [3]float64{}, [3]float64{1, 0, 0}, [3]float64{0, 1, 0}, [3]float64{0, 0, 1}
	for _, bad := range []float64{math.NaN(), math.Inf(1), math.Inf(-1)} {
		for i := 0; i < 3; i++ {
			p := [3]float64{}
			p[i] = bad
			for _, got := range []float64{
				Orient3D(p, x, y, z), Orient3D(o, x, y, p), Orient3DExact(x, p, y, z),
				InSphere(p, x, y, z, o), InSphere(o, x, y, z, p), InSphereExact(o, x, p, y, z),
			} {
				if !math.IsNaN(got) {
					t.Fatalf("non-finite input %v: got %v, want NaN", p, got)
				}
			}
		}
	}
}

func TestPlantedBreak3D(t *testing.T) {
	cs := corpus3(t, 200)
	run := func(name string, zero func(), restore func()) {
		zero()
		defer restore()
		wo, wi := 0, 0
		for _, c := range cs {
			if sign(Orient3D(c.p[0], c.p[1], c.p[2], c.p[3])) != oracleOrient3D(c.p[0], c.p[1], c.p[2], c.p[3]) {
				wo++
			}
			if sign(InSphere(c.p[0], c.p[1], c.p[2], c.p[3], c.p[4])) != oracleInSphere(c.p[0], c.p[1], c.p[2], c.p[3], c.p[4]) {
				wi++
			}
		}
		t.Logf("%s removed: Orient3D wrong %d, InSphere wrong %d of %d", name, wo, wi, len(cs))
		if wo == 0 || wi == 0 {
			t.Fatalf("planted break not detected: the corpus cannot see the %s", name)
		}
	}
	oa, ia, ob, ib := o3dErrBoundA, ispErrBoundA, o3dErrBoundB, ispErrBoundB
	run("stage A bound", func() { o3dErrBoundA, ispErrBoundA = 0, 0 }, func() { o3dErrBoundA, ispErrBoundA = oa, ia })
	run("stage B bound", func() { o3dErrBoundB, ispErrBoundB = 0, 0 }, func() { o3dErrBoundB, ispErrBoundB = ob, ib })
}
