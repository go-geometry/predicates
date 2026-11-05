package predicates

import (
	"errors"
	"math"
	"math/big"
	"math/rand/v2"
	"testing"
)

func oracleOrientN(pts [][]float64) int {
	m := make([][]*big.Rat, len(pts))
	for i, p := range pts {
		row := make([]*big.Rat, 0, len(p)+1)
		for _, x := range p {
			row = append(row, rat(x))
		}
		m[i] = append(row, big.NewRat(1, 1))
	}
	return ratDet(m).Sign()
}

func oracleInSphereN(pts [][]float64) int {
	m := make([][]*big.Rat, len(pts))
	for i, p := range pts {
		row := make([]*big.Rat, 0, len(p)+2)
		lift := new(big.Rat)
		for _, x := range p {
			v := rat(x)
			row = append(row, v)
			lift.Add(lift, new(big.Rat).Mul(v, v))
		}
		m[i] = append(row, lift, big.NewRat(1, 1))
	}
	return ratDet(m).Sign()
}

// corpusN returns point sets of d+extra points in d dimensions.
func corpusN(r *rand.Rand, d, extra, n int) (sets [][][]float64, kinds []kind) {
	add := func(k kind, s [][]float64) {
		sets = append(sets, s)
		kinds = append(kinds, k)
	}
	np := d + extra
	randSet := func(gen func() float64) [][]float64 {
		s := make([][]float64, np)
		for i := range s {
			s[i] = make([]float64, d)
			for j := range s[i] {
				s[i][j] = gen()
			}
		}
		return s
	}
	for i := 0; i < n; i++ {
		add(kRandom, randSet(func() float64 { return 2*r.Float64() - 1 }))
		// Near-degenerate: the last point rounded onto the affine hull of d
		// others (extra = 1) or onto the sphere through them (extra = 2,
		// approximated by a point on a sphere around a common centre).
		s := randSet(func() float64 { return 2*r.Float64() - 1 })
		if extra == 1 {
			for j := 0; j < d; j++ {
				x := s[0][j]
				for k := 1; k < d; k++ {
					x += 0.3 * float64(k) * (s[k][j] - s[0][j])
				}
				s[d][j] = nudge(x, r.IntN(3)-1)
			}
		} else {
			for _, p := range s {
				norm := 0.0
				for _, x := range p {
					norm += x * x
				}
				norm = math.Sqrt(norm)
				for j := range p {
					p[j] = 0.25 + p[j]/norm
				}
			}
			s[np-1][0] = nudge(s[np-1][0], r.IntN(5)-2)
		}
		add(kNear, s)
		// Exactly degenerate: integer points on the hyperplane x_0 = c, or a
		// repeated point.
		c := float64(r.IntN(1 << 20))
		s = randSet(func() float64 { return float64(r.IntN(2001) - 1000) })
		for _, p := range s {
			p[0] = c
		}
		add(kExact, s)
		s = randSet(func() float64 { return randExp(r, -20, 50) })
		copy(s[np-1], s[0])
		add(kExact, s)
		add(kHuge, randSet(func() float64 { return randExp(r, -600, 600) }))
	}
	return
}

func TestOrientNSignMatchesOracle(t *testing.T) {
	r := newRand(5)
	for d := 1; d <= 9; d++ {
		sets, kinds := corpusN(r, d, 1, 40)
		for i, s := range sets {
			got, err := OrientN(s)
			if err != nil {
				t.Fatal(err)
			}
			if want := oracleOrientN(s); sign(got) != want {
				t.Fatalf("d=%d OrientN(%v) sign %d, exact %d (%s)", d, s, sign(got), want, kinds[i])
			}
		}
	}
}

func TestInSphereNSignMatchesOracle(t *testing.T) {
	r := newRand(6)
	for d := 1; d <= 7; d++ {
		sets, kinds := corpusN(r, d, 2, 30)
		for i, s := range sets {
			got, err := InSphereN(s)
			if err != nil {
				t.Fatal(err)
			}
			if want := oracleInSphereN(s); sign(got) != want {
				t.Fatalf("d=%d InSphereN(%v) sign %d, exact %d (%s)", d, s, sign(got), want, kinds[i])
			}
		}
	}
}

// In 2-D and 3-D the n-D predicates are the fixed-dimension ones.
func TestNDimAgreesWithFixed(t *testing.T) {
	for _, c := range corpus2(t, 100) {
		a, b, cc, d := c.p[0], c.p[1], c.p[2], c.p[3]
		o, _ := OrientN([][]float64{a[:], b[:], cc[:]})
		i, _ := InSphereN([][]float64{a[:], b[:], cc[:], d[:]})
		if sign(o) != sign(Orient2D(a, b, cc)) || sign(i) != sign(InCircle(a, b, cc, d)) {
			t.Fatalf("2-D disagreement on %v (%s)", c.p, c.kind)
		}
	}
	for _, c := range corpus3(t, 60) {
		a, b, cc, d, e := c.p[0], c.p[1], c.p[2], c.p[3], c.p[4]
		o, _ := OrientN([][]float64{a[:], b[:], cc[:], d[:]})
		i, _ := InSphereN([][]float64{a[:], b[:], cc[:], d[:], e[:]})
		if sign(o) != sign(Orient3D(a, b, cc, d)) || sign(i) != sign(InSphere(a, b, cc, d, e)) {
			t.Fatalf("3-D disagreement on %v (%s)", c.p, c.kind)
		}
	}
}

// Large dimensions go through big-integer elimination; a regular simplex
// with the identity matrix has a known sign.
func TestNDimLarge(t *testing.T) {
	for _, d := range []int{8, 12, MaxDim} {
		pts := make([][]float64, d+1)
		for i := range pts {
			pts[i] = make([]float64, d)
			if i < d {
				pts[i][i] = 1
			}
		}
		// det [[I, 1], [0, 1]] = 1.
		if got, err := OrientN(pts); err != nil || got != 1 {
			t.Fatalf("d=%d OrientN(identity simplex) = %v, %v; want 1", d, got, err)
		}
		// The origin lies inside the sphere through e_1..e_d and -e_1
		// ... use the sphere through the simplex and test its centroid.
		sp := append(pts[:d+1:d+1], make([]float64, d))
		for j := range sp[d+1] {
			sp[d+1][j] = 1 / float64(d+1)
		}
		o, _ := OrientN(pts)
		got, err := InSphereN(sp)
		if err != nil || sign(got) != sign(o) {
			t.Fatalf("d=%d centroid InSphereN = %v, %v; want the sign of the orientation %v", d, got, err, o)
		}
	}
}

func TestNDimErrors(t *testing.T) {
	inf := math.Inf(1)
	cases := []struct {
		pts  [][]float64
		want error
	}{
		{nil, ErrShape},
		{[][]float64{{}}, ErrDimension},
		{[][]float64{make([]float64, MaxDim+1)}, ErrDimension},
		{[][]float64{{0, 0}, {1, 0}}, ErrShape},
		{[][]float64{{0, 0}, {1, 0}, {1}}, ErrShape},
		{[][]float64{{0, 0}, {1, inf}, {0, 1}}, ErrNotFinite},
		{[][]float64{{0, 0}, {1, math.NaN()}, {0, 1}}, ErrNotFinite},
	}
	for _, c := range cases {
		if _, err := OrientN(c.pts); !errors.Is(err, c.want) {
			t.Errorf("OrientN(%v) error %v, want %v", c.pts, err, c.want)
		}
	}
	if _, err := InSphereN([][]float64{{0}, {1}}); !errors.Is(err, ErrShape) {
		t.Errorf("InSphereN with d+1 points: error %v, want ErrShape", err)
	}
}
