package predicates

import (
	"errors"
	"math"
)

// MaxDim is the largest dimension OrientN and InSphereN accept. It bounds
// the cost of the exact evaluation: with coordinates spread over the whole
// float64 exponent range, the integers of the elimination grow to
// d·2100 bits, and d = 16 already takes about 40 ms (Apple M4, Go 1.27);
// d = 32 took 0.9 s and d = 64 took 56 s, an unacceptable worst case for a
// predicate fed untrusted input.
const MaxDim = 16

// Errors returned by OrientN and InSphereN.
var (
	ErrShape     = errors.New("predicates: wrong number of points or coordinates")
	ErrDimension = errors.New("predicates: dimension out of range [1, MaxDim]")
	ErrNotFinite = errors.New("predicates: coordinate is NaN or infinite")
)

// expansionMaxN is the largest matrix order evaluated with floating-point
// expansions; larger matrices use exact big-integer elimination, whose
// O(n³) cost beats the O(2^n·n) Laplace expansion there.
const expansionMaxN = 8

// OrientN reports the orientation of d+1 points in d dimensions: the sign of
//
//	| p0[0] ... p0[d-1] 1 |
//	|  ...                |
//	| pd[0] ... pd[d-1] 1 |
//
// which equals Orient2D for d = 2 and Orient3D for d = 3. The sign is exact;
// the magnitude approximates the determinant (d! times the signed volume of
// the simplex). points must hold d+1 points of d finite coordinates each,
// 1 <= d <= MaxDim.
func OrientN(points [][]float64) (float64, error) {
	d, err := checkPoints(points, 1)
	if err != nil {
		return 0, err
	}
	n := d + 1
	if n <= expansionMaxN && expansionSafe(points, d) {
		ar := &arena{}
		m := make([][][]float64, n)
		for i, p := range points {
			row := make([][]float64, n)
			for j, x := range p {
				row[j] = one(ar, x)
			}
			row[d] = one(ar, 1)
			m[i] = row
		}
		det := detExp(ar, m)
		return det[len(det)-1], nil
	}
	s := scaleOfPoints(points)
	m := make([][]*bigInt, n)
	for i, p := range points {
		row := make([]*bigInt, n)
		for j, x := range p {
			row[j] = toBig(x, s)
		}
		row[d] = new(bigInt).SetInt64(1) // the column of ones, unscaled
		m[i] = row
	}
	return bigDet(m, d*s), nil
}

// InSphereN reports whether the last of d+2 points in d dimensions lies
// inside the sphere through the first d+1: the sign of
//
//	| p0[0]   ... p0[d-1]   |p0|²   1 |
//	|  ...                            |
//	| pd+1[0] ... pd+1[d-1] |pd+1|² 1 |
//
// which equals InCircle for d = 2 and InSphere for d = 3. It is the
// orientation of the points lifted onto the paraboloid x ↦ (x, |x|²), so the
// coefficient of the last point's |x|² is minus OrientN of the first d+1
// points: the result is positive if the last point lies inside the sphere
// and OrientN(points[:d+1]) > 0, negative if it lies outside, zero if the
// points are cospherical; the sign is reversed when OrientN(points[:d+1]) < 0.
// The sign is exact. points must hold d+2 points of d finite coordinates
// each, 1 <= d <= MaxDim.
func InSphereN(points [][]float64) (float64, error) {
	d, err := checkPoints(points, 2)
	if err != nil {
		return 0, err
	}
	n := d + 2
	if n <= expansionMaxN && expansionSafe(points, d+2) {
		ar := &arena{}
		m := make([][][]float64, n)
		for i, p := range points {
			row := make([][]float64, n)
			for j, x := range p {
				row[j] = one(ar, x)
			}
			row[d] = squareSumExp(ar, row[:d]...)
			row[d+1] = one(ar, 1)
			m[i] = row
		}
		det := detExp(ar, m)
		return det[len(det)-1], nil
	}
	s := scaleOfPoints(points)
	m := make([][]*bigInt, n)
	for i, p := range points {
		row := make([]*bigInt, n)
		for j, x := range p {
			row[j] = toBig(x, s)
		}
		row[d] = bigSquareSum(row[:d]...)
		row[d+1] = new(bigInt).SetInt64(1)
		m[i] = row
	}
	return bigDet(m, (d+2)*s), nil
}

// scaleOfPoints is scaleOf over every coordinate of points.
func scaleOfPoints(points [][]float64) int {
	s := scaleOf(points[0]...)
	for _, p := range points[1:] {
		s = max(s, scaleOf(p...))
	}
	return s
}

// checkPoints validates d+extra points of d finite coordinates and returns d.
func checkPoints(points [][]float64, extra int) (int, error) {
	if len(points) == 0 {
		return 0, ErrShape
	}
	d := len(points[0])
	if d < 1 || d > MaxDim {
		return 0, ErrDimension
	}
	if len(points) != d+extra {
		return 0, ErrShape
	}
	for _, p := range points {
		if len(p) != d {
			return 0, ErrShape
		}
		if !finite(p...) {
			return 0, ErrNotFinite
		}
	}
	return d, nil
}

// expansionSafe reports whether a determinant of degree k in these
// coordinates can be evaluated with floating-point expansions without
// overflow or underflow: every nonzero coordinate x satisfies
// 2^(52 - 1000/k) <= |x| <= 2^(900/k), so that every product of k
// coordinates is a normal number (a multiple of 2^-1000 at worst) and the
// largest sum stays far below 2^1024 even after Split's 2^27 factor and the
// at most 8! terms of the Laplace expansion.
func expansionSafe(points [][]float64, k int) bool {
	lo := math.Ldexp(1, 52-1000/k)
	hi := math.Ldexp(1, 900/k)
	for _, p := range points {
		for _, x := range p {
			if !inRange(x, lo, hi) {
				return false
			}
		}
	}
	return true
}
