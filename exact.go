package predicates

import (
	"math"
	"math/big"
)

// Input range of the floating-point paths.
//
// The error bounds of the filters and the error-free transformations of the
// expansions hold only when no intermediate result overflows or underflows
// (paper, section 5, "Caveats"). For a predicate whose determinant has degree
// k in the coordinates, every nonzero coordinate must lie in [loK, hiK]:
//
//   - a nonzero coordinate >= 2^-L is a multiple of 2^-(L+52), so every
//     difference, tail, product and expansion component is a multiple of
//     2^-k(L+52); with k(L+52) + 106 <= 1022 even the stage-C error bound,
//     which carries an epsilon^2 factor, is a normal number;
//   - with |x| <= 2^H, every intermediate value stays below 2^(kH+40),
//     including Split's 2^27 factor, and kH + 40 < 1024.
//
// Coordinates outside the range, or non-finite, go to an exact big-integer
// evaluation (or yield NaN).
const (
	lo2, hi2 = 0x1p-400, 0x1p480 // degree 2: Orient2D
	lo3, hi3 = 0x1p-250, 0x1p320 // degree 3: Orient3D
	lo4, hi4 = 0x1p-170, 0x1p240 // degree 4: InCircle
	lo5, hi5 = 0x1p-130, 0x1p190 // degree 5: InSphere
)

var nan = math.NaN()

// inRange reports whether x is zero or lo <= |x| <= hi. NaN and ±Inf are not.
func inRange(x, lo, hi float64) bool {
	x = abs(x)
	return x == 0 || (x >= lo && x <= hi)
}

func inRange2(a, b, c [2]float64, lo, hi float64) bool {
	return inRange(a[0], lo, hi) && inRange(a[1], lo, hi) &&
		inRange(b[0], lo, hi) && inRange(b[1], lo, hi) &&
		inRange(c[0], lo, hi) && inRange(c[1], lo, hi)
}

func finite(xs ...float64) bool {
	for _, x := range xs {
		if math.IsNaN(x) || math.IsInf(x, 0) {
			return false
		}
	}
	return true
}

// diffExp returns a-b as an exact expansion of one or two components.
func diffExp(a, b float64) []float64 {
	x, y := twoDiff(a, b)
	if y == 0 {
		return []float64{x}
	}
	return []float64{y, x}
}

// squareSumExp returns the exact expansion of the sum of the squares of the
// given expansions.
func squareSumExp(es ...[]float64) []float64 {
	acc := mulExp(es[0], es[0])
	for _, e := range es[1:] {
		acc = addExp(acc, mulExp(e, e))
	}
	return acc
}

// detExp returns the exact determinant of the n×n matrix of expansions m
// (row-major), by Laplace expansion along successive rows over column
// subsets: minor[mask] is the determinant of the first popcount(mask) rows
// restricted to the columns in mask. It costs O(2^n · n) expansion products,
// which is why OrientN and InSphereN use it only for small n.
func detExp(m [][][]float64) []float64 {
	n := len(m)
	minor := make([][]float64, 1<<n)
	for j := 0; j < n; j++ {
		minor[1<<j] = m[0][j]
	}
	for mask := 1; mask < 1<<n; mask++ {
		r := popcount(mask) - 1 // the row this mask's minor ends with
		if r == 0 {
			continue
		}
		var acc []float64
		pos := 0
		for j := 0; j < n; j++ {
			if mask&(1<<j) == 0 {
				continue
			}
			t := mulExp(m[r][j], minor[mask&^(1<<j)])
			if (r+pos)%2 == 1 {
				t = negExp(t)
			}
			if acc == nil {
				acc = t
			} else {
				acc = addExp(acc, t)
			}
			pos++
		}
		minor[mask] = acc
	}
	return minor[1<<n-1]
}

func popcount(x int) int {
	c := 0
	for ; x != 0; x &= x - 1 {
		c++
	}
	return c
}

// bigInt is the integer type of the exact fallback.
type bigInt = big.Int

// scaleShift is the power of two that turns every finite float64 into an
// integer: the smallest subnormal is 2^-1074.
const scaleShift = 1074

// toBig returns x·2^1074 as an exact integer. x must be finite.
func toBig(x float64) *big.Int {
	if x == 0 {
		return new(big.Int)
	}
	frac, exp := math.Frexp(x) // x = frac·2^exp, 0.5 <= |frac| < 1
	m := int64(frac * (1 << 53))
	shift := exp - 53 + scaleShift
	z := big.NewInt(m)
	if shift >= 0 {
		return z.Lsh(z, uint(shift))
	}
	// Only subnormals land here, and their low bits are zero.
	return z.Rsh(z, uint(-shift))
}

// bigDiff returns (a-b)·2^1074 exactly.
func bigDiff(a, b float64) *big.Int {
	return new(big.Int).Sub(toBig(a), toBig(b))
}

// bigSquareSum returns the sum of the squares of xs.
func bigSquareSum(xs ...*big.Int) *big.Int {
	s := new(big.Int)
	t := new(big.Int)
	for _, x := range xs {
		s.Add(s, t.Mul(x, x))
	}
	return s
}

// bigDet returns the determinant of the integer matrix m, computed exactly by
// fraction-free Gaussian elimination (Bareiss), converted to float64 after
// dividing by 2^(1074·degree). The result has the exact sign: a nonzero
// determinant too small for float64 returns ±the smallest subnormal, one too
// large returns ±Inf. m is overwritten.
func bigDet(m [][]*big.Int, degree int) float64 {
	n := len(m)
	sign := 1
	prev := big.NewInt(1)
	t := new(big.Int)
	for k := 0; k < n-1; k++ {
		if m[k][k].Sign() == 0 {
			p := -1
			for i := k + 1; i < n; i++ {
				if m[i][k].Sign() != 0 {
					p = i
					break
				}
			}
			if p < 0 {
				return 0
			}
			m[k], m[p] = m[p], m[k]
			sign = -sign
		}
		for i := k + 1; i < n; i++ {
			for j := k + 1; j < n; j++ {
				// m[i][j] = (m[k][k]·m[i][j] - m[i][k]·m[k][j]) / prev, exact.
				v := new(big.Int).Mul(m[k][k], m[i][j])
				v.Sub(v, t.Mul(m[i][k], m[k][j]))
				m[i][j] = v.Quo(v, prev)
			}
		}
		prev = m[k][k]
	}
	d := m[n-1][n-1]
	if sign < 0 {
		d = new(big.Int).Neg(d)
	}
	return bigToFloat(d, degree)
}

// bigToFloat returns d·2^(-1074·degree) rounded to float64, keeping the sign
// of d even when the magnitude underflows.
func bigToFloat(d *big.Int, degree int) float64 {
	if d.Sign() == 0 {
		return 0
	}
	f := new(big.Float).SetInt(d)
	f.SetMantExp(f, -scaleShift*degree)
	v, _ := f.Float64()
	if v == 0 {
		return math.Copysign(math.SmallestNonzeroFloat64, float64(d.Sign()))
	}
	return v
}
