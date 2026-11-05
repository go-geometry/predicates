package predicates

import (
	"math"
	"math/big"
	"math/bits"
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

// The hot paths check the range without branches on the bit patterns: for
// x >= 0, Float64bits is monotonic, NaN and ±Inf lie above every finite
// value, and u-1 maps the bits of zero to the largest uint64, so the minimum
// of u-1 over the coordinates ignores zeros. A set of coordinates is in range
// iff max(u) <= bits(hi) and min(u-1) >= bits(lo)-1.

const signMask = 1 << 63

// rangeAcc folds one coordinate into the running maximum and minimum.
func rangeAcc(mx, mn uint64, x float64) (uint64, uint64) {
	u := math.Float64bits(x) &^ signMask
	return max(mx, u), min(mn, u-1)
}

func rangeOK(mx, mn uint64, lo, hi float64) bool {
	return mx <= math.Float64bits(hi) && mn >= math.Float64bits(lo)-1
}

func inRange2(a, b, c [2]float64, lo, hi float64) bool {
	mx, mn := uint64(0), uint64(math.MaxUint64)
	mx, mn = rangeAcc(mx, mn, a[0])
	mx, mn = rangeAcc(mx, mn, a[1])
	mx, mn = rangeAcc(mx, mn, b[0])
	mx, mn = rangeAcc(mx, mn, b[1])
	mx, mn = rangeAcc(mx, mn, c[0])
	mx, mn = rangeAcc(mx, mn, c[1])
	return rangeOK(mx, mn, lo, hi)
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
func diffExp(ar *arena, a, b float64) []float64 {
	x, y := twoDiff(a, b)
	h := ar.get(2)
	if y != 0 {
		h = append(h, y)
	}
	return append(h, x)
}

// one returns the single-component expansion x.
func one(ar *arena, x float64) []float64 {
	return append(ar.get(1), x)
}

// squareSumExp returns the exact expansion of the sum of the squares of the
// given expansions.
func squareSumExp(ar *arena, es ...[]float64) []float64 {
	acc := mulExp(ar, es[0], es[0])
	for _, e := range es[1:] {
		acc = addExp(ar, acc, mulExp(ar, e, e))
	}
	return acc
}

// detExp returns the exact determinant of the n×n matrix of expansions m
// (row-major), by Laplace expansion along successive rows over column
// subsets: minor[mask] is the determinant of the first popcount(mask) rows
// restricted to the columns in mask. It costs O(2^n · n) expansion products,
// which is why OrientN and InSphereN use it only for small n.
func detExp(ar *arena, m [][][]float64) []float64 {
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
			x := m[r][j]
			if (r+pos)%2 == 1 {
				x = negExp(ar, x) // negate the small factor, not the product
			}
			t := mulExp(ar, x, minor[mask&^(1<<j)])
			if acc == nil {
				acc = t
			} else {
				acc = addExp(ar, acc, t)
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

// mantExp returns the integer significand m and exponent e of a finite
// nonzero x = m·2^e, with m odd: e is the weight of x's lowest set bit.
func mantExp(x float64) (int64, int) {
	frac, exp := math.Frexp(x) // x = frac·2^exp, 0.5 <= |frac| < 1
	m := int64(frac * (1 << 53))
	tz := bits.TrailingZeros64(uint64(m))
	return m >> tz, exp - 53 + tz
}

// scaleOf returns the power of two 2^s that makes every coordinate an
// integer with the fewest bits: s = -min over nonzero x of the weight of x's
// lowest set bit. Scaling every coordinate by the same 2^s multiplies a
// determinant of degree k by 2^(k·s), which keeps its sign; the integers
// are as short as the input's exponent spread allows (a few dozen bits for
// coordinates of similar magnitude, not the 2^1074 of a fixed scale).
func scaleOf(xs ...float64) int {
	s := math.MinInt
	for _, x := range xs {
		if x != 0 {
			_, e := mantExp(x)
			s = max(s, -e)
		}
	}
	if s == math.MinInt {
		return 0
	}
	return s
}

// toBig returns x·2^s as an exact integer. x must be finite and x·2^s an
// integer (s >= scaleOf of a set containing x).
func toBig(x float64, s int) *big.Int {
	if x == 0 {
		return new(big.Int)
	}
	m, e := mantExp(x)
	z := big.NewInt(m)
	return z.Lsh(z, uint(e+s))
}

// bigDiff returns (a-b)·2^s exactly.
func bigDiff(a, b float64, s int) *big.Int {
	return new(big.Int).Sub(toBig(a, s), toBig(b, s))
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
// multiplying by 2^-scale (the entries were coordinates scaled by 2^s, so a
// determinant of degree k has scale k·s). The result has the exact sign: a
// nonzero determinant too small for float64 returns ±the smallest subnormal,
// one too large returns ±Inf. m is overwritten.
func bigDet(m [][]*big.Int, scale int) float64 {
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
	return bigToFloat(d, scale)
}

// bigToFloat returns d·2^-scale rounded to float64, keeping the sign of d
// even when the magnitude underflows.
func bigToFloat(d *big.Int, scale int) float64 {
	if d.Sign() == 0 {
		return 0
	}
	f := new(big.Float).SetInt(d)
	f.SetMantExp(f, -scale)
	v, _ := f.Float64()
	if v == 0 {
		return math.Copysign(math.SmallestNonzeroFloat64, float64(d.Sign()))
	}
	return v
}
