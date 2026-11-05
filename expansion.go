package predicates

// Floating-point expansion arithmetic (Shewchuk 1997, section 2).
//
// An expansion is a slice of float64 components, sorted by increasing
// magnitude, pairwise nonoverlapping, whose exact sum is the value it
// represents. Zero components are eliminated, except that the expansion of
// zero is the single component 0.
//
// ⛔ No fused multiply-add. The Go compiler may fuse x*y + z (and may do so
// across statements) on arm64, ppc64le, s390x and riscv64; a fused product
// is not rounded, which breaks every error-free transformation below. Every
// product is therefore wrapped in an explicit float64(...) conversion, which
// the Go specification defines as a rounding that prevents fusion.

// epsilon is 2^-53, half an ulp of 1: the largest relative rounding error of
// a float64 operation under round-to-nearest. splitter is 2^27 + 1, the
// constant that splits a 53-bit significand into two 26-bit halves. The
// paper's exactinit() derives both by a loop; for IEEE 754 binary64 they are
// these constants.
const (
	epsilon  = 1.0 / (1 << 53)
	splitter = 1<<27 + 1
)

// Error-bound coefficients (predicates.c exactinit(); paper section 4).
// They are variables, not constants, only so that the planted-break tests can
// zero them and prove the suite notices.
var (
	resultErrBound = (3.0 + 8.0*epsilon) * epsilon
	ccwErrBoundA   = (3.0 + 16.0*epsilon) * epsilon
	ccwErrBoundB   = (2.0 + 12.0*epsilon) * epsilon
	ccwErrBoundC   = (9.0 + 64.0*epsilon) * epsilon * epsilon
	o3dErrBoundA   = (7.0 + 56.0*epsilon) * epsilon
	o3dErrBoundB   = (3.0 + 28.0*epsilon) * epsilon
	o3dErrBoundC   = (26.0 + 288.0*epsilon) * epsilon * epsilon
	iccErrBoundA   = (10.0 + 96.0*epsilon) * epsilon
	iccErrBoundB   = (4.0 + 48.0*epsilon) * epsilon
	iccErrBoundC   = (44.0 + 576.0*epsilon) * epsilon * epsilon
	ispErrBoundA   = (16.0 + 224.0*epsilon) * epsilon
	ispErrBoundB   = (5.0 + 72.0*epsilon) * epsilon
	ispErrBoundC   = (71.0 + 1408.0*epsilon) * epsilon * epsilon
)

// fastTwoSum returns x = fl(a+b) and the roundoff y, with a+b = x+y exactly.
// It requires |a| >= |b| (or a == 0).
func fastTwoSum(a, b float64) (x, y float64) {
	x = a + b
	y = b - (x - a)
	return
}

// twoSum returns x = fl(a+b) and the roundoff y, with a+b = x+y exactly, for
// any a and b (Knuth's six-operation sum).
func twoSum(a, b float64) (x, y float64) {
	x = a + b
	bv := x - a
	av := x - bv
	y = (a - av) + (b - bv)
	return
}

// twoDiff returns x = fl(a-b) and the roundoff y, with a-b = x+y exactly.
func twoDiff(a, b float64) (x, y float64) {
	x = a - b
	return x, twoDiffTail(a, b, x)
}

// twoDiffTail returns the roundoff of x = fl(a-b).
func twoDiffTail(a, b, x float64) float64 {
	bv := a - x
	av := x + bv
	return (a - av) + (bv - b)
}

// split returns hi and lo, each with at most 26 significant bits, such that
// a = hi + lo (Dekker's splitting; the paper's Theorem 17).
func split(a float64) (hi, lo float64) {
	c := float64(splitter * a)
	big := c - a
	hi = c - big
	lo = a - hi
	return
}

// twoProduct returns x = fl(a*b) and the roundoff y, with a*b = x+y exactly.
func twoProduct(a, b float64) (x, y float64) {
	x = float64(a * b)
	ahi, alo := split(a)
	bhi, blo := split(b)
	return x, productTail(x, ahi, alo, bhi, blo)
}

// productTail is the roundoff of x = fl(a*b) given the splits of a and b.
func productTail(x, ahi, alo, bhi, blo float64) float64 {
	e1 := x - float64(ahi*bhi)
	e2 := e1 - float64(alo*bhi)
	e3 := e2 - float64(ahi*blo)
	return float64(alo*blo) - e3
}

// twoTwoDiff returns the four-component expansion of (a1+a0) - (b1+b0), where
// each operand is a two-component expansion (the paper's Two-Two-Diff).
func twoTwoDiff(a1, a0, b1, b0 float64) [4]float64 {
	i, x0 := twoDiff(a0, b0)
	j, r0 := twoSum(a1, i)
	k, x1 := twoDiff(r0, b1)
	x3, x2 := twoSum(j, k)
	return [4]float64{x0, x1, x2, x3}
}

// cross2 returns the exact four-component expansion of a*d - b*c.
func cross2(a, d, b, c float64) [4]float64 {
	l1, l0 := twoProduct(a, d)
	r1, r0 := twoProduct(b, c)
	return twoTwoDiff(l1, l0, r1, r0)
}

// sumExp sets h = e + f and returns h. e and f are nonempty expansions; h must
// not alias either and its capacity should be at least len(e)+len(f) (append
// grows it otherwise). This is the zero-eliminating fast expansion sum (the
// paper's Theorem 13; it relies on round-to-even tie breaking).
func sumExp(h, e, f []float64) []float64 {
	h = h[:0]
	ei, fi := 0, 0
	enow, fnow := e[0], f[0]
	var q float64
	// Take the smaller-magnitude component first.
	if (fnow > enow) == (fnow > -enow) {
		q = enow
		ei++
	} else {
		q = fnow
		fi++
	}
	if ei < len(e) && fi < len(f) {
		enow, fnow = e[ei], f[fi]
		var hh float64
		if (fnow > enow) == (fnow > -enow) {
			q, hh = fastTwoSum(enow, q)
			ei++
		} else {
			q, hh = fastTwoSum(fnow, q)
			fi++
		}
		if hh != 0 {
			h = append(h, hh)
		}
		for ei < len(e) && fi < len(f) {
			enow, fnow = e[ei], f[fi]
			if (fnow > enow) == (fnow > -enow) {
				q, hh = twoSum(q, enow)
				ei++
			} else {
				q, hh = twoSum(q, fnow)
				fi++
			}
			if hh != 0 {
				h = append(h, hh)
			}
		}
	}
	for ; ei < len(e); ei++ {
		var hh float64
		q, hh = twoSum(q, e[ei])
		if hh != 0 {
			h = append(h, hh)
		}
	}
	for ; fi < len(f); fi++ {
		var hh float64
		q, hh = twoSum(q, f[fi])
		if hh != 0 {
			h = append(h, hh)
		}
	}
	if q != 0 || len(h) == 0 {
		h = append(h, q)
	}
	return h
}

// scaleExp sets h = b*e and returns h. h must not alias e; its capacity should
// be at least 2*len(e). Zero-eliminating expansion scaling (Theorem 19).
func scaleExp(h, e []float64, b float64) []float64 {
	h = h[:0]
	bhi, blo := split(b)
	q := float64(e[0] * b)
	ahi, alo := split(e[0])
	hh := productTail(q, ahi, alo, bhi, blo)
	if hh != 0 {
		h = append(h, hh)
	}
	for _, en := range e[1:] {
		p1 := float64(en * b)
		ahi, alo := split(en)
		p0 := productTail(p1, ahi, alo, bhi, blo)
		var sum float64
		sum, hh = twoSum(q, p0)
		if hh != 0 {
			h = append(h, hh)
		}
		q, hh = fastTwoSum(p1, sum)
		if hh != 0 {
			h = append(h, hh)
		}
	}
	if q != 0 || len(h) == 0 {
		h = append(h, q)
	}
	return h
}

// arena hands out expansion buffers from large chunks, so that an exact
// evaluation makes a handful of allocations instead of one per operation.
type arena struct{ buf []float64 }

// get returns an empty slice with capacity n.
func (a *arena) get(n int) []float64 {
	if cap(a.buf)-len(a.buf) < n {
		a.buf = make([]float64, 0, max(2*cap(a.buf), n, 512))
	}
	l := len(a.buf)
	a.buf = a.buf[:l+n]
	return a.buf[l : l : l+n]
}

// mulExp returns the exact product e*f: the sum over the components of f of
// e scaled by each. The smaller operand is the one iterated over.
func mulExp(ar *arena, e, f []float64) []float64 {
	if len(f) > len(e) {
		e, f = f, e
	}
	acc := scaleExp(ar.get(2*len(e)), e, f[0])
	if len(f) == 1 {
		return acc
	}
	tmp := make([]float64, 0, 2*len(e))
	for _, c := range f[1:] {
		tmp = scaleExp(tmp, e, c)
		acc = sumExp(ar.get(len(acc)+len(tmp)), acc, tmp)
	}
	return acc
}

// negExp returns -e.
func negExp(ar *arena, e []float64) []float64 {
	h := ar.get(len(e))[:len(e)]
	for i, c := range e {
		h[i] = -c
	}
	return h
}

// addExp returns e+f.
func addExp(ar *arena, e, f []float64) []float64 {
	return sumExp(ar.get(len(e)+len(f)), e, f)
}

// estimate returns the floating-point sum of the components of e, an
// approximation of its value whose sign is the sign of e.
func estimate(e []float64) float64 {
	q := e[0]
	for _, c := range e[1:] {
		q += c
	}
	return q
}

// abs is |a|, without a call to math.Abs so that the filters stay inlinable
// in a predictable way. -0 maps to -0, which compares equal to 0.
func abs(a float64) float64 {
	if a >= 0 {
		return a
	}
	return -a
}
