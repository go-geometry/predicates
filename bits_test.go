package predicates

import (
	"hash/fnv"
	"math"
	"runtime"
	"testing"
)

// goldenBits is the FNV-64a hash of the bits of every value the predicates
// return on bitsCorpus, recorded on darwin/arm64 (where Go may fuse
// multiply-adds). CI runs this test natively on amd64 and arm64 and under
// emulation on riscv64, loong64, ppc64le, s390x, 386 and arm: the same hash
// proves that no product was fused into an FMA anywhere and that every stage
// rounds identically. The corpus uses only +, -, *, Ldexp and Nextafter, all
// exact or correctly rounded on every architecture (no trigonometry, whose
// assembly implementations differ).
const goldenBits = 0x7829fdcfc9d3738f

// bitsCorpus feeds every predicate, every variant, inputs that reach every
// stage: random, near-degenerate, exactly degenerate and out of range.
func bitsCorpus(emit func(float64)) {
	r := newRand(11)
	p2 := func(k int) [2]float64 {
		switch k % 4 {
		case 0:
			return randPt2(r, 1)
		case 1:
			return [2]float64{nudge(0.5, r.IntN(64)), nudge(0.5, r.IntN(64))}
		case 2:
			return [2]float64{randExp(r, -30, 60), randExp(r, -30, 60)}
		}
		return [2]float64{randExp(r, -700, 700), randExp(r, -700, 700)}
	}
	p3 := func(k int) [3]float64 {
		a, b := p2(k), p2(k)
		return [3]float64{a[0], a[1], b[0]}
	}
	for i := 0; i < 400; i++ {
		k := i
		a, b, c, d := p2(k), p2(k), p2(k), p2(k)
		if i%5 == 0 { // collinear / repeated
			c = [2]float64{a[0] + float64(3*(b[0]-a[0])), a[1] + float64(3*(b[1]-a[1]))}
			d = b
		}
		emit(Orient2D(a, b, c))
		emit(Orient2DFast(a, b, c))
		emit(Orient2DExact(a, b, c))
		emit(InCircle(a, b, c, d))
		emit(InCircleFast(a, b, c, d))
		emit(InCircleExact(a, b, c, d))
		e3, f3, g3, h3, j3 := p3(k), p3(k), p3(k), p3(k), p3(k)
		if i%5 == 0 {
			h3 = [3]float64{e3[0] + f3[0] - g3[0], e3[1] + f3[1] - g3[1], e3[2] + f3[2] - g3[2]}
			j3 = f3
		}
		emit(Orient3D(e3, f3, g3, h3))
		emit(Orient3DFast(e3, f3, g3, h3))
		emit(Orient3DExact(e3, f3, g3, h3))
		emit(InSphere(e3, f3, g3, h3, j3))
		emit(InSphereFast(e3, f3, g3, h3, j3))
		emit(InSphereExact(e3, f3, g3, h3, j3))
		o, _ := OrientN([][]float64{e3[:], f3[:], g3[:], h3[:]})
		s, _ := InSphereN([][]float64{e3[:], f3[:], g3[:], h3[:], j3[:]})
		emit(o)
		emit(s)
	}
}

func TestBitsIdenticalAcrossArchitectures(t *testing.T) {
	h := fnv.New64a()
	var buf [8]byte
	n := 0
	bitsCorpus(func(x float64) {
		u := math.Float64bits(x)
		for i := range buf {
			buf[i] = byte(u >> (8 * i))
		}
		h.Write(buf[:])
		n++
	})
	got := h.Sum64()
	t.Logf("%s/%s: %d values, hash %#x", runtime.GOOS, runtime.GOARCH, n, got)
	if got != goldenBits {
		t.Fatalf("hash %#x, golden %#x: a result differs in its bits from darwin/arm64", got, uint64(goldenBits))
	}
}
