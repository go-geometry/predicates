package bench

import (
	"testing"
	"time"

	"github.com/go-scikit/refbench"
	"github.com/go-scikit/refbench/rbtest"
)

// cases registers every benchmark row: each predicate on "random" records
// (decided by the floating-point filter) and on "degenerate" records (which run
// to the exact stage D), Go against predicates.c, timed by refbench.
//
//	cd bench
//	go test -run TestPerf .                          CI mode: allocations, table, denominator
//	go test -run TestPerf -cpu 1 -perf.run .         measure (refuses at a 1-minute load of 8 or more)
//	go test -run TestPerf -perf.record .             allocation counts only (a deferred report)
func cases() []refbench.Case {
	var out []refbench.Case
	for _, pr := range preds {
		for _, deg := range []bool{false, true} {
			kind := "random"
			if deg {
				kind = "degenerate"
			}
			pts := inputs(pr.per, deg, pr.repeat)
			pr := pr
			out = append(out, refbench.Case{
				Name:      pr.name + "/" + kind,
				Func:      pr.name,
				Reference: "predicates.c (Shewchuk, cgo, no FMA contraction)",
				Parity:    refbench.ParityExact, // the signs agree on every record: TestSignsAgreeWithPredicatesC
				Go:        goBody(pr.kind, pts),
				// The body cycles through the records: count every one.
				AllocIters: records,
				Ref: func(iters int) (time.Duration, error) {
					t0 := time.Now()
					sink = Loop(pr.kind, pts, records, iters)
					return time.Since(t0), nil
				},
			})
		}
	}
	return out
}

func TestPerf(t *testing.T) {
	rbtest.Perf(t, refbench.CIConfig{Root: "..", Cases: cases()}, refbench.Options{
		Repo:        "go-geometry/predicates",
		RefVersions: map[string]string{"predicates.c": "Shewchuk, public domain (copied verbatim in bench/shewchuk.c)", "cgo flags": "-O2 -w; no FMA contraction (pragmas in shewchuk.c)"},
	})
}
