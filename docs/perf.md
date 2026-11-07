# Speed against predicates.c

The Go predicates are measured against Shewchuk's public-domain `predicates.c`
(copied verbatim in `bench/`, compiled without FMA contraction, which its
error-free transformations forbid) by [refbench](https://github.com/go-scikit/refbench),
the shared tool of the go-scikit performance campaign. The numbers are in
[perf-table.md](perf-table.md), generated from [perf.json](perf.json); do not
edit the table.

## Method

- One core on both sides (`-cpu 1`, `GOMAXPROCS=1` set by refbench; the C loop is
  single-threaded).
- Each side is a loop calibrated to the same minimum duration, timed inside its
  own process, over the same 256 records (`inputs` in `bench/cref_test.go`):
  *random* points in [-1, 1), which the floating-point filter decides, and
  *degenerate* records, one point repeated under inexact translations of mixed
  exponents, which run to the exact stage D.
- Rounds alternate: Go then C, C then Go; the 1-minute load is read before and
  after every round and stored. The tool refuses to measure at a load of 8 or
  more. The previous harness ran all the C rounds, then all the Go rounds, with
  no gate.
- The verdict is on the median of the rounds with the quartiles kept: *slower*
  when the Go first quartile is above the C third quartile.
- `TestSignsAgreeWithPredicatesC` checks that both sides return the same sign on
  every benchmark record, so the cases time the same work.

## Status

The allocation counts per call are recorded in `perf.json` and checked on every
pull request (`cd bench && go test -run TestPerf .`: allocations, table,
and the list of public functions with their cases or [exemptions](perf-exempt.md)).

**Timing is deferred**: it has not been run with the tool yet, because the
host's 1-minute load was above the gate (29 when the tool was added). The
numbers measured earlier with the old harness, at a load of 7.3 to 8.8 and so
indicative only, are in
[go-geometry/predicates#4](https://github.com/go-geometry/predicates/issues/4):
Go 1.6 to 3.3 times slower than C on random records, 1.4 to 3.6 times on the
degenerate records that reach stage D. To measure:

```text
cd bench
go test -run TestPerf -cpu 1 -perf.run .
```

which rewrites `perf.json` and `perf-table.md` only if the gate opens, and
`refbench issues docs/perf.json` then drafts one issue per slower row.
