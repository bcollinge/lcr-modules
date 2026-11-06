#!/usr/bin/env python3
"""Call the capture kit of one library from its bedcov over the kit panel.

Usage: call_kit.py <bedcov.tsv> <sample_id> <out_prefix> [--hi 0.35] [--lo 0.2]
       [--min-margin 0.3] [--min-core-depth 3] [--min-set-bp 20000] [--genome-floor 0.5]

Depth per set = summed base depth / summed bp; every set is divided by core depth.
For kit A: own_A = median over B of ratio(A--not--B), against_A = max over B of
ratio(B--not--A), each over sets with at least min_set_bp. Margin_A = own_A - against_A.
A kit passes when own_A >= hi, against_A <= lo and margin_A >= min_margin. A passing kit
is dropped when another passing kit B covers it (ratio(A--not--B) <= lo), which resolves
nested designs whose separating set is below min_set_bp.
Label, in order: too_shallow (core depth below min_core_depth), genome_like
(floor/core >= genome_floor), the kit when exactly one passes, mixed when several pass,
or when none passes but the best kit A has own_A >= hi and the kit B that sets
against_A has ratio(B--not--A) >= lo (call "A+B"); otherwise none_known.
The call's margin is that of the best kit by margin.

Writes <out_prefix>.capture_kit.tsv (one row per set) and
<out_prefix>.capture_kit.call.tsv (one row).
"""
import argparse
import statistics
from collections import defaultdict


def main():
    p = argparse.ArgumentParser()
    p.add_argument("bedcov")
    p.add_argument("sample_id")
    p.add_argument("out_prefix")
    p.add_argument("--hi", type=float, default=0.35)
    p.add_argument("--lo", type=float, default=0.2)
    p.add_argument("--min-margin", type=float, default=0.3)
    p.add_argument("--min-core-depth", type=float, default=3)
    p.add_argument("--min-set-bp", type=int, default=20000)
    p.add_argument("--genome-floor", type=float, default=0.5)
    a = p.parse_args()

    bases = defaultdict(int)
    bp = defaultdict(int)
    with open(a.bedcov) as f:
        for line in f:
            chrom, start, end, label, depth = line.rstrip("\n").split("\t")[:5]
            bp[label] += int(end) - int(start)
            bases[label] += int(depth)

    depth = {s: bases[s] / bp[s] for s in bp}
    core = depth.get("core", 0.0)
    ratio = {s: (depth[s] / core if core > 0 else float("nan")) for s in depth}

    pairs = {}
    for s in bp:
        if "--not--" in s and bp[s] >= a.min_set_bp:
            x, y = s.split("--not--")
            pairs[(x, y)] = ratio[s]
    kits = sorted({k for s in bp if "--not--" in s for k in s.split("--not--")})

    rows = []
    for k in kits:
        own = [r for (x, _), r in pairs.items() if x == k]
        against = [(r, x) for (x, y), r in pairs.items() if y == k]
        o = statistics.median(own) if own else float("nan")
        g, gk = max(against) if against else (float("nan"), "")
        rows.append((k, o, g, o - g, gk))

    passing = [k for k, o, g, m, _ in rows if o >= a.hi and g <= a.lo and m >= a.min_margin]
    passing = [k for k in passing
               if not any(b != k and pairs.get((k, b), 1) <= a.lo for b in passing)]
    best = max(rows, key=lambda r: r[3]) if rows else ("", float("nan"), float("nan"), float("nan"), "")
    floor_core = ratio.get("floor", float("nan"))

    if core < a.min_core_depth:
        label, call = "too_shallow", ""
    elif floor_core >= a.genome_floor:
        label, call = "genome_like", ""
    elif len(passing) == 1:
        label, call = "kit", passing[0]
    elif passing:
        label, call = "mixed", "+".join(passing)
    elif best[1] >= a.hi and pairs.get((best[4], best[0]), 0) >= a.lo:
        label, call = "mixed", best[0] + "+" + best[4]
    else:
        label, call = "none_known", ""

    with open(a.out_prefix + ".capture_kit.tsv", "w") as f:
        f.write("sample_id\tset\tbp\tdepth\tratio_to_core\n")
        for s in sorted(bp):
            f.write(f"{a.sample_id}\t{s}\t{bp[s]}\t{depth[s]:.3f}\t{ratio[s]:.4f}\n")
        for k, o, g, m, _ in rows:
            f.write(f"{a.sample_id}\town:{k}\t\t\t{o:.4f}\n")
            f.write(f"{a.sample_id}\tagainst:{k}\t\t\t{g:.4f}\n")

    with open(a.out_prefix + ".capture_kit.call.tsv", "w") as f:
        f.write("sample_id\tlabel\tcall\tbest_kit\town\tagainst\tmargin\tcore_depth\tfloor_to_core\n")
        f.write(f"{a.sample_id}\t{label}\t{call}\t{best[0]}\t{best[1]:.4f}\t{best[2]:.4f}"
                f"\t{best[3]:.4f}\t{core:.2f}\t{floor_core:.4f}\n")


if __name__ == "__main__":
    main()
