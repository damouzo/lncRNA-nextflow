#!/usr/bin/env python3
"""F1 threshold derivation (cheap process, Phase A).

Reads the metric tables (no BAMs, no counts) and derives the flag thresholds from
the controls with a fixed rule. Overrides take priority; if a control set is too
small or the negatives are not separated, it falls back to the config defaults and
flags `thresholds_unreliable`.

  flank_ratio: t = percentile 10 of pos_flank.flank_ratio_median
               (keeps >= 90% of positives); report neg_flank rejection at t.
               With --no-flank (F1.3 disabled) the flank row is omitted entirely.
  antisense:   r = percentile 95 of null_antisense.ratio_median
               (a candidate is suspect when ratio_median <= r). The pos_antisense
               flag rate at r is reported as `separation`; if positives are not
               above the null the status is `thresholds_inverted`, and if too many
               positives are flagged it is `thresholds_not_specific`.

Writes:
  thresholds_used.tsv         metric, value, origin, n_controls, separation, status
  sensitivity_thresholds.tsv  alternative percentiles for F1.6
  f1_controls/<set>.tsv       the metric distribution per control set
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f1_common as C  # noqa: E402


def read_metric(path, value_col):
    """-> (header, rows) keeping elements of any set."""
    header, rows = C.read_tsv(path)
    if value_col not in header:
        raise SystemExit("ERROR: column %s not in %s" % (value_col, path))
    return header, rows


def group_by_set(rows):
    out = {}
    for r in rows:
        out.setdefault(r["set"], []).append(r)
    return out


def element_medians(rows, value_col):
    """Collapse per-norm_group rows to one value per element (median over groups)."""
    by_elem = {}
    for r in rows:
        v = _num(r.get(value_col))
        if v is not None:
            by_elem.setdefault(r["element_id"], []).append(v)
    return {e: C.median(v) for e, v in by_elem.items()}


def rejection_rate(neg_values, thr):
    """Fraction of negatives below the threshold (correctly rejected)."""
    vals = [v for v in neg_values if v is not None]
    if not vals:
        return None
    return sum(1 for v in vals if v < thr) / len(vals)


def flag_rate(values, thr):
    """Fraction of elements flagged as suspect (value <= thr). None if empty."""
    vals = [v for v in values if v is not None]
    if not vals or thr is None:
        return None
    return sum(1 for v in vals if v <= thr) / len(vals)


def decide(derived, override, default, n_controls, min_controls, separation,
           separation_min):
    """Return (value, origin, status). origin in derived/override/default."""
    if override is not None:
        return override, "override", ""
    if n_controls < min_controls or (separation is not None and separation < separation_min):
        return default, "default", "thresholds_unreliable"
    if derived is None:
        return default, "default", "thresholds_unreliable"
    return derived, "derived", ""


def _num(x):
    if x is None or x == "" or x == "NA":
        return None
    try:
        return float(x)
    except ValueError:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flank-metrics", required=True)
    ap.add_argument("--antisense-metrics", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--no-flank", action="store_true",
                    help="F1.3 disabled: skip the flank threshold entirely (exploratory code).")
    ap.add_argument("--flank-override", type=float, default=None)
    ap.add_argument("--antisense-override", type=float, default=None)
    ap.add_argument("--flank-default", type=float, default=2.0)
    ap.add_argument("--antisense-default", type=float, default=0.15)
    ap.add_argument("--flank-percentile", type=float, default=10.0)
    ap.add_argument("--antisense-percentile", type=float, default=95.0)
    ap.add_argument("--min-controls", type=int, default=200)
    ap.add_argument("--neg-reject-min", type=float, default=0.90)
    ap.add_argument("--pos-flag-max", type=float, default=0.30,
                    help="Max pos_antisense flag rate at the derived antisense "
                         "threshold before it is reported as not specific.")
    args = ap.parse_args()

    C.ensure_file(args.flank_metrics, "flank metrics")
    C.ensure_file(args.antisense_metrics, "antisense metrics")
    os.makedirs(os.path.join(args.outdir, "f1_controls"), exist_ok=True)

    _, anti_rows = read_metric(args.antisense_metrics, "ratio_median")
    aset = group_by_set(anti_rows)
    null_anti = element_medians(aset.get("null_antisense", []), "ratio_median")
    pos_anti = element_medians(aset.get("pos_antisense", []), "ratio_median")

    r_anti_derived = C.percentile(list(null_anti.values()), args.antisense_percentile)
    # antisense has no rejection separation; rely on control-set size only
    r_anti, r_origin, r_status = decide(
        r_anti_derived, args.antisense_override, args.antisense_default,
        len(null_anti), args.min_controls, None, 0.0)
    # Diagnostic: how many annotated antisense pairs the derived threshold would
    # flag, and whether positives sit above the null at all (a sign that the
    # strand convention is inverted). Does not change the value.
    pos_flag = flag_rate(list(pos_anti.values()), r_anti)
    pos_med = C.median(list(pos_anti.values()))
    null_med = C.median(list(null_anti.values()))
    if pos_anti and pos_med is not None and null_med is not None and pos_med <= null_med:
        r_status = "thresholds_inverted"
    elif pos_flag is not None and pos_flag > args.pos_flag_max and not r_status:
        r_status = "thresholds_not_specific"

    if args.no_flank:
        t_flank, t_origin, t_status, neg_rej, n_pos = None, "disabled", "", None, 0
    else:
        _, flank_rows = read_metric(args.flank_metrics, "flank_ratio_median")
        fset = group_by_set(flank_rows)
        pos_flank = element_medians(fset.get("pos_flank", []), "flank_ratio_median")
        neg_flank = element_medians(fset.get("neg_flank", []), "flank_ratio_median")
        t_flank_derived = C.percentile(list(pos_flank.values()), args.flank_percentile)
        neg_rej = rejection_rate(list(neg_flank.values()),
                                 t_flank_derived if t_flank_derived is not None else 0.0)
        t_flank, t_origin, t_status = decide(
            t_flank_derived, args.flank_override, args.flank_default,
            len(pos_flank), args.min_controls, neg_rej, args.neg_reject_min)
        n_pos = len(pos_flank)

    status = ""
    if t_status:
        status = t_status
    if r_status in ("thresholds_unreliable",):
        status = "thresholds_unreliable"
    elif r_status and not status:
        status = r_status

    rows = []
    if not args.no_flank:
        rows.append(("flank_ratio", _fmt(t_flank), t_origin,
                     n_pos, _fmt(neg_rej), t_status))
    rows.append(("antisense_ratio", _fmt(r_anti), r_origin,
                 len(null_anti),
                 _fmt(pos_flag) if pos_anti else "NA", r_status))
    C.write_tsv(rows, os.path.join(args.outdir, "thresholds_used.tsv"),
                ["metric", "value", "origin", "n_controls", "separation", "status"])

    # sensitivity candidates (F1.6)
    sens = []
    for q in (90.0, 95.0, 99.0):
        sens.append(("antisense_ratio", "p%g" % q, _fmt(C.percentile(list(null_anti.values()), q))))
    if not args.no_flank:
        _, flank_rows = read_metric(args.flank_metrics, "flank_ratio_median")
        fset = group_by_set(flank_rows)
        pos_flank = element_medians(fset.get("pos_flank", []), "flank_ratio_median")
        for q in (5.0, 10.0, 20.0):
            sens.append(("flank_ratio", "p%g" % q, _fmt(C.percentile(list(pos_flank.values()), q))))
    C.write_tsv(sens, os.path.join(args.outdir, "sensitivity_thresholds.tsv"),
                ["metric", "percentile", "value"])

    # control distributions per set
    for set_name, rows in aset.items():
        C.write_tsv([(r["element_id"], r.get("norm_group"), r.get("ratio_median"),
                      r.get("sense_reads"), r.get("antisense_reads"), r.get("rho"), r.get("n"))
                     for r in rows],
                    os.path.join(args.outdir, "f1_controls", "antisense_%s.tsv" % set_name),
                    ["element_id", "norm_group", "ratio_median", "sense_reads",
                     "antisense_reads", "rho", "n"])
    if not args.no_flank:
        _, flank_rows = read_metric(args.flank_metrics, "flank_ratio_median")
        for set_name, rows in group_by_set(flank_rows).items():
            C.write_tsv([(r["element_id"], r.get("norm_group"), r.get("flank_ratio_median"),
                          r.get("opp_strand_ratio_median"), r.get("n_samples_used"))
                         for r in rows],
                        os.path.join(args.outdir, "f1_controls", "flank_%s.tsv" % set_name),
                        ["element_id", "norm_group", "flank_ratio_median",
                         "opp_strand_ratio_median", "n_samples_used"])

    with open(os.path.join(args.outdir, "thresholds_used.txt"), "w") as fh:
        if not args.no_flank:
            fh.write("t_flank=%s (origin=%s, n_pos=%d, neg_rejection=%s)\n"
                     % (_fmt(t_flank), t_origin, n_pos, _fmt(neg_rej)))
        else:
            fh.write("t_flank=disabled (run_flank_check=false; F1.3 exploratory, not applied)\n")
        fh.write("r_antisense=%s (origin=%s, n_null=%d, n_pos=%d, pos_median=%s, "
                 "null_median=%s, pos_flag_rate=%s)\n"
                 % (_fmt(r_anti), r_origin, len(null_anti), len(pos_anti),
                    _fmt(pos_med), _fmt(null_med), _fmt(pos_flag)))
        fh.write("status=%s\n" % (status or "ok"))
    if status:
        print("WARNING: %s" % status)
    print("F1 thresholds: t_flank=%s (%s), r_antisense=%s (%s)"
          % (_fmt(t_flank), t_origin, _fmt(r_anti), r_origin))


def _fmt(x):
    return "NA" if x is None else ("%.6g" % x)


if __name__ == "__main__":
    main()
