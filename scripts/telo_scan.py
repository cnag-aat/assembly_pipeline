#!/usr/bin/env python3
"""
telo_scan.py — Telomeric tandem repeat scanner
================================================
Scans genome assemblies for telomeric motifs with strand separation and
tandem-specificity filtering (only counts motif copies that appear in runs
of ≥ min_tandem consecutive copies).

Softmasked (lowercase) sequences are handled transparently; matching is
case-insensitive and the original case is preserved in the hits TSV so
masking quality within telomeric arrays remains visible.

Outputs (per motif):
  <prefix>.<motif>.5p_telomere.bg     — forward-strand copy counts per window
  <prefix>.<motif>.3p_telomere.bg     — reverse-complement copy counts per window
  <prefix>.<motif>.5p+3p_telomere.bg  — combined counts per window
  <prefix>.<motif>.hits.tsv                 — every tandem run: scaffold, start, end,
                                              strand, n_copies, sequence

Track names follow the PretextGraph convention so the 5 and 3 keyboard shortcuts
work out of the box. Note: if a scaffold is inverted, what was 5p becomes 3p.

Usage:
  python telo_scan.py -i assembly.fasta -m TTAGGG [TTAGGG ...] \\
      [--window 10000] [--min-tandem 2] [--max-gap 0] [-o prefix]
      [--threads N]

Arguments:
  -i / --input       Input FASTA (can be gzipped)
  -m / --motifs      One or more telomeric motifs (5'→3', forward strand)
  --window           Window size in bp (default: 10000)
  --min-tandem       Minimum consecutive copies to count (default: 2)
  --max-gap          Allowed extra bases between copies (default: 0 = strict
                     tandem; set >0 for imperfect arrays)
  -o / --output      Output prefix (default: telo_out)
  --no-hits          Skip writing the per-run TSV (faster for large genomes)
  --threads          Worker processes (default: all available cores)
"""

import argparse
import gzip
import multiprocessing
import os
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def open_fasta(path: str):
    """Yield (name, sequence) tuples from a FASTA file (plain or gzipped).

    Original case is preserved so softmasking information is not lost.
    """
    opener = gzip.open if path.endswith(".gz") else open
    name, parts = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)  # preserve case for softmask support
    if name is not None:
        yield name, "".join(parts)


def find_tandem_runs(seq: str, motif: str, min_tandem: int, max_gap: int):
    """
    Locate all positions where `motif` appears in tandem runs of ≥ min_tandem
    copies (allowing up to max_gap extra bases between copies).

    Matching is case-insensitive so softmasked sequences are handled correctly.
    `motif` must be uppercase.

    Returns list of (start, end, n_copies) — all 0-based, half-open coords.
    `n_copies` is the number of minimal motif units in the run.
    """
    if max_gap == 0:
        pattern = re.compile(
            f"(?:{re.escape(motif)}){{{min_tandem},}}",
            re.IGNORECASE,
        )
    else:
        spacer = f".{{0,{max_gap}}}"
        pattern = re.compile(
            f"(?:(?:{re.escape(motif)}){spacer}){{{min_tandem - 1},}}{re.escape(motif)}",
            re.IGNORECASE,
        )

    runs = []
    pos = 0
    while pos < len(seq):
        m_obj = pattern.search(seq, pos)
        if m_obj is None:
            break
        start, end = m_obj.start(), m_obj.end()
        run_seq = seq[start:end]
        # Count on uppercased slice; motif is already uppercase
        n_copies = run_seq.upper().count(motif)
        if n_copies < min_tandem:
            pos = start + 1
            continue
        runs.append((start, end, n_copies))
        pos = end  # no overlapping runs
    return runs


def build_windows(seq_len: int, window: int):
    """Yield (win_start, win_end) pairs covering the sequence."""
    for start in range(0, seq_len, window):
        yield start, min(start + window, seq_len)


# ---------------------------------------------------------------------------
# Per-window counting
# ---------------------------------------------------------------------------

def count_in_windows(runs, seq_len, window):
    """
    Given a list of (start, end, n_copies) runs, sum n_copies per window.
    Runs that span a window boundary are split proportionally by bp overlap.
    """
    wins = list(build_windows(seq_len, window))
    counts = [0.0] * len(wins)

    for (rstart, rend, ncopies) in runs:
        run_len = rend - rstart
        if run_len == 0:
            continue
        for wi, (ws, we) in enumerate(wins):
            overlap = min(rend, we) - max(rstart, ws)
            if overlap <= 0:
                continue
            counts[wi] += ncopies * overlap / run_len

    return wins, [round(c) for c in counts]


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

def write_bedgraph(path, scaffold_windows):
    """
    scaffold_windows: list of (scaffold, wins_list, counts_list)
    """
    with open(path, "w") as fh:
        for scaffold, wins, counts in scaffold_windows:
            for (ws, we), cnt in zip(wins, counts):
                fh.write(f"{scaffold}\t{ws}\t{we}\t{cnt}\n")


def write_hits_tsv(path, all_hits):
    """
    all_hits: list of (scaffold, start, end, strand, n_copies, sequence)
    """
    with open(path, "w") as fh:
        fh.write("scaffold\tstart\tend\tstrand\tn_copies\tsequence\n")
        for row in all_hits:
            fh.write("\t".join(str(x) for x in row) + "\n")


# ---------------------------------------------------------------------------
# Multiprocessing worker
# ---------------------------------------------------------------------------

def _scan_scaffold(args_tuple):
    """Worker: scan one scaffold for all motifs. Returns structured results."""
    scaffold, seq, canonical, min_tandem, max_gap, window, no_hits = args_tuple
    seq_len = len(seq)
    result = {}
    for fwd_motif, rev_motif in canonical.items():
        fwd_runs = find_tandem_runs(seq, fwd_motif, min_tandem, max_gap)
        rev_runs = find_tandem_runs(seq, rev_motif, min_tandem, max_gap)

        wins_f, cnts_f = count_in_windows(fwd_runs, seq_len, window)
        wins_r, cnts_r = count_in_windows(rev_runs, seq_len, window)
        cnts_c = [a + b for a, b in zip(cnts_f, cnts_r)]

        hits = []
        if not no_hits:
            for start, end, nc in fwd_runs:
                hits.append((scaffold, start, end, "+", nc, seq[start:end]))
            for start, end, nc in rev_runs:
                hits.append((scaffold, start, end, "-", nc, seq[start:end]))

        result[fwd_motif] = {
            "fwd": (scaffold, wins_f, cnts_f),
            "rev": (scaffold, wins_r, cnts_r),
            "comb": (scaffold, wins_f, cnts_c),
            "hits": hits,
        }
    return scaffold, seq_len, result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Telomeric tandem repeat scanner with strand separation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("-i", "--input", required=True, help="Input FASTA")
    ap.add_argument(
        "-m", "--motifs", nargs="+", required=True,
        help="Telomeric motif(s), 5'→3' on the forward strand (e.g. TTAGGG)"
    )
    ap.add_argument("--window", type=int, default=10000, help="Window size bp (default 10000)")
    ap.add_argument(
        "--min-tandem", type=int, default=2,
        help="Min consecutive copies to qualify as telomeric (default 2)"
    )
    ap.add_argument(
        "--max-gap", type=int, default=0,
        help="Max extra bases allowed between copies (default 0 = strict tandem)"
    )
    ap.add_argument("-o", "--output", default="telo_out", help="Output prefix")
    ap.add_argument("--no-hits", action="store_true", help="Skip per-run TSV output")
    ap.add_argument(
        "--threads", type=int, default=os.cpu_count(),
        help="Worker processes for parallel scaffold scanning (default: all cores)"
    )
    args = ap.parse_args()

    motifs_input = [m.upper() for m in args.motifs]

    # Deduplicate: if user supplies both TTAGGG and its revcomp CCCTAA,
    # treat them as the same canonical motif pair.
    canonical = {}  # fwd_motif -> rev_motif
    seen = set()
    for m in motifs_input:
        rc = revcomp(m)
        if m in seen or rc in seen:
            print(f"[info] Skipping {m} — already covered as revcomp of a prior motif", file=sys.stderr)
            continue
        seen.add(m)
        seen.add(rc)
        canonical[m] = rc
        print(f"[info] Motif: {m}  revcomp: {rc}", file=sys.stderr)

    # Storage: motif -> strand -> list of (scaffold, wins, counts)  [ordered]
    fwd_data = {m: [] for m in canonical}
    rev_data = {m: [] for m in canonical}
    comb_data = {m: [] for m in canonical}
    hit_data = {m: [] for m in canonical}

    print(
        f"[info] Scanning {args.input}  window={args.window}  "
        f"min_tandem={args.min_tandem}  max_gap={args.max_gap}  "
        f"threads={args.threads}",
        file=sys.stderr,
    )

    # Build task iterator — each item is one scaffold's worth of work
    def task_iter():
        for scaffold, seq in open_fasta(args.input):
            yield scaffold, seq, canonical, args.min_tandem, args.max_gap, args.window, args.no_hits

    scaffold_lengths = {}  # used later by the strand-label heuristic

    with multiprocessing.Pool(processes=args.threads) as pool:
        for scaffold, seq_len, result in pool.imap(_scan_scaffold, task_iter()):
            print(f"  scaffold {scaffold}  len={seq_len:,}", file=sys.stderr)
            scaffold_lengths[scaffold] = seq_len
            for fwd_motif, data in result.items():
                fwd_data[fwd_motif].append(data["fwd"])
                rev_data[fwd_motif].append(data["rev"])
                comb_data[fwd_motif].append(data["comb"])
                hit_data[fwd_motif].extend(data["hits"])

    # ---------------------------------------------------------------------------
    # Heuristic: decide which strand is 5p and which is 3p per motif.
    #
    # For each scaffold with ≥ 10 windows, sum fwd counts in the first 5 windows
    # and the last 5 windows.  If the fwd motif is predominantly at starts across
    # all qualifying scaffolds it is labelled 5p_telomere; if predominantly at
    # ends it is 3p_telomere.  Default (tie / no data) is fwd = 3p_telomere,
    # because the canonical telomeric motif (e.g. TTAGGG) is the sequence
    # synthesised by telomerase and is added to the 3' end of chromosomes.
    # ---------------------------------------------------------------------------
    TERM_WINS = 5   # number of terminal windows to inspect at each end
    MIN_WINS  = 10  # minimum windows a scaffold must have to be included

    motif_labels = {}  # fwd_motif -> (fwd_label, rev_label)
    for fwd_motif in canonical:
        start_score = end_score = 0
        for scaffold, wins, counts in fwd_data[fwd_motif]:
            n = len(wins)
            if n < MIN_WINS:
                continue
            start_score += sum(counts[:TERM_WINS])
            end_score   += sum(counts[-TERM_WINS:])

        if start_score > end_score:
            fwd_label, rev_label = "5p", "3p"
            reason = f"fwd skews to starts (start={start_score} end={end_score})"
        else:
            fwd_label, rev_label = "3p", "5p"
            if end_score > start_score:
                reason = f"fwd skews to ends (start={start_score} end={end_score})"
            else:
                reason = f"no signal or tie (start={start_score} end={end_score}) — using default"

        motif_labels[fwd_motif] = (fwd_label, rev_label)
        print(
            f"[info] {fwd_motif}: fwd→{fwd_label}_telomere  "
            f"rev→{rev_label}_telomere  ({reason})",
            file=sys.stderr,
        )

    # Write outputs
    prefix = args.output
    for fwd_motif in canonical:
        tag = fwd_motif
        fwd_label, rev_label = motif_labels[fwd_motif]

        fwd_path  = f"{prefix}.{tag}.{fwd_label}_telomere.bg"
        rev_path  = f"{prefix}.{tag}.{rev_label}_telomere.bg"
        comb_path = f"{prefix}.{tag}.5p+3p_telomere.bg"
        hits_path = f"{prefix}.{tag}.hits.tsv"

        write_bedgraph(fwd_path, fwd_data[fwd_motif])
        write_bedgraph(rev_path, rev_data[fwd_motif])
        write_bedgraph(comb_path, comb_data[fwd_motif])
        print(f"[out] {fwd_path}", file=sys.stderr)
        print(f"[out] {rev_path}", file=sys.stderr)
        print(f"[out] {comb_path}", file=sys.stderr)

        if not args.no_hits:
            hit_data[fwd_motif].sort(key=lambda x: (x[0], x[1]))
            write_hits_tsv(hits_path, hit_data[fwd_motif])
            print(f"[out] {hits_path}", file=sys.stderr)

    print("[done]", file=sys.stderr)


if __name__ == "__main__":
    main()
