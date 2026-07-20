#!/usr/bin/env python3
"""
Compara variantele "pipeline-only" (0000.vcf din bcftools isec) intre metodele
de variant calling, peste TOATE datele, dintr-o singura pornire.

La pornire scriptul:
  1. gaseste singur directorul .../06_comparison/bcftools_isec (relativ la el
     insusi, sau il poti da explicit ca argument),
  2. descopera automat toate metodele prezente (subdirectoarele: deepvariant,
     gatk, ...),
  3. compara FIECARE pereche de metode, per esantion, pe cheia CHROM:POS:REF:ALT
     (VCF-urile sunt deja normalizate de modulul BCFTOOLS_ISEC),
  4. scrie cate un TSV agregat per pereche + afiseaza un rezumat global.

0000.vcf = variante apelate de pipeline, absente in referinta de laborator
("pipeline-only" / FP-like). Intrebarea la care raspunde: aceleasi variante
sunt depistate de ambele metode, sau difera, si care anume?

Rulare (zero-config, analizeaza tot):
    python3 compare_isec_callers.py

Optional:
    python3 compare_isec_callers.py /alt/bcftools_isec --isec-file 0002.vcf
"""

import argparse
import csv
import sys
from itertools import combinations
from pathlib import Path

# Directorul implicit, relativ la locatia scriptului (nextflow/bin/ -> nextflow/results/...)
DEFAULT_ISEC_DIR = (Path(__file__).resolve().parent.parent
                    / "results" / "06_comparison" / "bcftools_isec")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("isec_dir", type=Path, nargs="?", default=DEFAULT_ISEC_DIR,
                   help=f"Directorul .../06_comparison/bcftools_isec "
                        f"(implicit: auto, {DEFAULT_ISEC_DIR})")
    p.add_argument("--isec-file", default="0000.vcf",
                   help="Care fisier isec se compara. Implicit: 0000.vcf "
                        "(pipeline-only). Foloseste 0002.vcf pentru concordante.")
    p.add_argument("-o", "--outdir", type=Path, default=None,
                   help="Directorul pentru TSV-urile de iesire. "
                        "Implicit: <isec_dir>/../compare_isec_callers")
    return p.parse_args(argv)


def load_variants(vcf_path):
    """Returneaza set de chei CHROM:POS:REF:ALT dintr-un VCF (poate lipsi)."""
    variants = set()
    if not vcf_path.is_file():
        return variants
    with vcf_path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 5:
                continue
            chrom, pos, _id, ref, alt = f[0], f[1], f[2], f[3], f[4]
            for a in alt.split(","):
                variants.add(f"{chrom}:{pos}:{ref}:{a}")
    return variants


def discover_callers(isec_dir):
    """Subdirectoarele de prim nivel = metodele de variant calling prezente."""
    return sorted(d.name for d in isec_dir.iterdir() if d.is_dir())


def collect_samples(caller_dir, isec_file):
    """
    Mapeaza sample_id -> set de variante, parcurgand
    <caller_dir>/<aligner>/<qc>/<reference>/<SAMPLE_..._caller_...>/<isec_file>.
    Cheia esantionului = numele folderului fara sufixul caller-ului incolo.
    """
    samples = {}
    caller = caller_dir.name
    marker = f"_{caller}_"
    for vcf in caller_dir.rglob(isec_file):
        comp_dir = vcf.parent.name  # ex: 001GN_S1_L001_deepvariant_bwamem_fastp_hg38
        sample_id = comp_dir.split(marker, 1)[0] if marker in comp_dir else comp_dir
        if sample_id in samples:
            print(f"AVERTISMENT: esantion duplicat '{sample_id}' in {caller_dir}",
                  file=sys.stderr)
        samples[sample_id] = load_variants(vcf)
    return samples


def compare_pair(a, set_a, b, set_b, out_path):
    """Scrie TSV-ul pereche a vs b si returneaza totalurile (a_only, b_only, shared)."""
    samples = sorted(set(set_a) | set(set_b))
    header = ["sample", f"{a}_only", f"{b}_only", "shared",
              "total_union", "jaccard", "status"]

    rows = []
    tot_a = tot_b = tot_shared = 0
    for s in samples:
        va, vb = set_a.get(s), set_b.get(s)
        if va is None or vb is None:
            rows.append([s, "", "", "", "", "",
                         f"lipseste_la_{a if va is None else b}"])
            continue
        shared, a_only, b_only = va & vb, va - vb, vb - va
        union = va | vb
        jac = f"{len(shared) / len(union):.3f}" if union else "NA"
        if not union:
            status = "zero_variante"
        elif not a_only and not b_only:
            status = "identice"
        else:
            status = "difera"
        rows.append([s, len(a_only), len(b_only), len(shared),
                     len(union), jac, status])
        tot_a += len(a_only)
        tot_b += len(b_only)
        tot_shared += len(shared)

    tot_union = tot_a + tot_b + tot_shared
    tot_jac = f"{tot_shared / tot_union:.3f}" if tot_union else "NA"
    rows.append(["TOTAL", tot_a, tot_b, tot_shared, tot_union, tot_jac, ""])

    with out_path.open("w", newline="") as out:
        # lineterminator: dialectul implicit 'excel' scrie \r\n, iar CR-ul ajunge in
        # ultima coloana si strica orice awk/grep pe ea.
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)

    return len(samples), tot_a, tot_b, tot_shared, tot_jac


def main(argv=None):
    args = parse_args(argv)
    isec_dir = args.isec_dir
    if not isec_dir.is_dir():
        sys.exit(f"EROARE: directorul nu exista: {isec_dir}")

    callers = discover_callers(isec_dir)
    if len(callers) < 2:
        sys.exit(f"EROARE: am gasit < 2 metode in {isec_dir}: {callers}")

    outdir = args.outdir or (isec_dir.parent / "compare_isec_callers")
    outdir.mkdir(parents=True, exist_ok=True)

    # Incarca o singura data setul de variante pentru fiecare metoda.
    print(f"Metode gasite: {', '.join(callers)}", file=sys.stderr)
    loaded = {c: collect_samples(isec_dir / c, args.isec_file) for c in callers}

    base = Path(args.isec_file).stem  # ex: 0000
    print(f"\nComparatie pe {args.isec_file} (pipeline-only). "
          f"shared = depistat de ambele metode.\n", file=sys.stderr)
    for a, b in combinations(callers, 2):
        out_path = outdir / f"compare_{base}_{a}_vs_{b}.tsv"
        n, ta, tb, ts, tj = compare_pair(a, loaded[a], b, loaded[b], out_path)
        print(f"  {a} vs {b}: {n} esantioane | "
              f"{a}-only={ta}  {b}-only={tb}  shared={ts}  jaccard={tj}",
              file=sys.stderr)
        print(f"    -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
