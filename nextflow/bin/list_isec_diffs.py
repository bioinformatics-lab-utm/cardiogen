#!/usr/bin/env python3
"""
Listeaza EXACT ce variante a depistat fiecare metoda de variant calling, per
esantion: si cele care difera, si cele COMUNE (gasite de ambele).

Unde compare_isec_callers.py da doar numaratorile, acest script da variantele in
sine: o linie per varianta, cu CHROM/POS/REF/ALT si o coloana "found_in" care
spune unde a fost depistata:
    - numele metodei (ex. "deepvariant" / "gatk") -> doar acea metoda o are,
    - "both"                                       -> comuna, gasita de ambele.

Se compara fisierul 0000.vcf (pipeline-only) al fiecarei metode, peste TOATE
datele, dintr-o singura pornire. La pornire scriptul:
  1. gaseste singur .../06_comparison/bcftools_isec,
  2. descopera automat metodele prezente (deepvariant, gatk, ...),
  3. pentru fiecare pereche scrie un TSV "lung" cu toate variantele.

Coloane iesire:
    sample  chrom  pos  ref  alt  found_in
Filtrare rapida:
    grep -P '\\tboth$'  variants_...tsv   # doar comunele
    grep -vP '\\tboth$' variants_...tsv   # doar diferentele

Rulare (zero-config, peste tot):
    python3 list_isec_diffs.py

Optional:
    python3 list_isec_diffs.py /alt/bcftools_isec --isec-file 0002.vcf
"""

import argparse
import csv
import sys
from itertools import combinations
from pathlib import Path

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
    """sample VCF -> dict cheie 'CHROM:POS:REF:ALT' -> (chrom, pos, ref, alt)."""
    variants = {}
    if not vcf_path.is_file():
        return variants
    with vcf_path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 5:
                continue
            chrom, pos, ref, alt = f[0], f[1], f[3], f[4]
            for a in alt.split(","):
                variants[f"{chrom}:{pos}:{ref}:{a}"] = (chrom, pos, ref, a)
    return variants


def discover_callers(isec_dir):
    return sorted(d.name for d in isec_dir.iterdir() if d.is_dir())


def collect_samples(caller_dir, isec_file):
    """sample_id -> dict de variante (vezi load_variants)."""
    samples = {}
    caller = caller_dir.name
    marker = f"_{caller}_"
    for vcf in caller_dir.rglob(isec_file):
        comp_dir = vcf.parent.name
        sample_id = comp_dir.split(marker, 1)[0] if marker in comp_dir else comp_dir
        samples[sample_id] = load_variants(vcf)
    return samples


def pos_key(rec):
    """Sortare naturala: (sample, chrom, pos_int, ref, alt)."""
    sample, chrom, pos, ref, alt, _found = rec
    try:
        pos_i = int(pos)
    except ValueError:
        pos_i = 0
    return (sample, chrom, pos_i, ref, alt)


def diff_pair(a, samples_a, b, samples_b, out_path):
    """
    Scrie TSV-ul lung cu toate variantele a vs b (diferente + comune).
    Returneaza (n_total, n_shared, n_samples).
    """
    all_samples = sorted(set(samples_a) | set(samples_b))
    rows = []
    n_shared = 0
    for s in all_samples:
        va = samples_a.get(s, {})
        vb = samples_b.get(s, {})
        for key in set(va) - set(vb):           # doar la metoda a
            chrom, pos, ref, alt = va[key]
            rows.append([s, chrom, pos, ref, alt, a])
        for key in set(vb) - set(va):           # doar la metoda b
            chrom, pos, ref, alt = vb[key]
            rows.append([s, chrom, pos, ref, alt, b])
        for key in set(va) & set(vb):           # comuna ambelor metode
            chrom, pos, ref, alt = va[key]
            rows.append([s, chrom, pos, ref, alt, "both"])
            n_shared += 1

    rows.sort(key=pos_key)

    with out_path.open("w", newline="") as out:
        # lineterminator: dialectul implicit 'excel' scrie \r\n, iar CR-ul ajunge in
        # ultima coloana si strica orice awk/grep pe ea (inclusiv exemplele de mai sus).
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "chrom", "pos", "ref", "alt", "found_in"])
        w.writerows(rows)

    n_samples = len({r[0] for r in rows})
    return len(rows), n_shared, n_samples


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

    print(f"Metode gasite: {', '.join(callers)}", file=sys.stderr)
    loaded = {c: collect_samples(isec_dir / c, args.isec_file) for c in callers}

    base = Path(args.isec_file).stem
    print(f"\nVariante pe {args.isec_file}, per esantion "
          f"(found_in: metoda = doar ea, 'both' = comuna):\n", file=sys.stderr)
    for a, b in combinations(callers, 2):
        out_path = outdir / f"variants_{base}_{a}_vs_{b}.tsv"
        n_rows, n_shared, n_samp = diff_pair(a, loaded[a], b, loaded[b], out_path)
        n_diff = n_rows - n_shared
        print(f"  {a} vs {b}: {n_rows} variante in {n_samp} esantioane "
              f"({n_diff} difera, {n_shared} comune)", file=sys.stderr)
        print(f"    -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
