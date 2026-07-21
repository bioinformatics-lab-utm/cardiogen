#!/usr/bin/env python3
"""
Listeaza variantele depistate de pipeline si ABSENTE in referinta de laborator,
pastrand doar apelurile reale si explicand de ce lipsesc de la laborator.

De ce e nevoie de scriptul asta si nu ajunge 0000.vcf din bcftools isec:

  1. DeepVariant scrie in VCF si inregistrari `RefCall` (FILTER=RefCall, GT=0/0,
     QUAL=0) = "aici NU e varianta". bcftools isec compara doar CHROM:POS:REF:ALT
     si le numara ca apeluri. Ele sunt ~96% din 0000.vcf. Scriptul le exclude.
  2. VCF-urile de laborator sunt hg19, liftate la hg38 cu CrossMap. CrossMap NU
     lifteaza tot: scrie esecurile in <output>.unmap si le lasa afara tacut. O
     varianta de laborator pierduta acolo face ca apelul corect al pipeline-ului
     sa apara fals ca "gasit doar de mine". Scriptul lifteaza pozitiile din
     .unmap (hg19 -> hg38) si marcheaza aceste cazuri.

Coloana `reason` spune ce e fiecare varianta:
    real                -> pipeline o apeleaza, laboratorul nu o are deloc. Astea
                           sunt diferentele adevarate, de analizat.
    liftover_ref_swap   -> laboratorul O ARE, dar hg38 si-a schimbat referinta la
                           acea pozitie (CrossMap: Fail(REF==ALT)). Acord real.
    liftover_unmapped   -> laboratorul O ARE, dar pozitia nu se mapeaza pe hg38
                           (CrossMap: Fail(Unmap)). Artefact de conversie.

Un rand per varianta. Coloana `found_in` spune cine a depistat-o (in aceeasi proba):
    both                -> apelata de toate metodele => cea mai credibila,
                           doua metode independente nu greseau la fel,
    deepvariant / gatk  -> doar acea metoda o are.
Metricile fiecarei metode stau in coloane proprii (deepvariant_gq, gatk_gq, ...),
goale daca metoda nu a apelat varianta.

Implicit se pastreaza doar apelurile care trec pragurile de calitate (coloana
`qc`): GQ>=20, DP>=20, si VAF coerent cu genotipul (0/1 -> 0.25-0.75, 1/1 ->
>=0.90). Un 0/1 cu VAF=0.08 nu e varianta, e zgomot. Pragurile sunt reglabile
(--min-gq, --min-dp, ...), iar --keep-all pastreaza tot si scrie motivul in `qc`.
`found_in` se calculeaza DUPA filtrare, deci "both" inseamna ca ambele metode o
apeleaza SI ambele apeluri sunt de calitate.

GATK nu emite campul VAF, dar are AD, deci VAF-ul se calculeaza de acolo -
altfel filtrul pe VAF ar taia doar DeepVariant.

Rulare (zero-config):
    python3 list_pipeline_only.py

Filtrare rapida:
    awk -F'\\t' '$7=="both"'    pipeline_only_variants.tsv   # confirmate de ambele
    awk -F'\\t' '$8=="real"'    pipeline_only_variants.tsv   # fara artefacte liftover

Optional:
    python3 list_pipeline_only.py --work-dir /alt/work -o /alt/iesire
"""

import argparse
import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path

NEXTFLOW_DIR = Path(__file__).resolve().parent.parent
DEFAULT_ISEC_DIR = NEXTFLOW_DIR / "results" / "06_comparison" / "bcftools_isec"
DEFAULT_CHAIN = NEXTFLOW_DIR / "reference" / "bed_file" / "hg19ToHg38.over.chain.gz"

# De unde se citesc fisierele CrossMap .unmap (pentru adnotarea artefactelor de
# liftover). Prima locatie care exista castiga. `unmap_backup` e copia de 2 MB a
# lui work/**/*.unmap, facuta ca analiza sa nu depinda de work/ (106 GB, nemutat
# intre sisteme). Vezi docs/HANDOFF_pipeline_only_analysis.md sectiunea 12.
DEFAULT_UNMAP_DIRS = [
    NEXTFLOW_DIR / "results" / "06_comparison" / "unmap_backup",
    NEXTFLOW_DIR / "work",
]

NON_VARIANT_GT = {"0/0", "0|0", "./.", ".|."}
HET_GT = {"0/1", "1/0", "0|1", "1|0"}
HOM_GT = {"1/1", "1|1"}
REASON_BY_FAIL = {"Fail(REF==ALT)": "liftover_ref_swap",
                  "Fail(Unmap)": "liftover_unmapped"}

# Praguri implicite de calitate. GQ/DP urmeaza filtrele laboratorului (LowGQ<30,
# LowDP<20), dar cu GQ>=20 (~99% incredere in genotip) ca minim rezonabil.
# VAF trebuie sa fie coerent cu genotipul: un 0/1 cu VAF=0.08 nu e varianta.
DEFAULT_MIN_GQ = 20
DEFAULT_MIN_DP = 20
DEFAULT_HET_MIN, DEFAULT_HET_MAX = 0.25, 0.75
DEFAULT_HOM_MIN = 0.90


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("isec_dir", type=Path, nargs="?", default=DEFAULT_ISEC_DIR,
                   help=f"Directorul .../06_comparison/bcftools_isec (implicit: {DEFAULT_ISEC_DIR})")
    p.add_argument("--work-dir", type=Path, default=None,
                   help="Directorul cu fisierele CrossMap .unmap (implicit: auto - "
                        "cauta results/06_comparison/unmap_backup, apoi work/)")
    p.add_argument("--chain", type=Path, default=DEFAULT_CHAIN,
                   help="Chain hg19->hg38 pentru liftarea pozitiilor din .unmap")
    p.add_argument("--crossmap", default="CrossMap",
                   help="Executabilul CrossMap (implicit: din PATH)")
    p.add_argument("-o", "--outdir", type=Path, default=None,
                   help="Directorul de iesire (implicit: <isec_dir>/../pipeline_only)")

    q = p.add_argument_group("praguri de calitate")
    q.add_argument("--min-gq", type=float, default=DEFAULT_MIN_GQ,
                   help=f"GQ minim (implicit: {DEFAULT_MIN_GQ})")
    q.add_argument("--min-dp", type=float, default=DEFAULT_MIN_DP,
                   help=f"DP minim (implicit: {DEFAULT_MIN_DP})")
    q.add_argument("--het-min", type=float, default=DEFAULT_HET_MIN,
                   help=f"VAF minim pentru 0/1 (implicit: {DEFAULT_HET_MIN})")
    q.add_argument("--het-max", type=float, default=DEFAULT_HET_MAX,
                   help=f"VAF maxim pentru 0/1 (implicit: {DEFAULT_HET_MAX})")
    q.add_argument("--hom-min", type=float, default=DEFAULT_HOM_MIN,
                   help=f"VAF minim pentru 1/1 (implicit: {DEFAULT_HOM_MIN})")
    q.add_argument("--keep-all", action="store_true",
                   help="Nu elimina variantele care pica pragurile; le pastreaza in "
                        "fisier cu motivul in coloana `qc`. Atentie: cu --keep-all, "
                        "`found_in` se calculeaza pe TOATE apelurile, deci 'both' nu "
                        "mai garanteaza ca ambele treceau pragurile.")
    return p.parse_args(argv)


def gt_of(fmt_sample):
    """Genotipul din prima subcoloana a campului de esantion."""
    return fmt_sample.split(":", 1)[0] if fmt_sample else "."


def field(fmt_keys, fmt_sample, key):
    """Valoarea unei chei FORMAT (ex. GQ, DP, VAF) pentru esantion."""
    keys = fmt_keys.split(":")
    vals = fmt_sample.split(":")
    return vals[keys.index(key)] if key in keys and keys.index(key) < len(vals) else ""


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def vaf_of(fmt_keys, fmt_sample):
    """
    Fractia read-urilor cu alela alt. DeepVariant o da direct in campul VAF;
    GATK nu il emite, dar are AD (adancimi per alela), deci o calculam de acolo.
    Altfel filtrul pe VAF ar taia doar DeepVariant, ceea ce ar fi nedrept.
    """
    vaf = num(field(fmt_keys, fmt_sample, "VAF"))
    if vaf is not None:
        return vaf
    ad = [num(x) for x in field(fmt_keys, fmt_sample, "AD").split(",") if x != ""]
    if len(ad) < 2 or None in ad or sum(ad) == 0:
        return None
    return sum(ad[1:]) / sum(ad)


def qc_of(gt, gq, dp, vaf, args):
    """
    Verdict de calitate: "pass" sau lista motivelor de esec, separate prin ";".
    Un camp lipsa nu penalizeaza - raportam doar ce putem evalua.
    """
    fails = []
    if gq is not None and gq < args.min_gq:
        fails.append("low_gq")
    if dp is not None and dp < args.min_dp:
        fails.append("low_dp")
    if vaf is not None:
        if gt in HET_GT and not (args.het_min <= vaf <= args.het_max):
            fails.append("vaf_het")
        elif gt in HOM_GT and vaf < args.hom_min:
            fails.append("vaf_hom")
    return ";".join(fails) if fails else "pass"


def load_unmap(work_dir):
    """
    Aduna toate inregistrarile din fisierele CrossMap .unmap (coordonate hg19).
    Returneaza lista de (sample, chrom, pos, ref, alt, gt, reason).
    """
    records = []
    for unmap in work_dir.rglob("*.lifted.vcf.unmap"):
        sample = unmap.name.split(".lifted")[0]      # ex: 001TC_S1
        for line in unmap.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 10:
                continue
            reason = REASON_BY_FAIL.get(f[-1].strip())
            if reason is None:
                continue
            records.append((sample, f[0], f[1], f[3], f[4], gt_of(f[9]), reason))
    return records


def lift_unmap(records, chain, crossmap_exe):
    """
    Lifteaza pozitiile hg19 din .unmap la hg38 printr-un BED intermediar.
    Returneaza dict (sample, chrom_hg38, pos_hg38) -> reason.
    """
    if not records:
        return {}
    lifted = {}
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        bed_in, bed_out = tmp / "unmap.hg19.bed", tmp / "unmap.hg38.bed"
        with bed_in.open("w") as fh:
            for i, (sample, chrom, pos, ref, alt, gt, reason) in enumerate(records):
                # name = index, ca sa regasim inregistrarea dupa liftare
                fh.write(f"{chrom}\t{int(pos) - 1}\t{pos}\t{i}\n")

        cmd = [crossmap_exe, "bed", str(chain), str(bed_in), str(bed_out)]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"AVERTISMENT: liftarea .unmap a esuat ({e}). Variantele pierdute "
                  f"la liftover NU vor fi marcate.", file=sys.stderr)
            return {}

        for line in bed_out.read_text().splitlines():
            f = line.split("\t")
            if len(f) < 4 or line.startswith("#"):
                continue
            chrom_38, pos_38, idx = f[0], f[2], f[3]
            if not idx.isdigit():
                continue
            sample, _, _, _, _, _, reason = records[int(idx)]
            lifted[(sample, chrom_38, pos_38)] = reason
    return lifted


def sample_key(comp_dir, caller):
    """001TC_S1_L001_deepvariant_... -> ('001TC_S1_L001', '001TC_S1')."""
    marker = f"_{caller}_"
    sample_id = comp_dir.split(marker, 1)[0] if marker in comp_dir else comp_dir
    # cheia din .unmap nu are sufixul de lane: 001TC_S1_L001 -> 001TC_S1
    return sample_id, re.sub(r"_L\d+$", "", sample_id)


def collapse(rows, all_callers, keep_qc):
    """
    Un singur rand per varianta, nu unul per metoda. Metricile fiecarei metode
    ajung in coloane proprii (<caller>_gq, <caller>_dp, ...), goale daca metoda
    nu a apelat varianta. Asa o varianta gasita de ambele apare o data, cu
    found_in="both", si vezi dintr-o privire ce a zis fiecare.

    Coloana `found_in`: "both" daca varianta e apelata de toate metodele prezente,
    altfel numele metodelor care o au. Se calculeaza DUPA filtrarea QC, deci "both"
    inseamna "ambele metode o apeleaza SI ambele apeluri trec pragurile". O varianta
    unde DeepVariant trece iar GATK pica devine "deepvariant" - altfel "both" ar
    promite o confirmare dubla care nu exista.

    Returneaza (randuri, antet).
    """
    by_variant = {}
    for r in rows:
        by_variant.setdefault((r[0], r[2], r[3], r[4], r[5]), {})[r[1]] = r

    callers = sorted(all_callers)
    header = ["sample", "chrom", "pos", "ref", "alt", "gt", "found_in", "reason"]
    for c in callers:
        header += [f"{c}_gq", f"{c}_dp", f"{c}_vaf"] + ([f"{c}_qc"] if keep_qc else [])

    out = []
    for key in sorted(by_variant, key=lambda k: (k[0], k[1], k[2])):
        per_caller = by_variant[key]
        found = set(per_caller)
        found_in = ("both" if found == all_callers and len(all_callers) > 1
                    else ",".join(sorted(found)))
        # GT si reason nu depind de metoda: verificat, cele doua metode nu difera
        # niciodata pe GT cand cad de acord pe varianta.
        any_row = next(iter(per_caller.values()))
        row = list(key[:1]) + list(key[1:]) + [any_row[6], found_in, any_row[12]]
        for c in callers:
            r = per_caller.get(c)
            row += ([r[7], r[8], r[9]] + ([r[13]] if keep_qc else []) if r
                    else ["", "", ""] + ([""] if keep_qc else []))
        out.append(row)
    return out, header


def collect(isec_dir, unmap_hg38, args):
    """Parcurge 0000.vcf, pastreaza apelurile reale, le adnoteaza cu motiv + QC."""
    rows = []
    for caller_dir in sorted(d for d in isec_dir.iterdir() if d.is_dir()):
        caller = caller_dir.name
        for vcf in caller_dir.rglob("0000.vcf"):
            sample_id, unmap_id = sample_key(vcf.parent.name, caller)
            for line in vcf.read_text().splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split("\t")
                if len(f) < 10:
                    continue
                chrom, pos, ref, alt, qual, filt = f[0], f[1], f[3], f[4], f[5], f[6]
                gt = gt_of(f[9])
                # Apel real: nu RefCall/non-PASS, si genotip variant.
                if filt not in ("PASS", ".") or gt in NON_VARIANT_GT:
                    continue
                gq, dp = num(field(f[8], f[9], "GQ")), num(field(f[8], f[9], "DP"))
                vaf = vaf_of(f[8], f[9])
                reason = unmap_hg38.get((unmap_id, chrom, pos), "real")
                rows.append([sample_id, caller, chrom, int(pos), ref, alt, gt,
                             "" if gq is None else f"{gq:g}",
                             "" if dp is None else f"{dp:g}",
                             "" if vaf is None else f"{vaf:.3f}",
                             qual, filt, reason,
                             qc_of(gt, gq, dp, vaf, args)])
    rows.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
    return rows


def resolve_unmap_dir(explicit):
    """Prima locatie care contine .unmap castiga; explicit --work-dir bate auto."""
    candidates = [explicit] if explicit else DEFAULT_UNMAP_DIRS
    for d in candidates:
        if d and d.is_dir() and next(d.rglob("*.lifted.vcf.unmap"), None):
            return d
    return explicit or (DEFAULT_UNMAP_DIRS[0] if DEFAULT_UNMAP_DIRS else None)


def main(argv=None):
    args = parse_args(argv)
    if not args.isec_dir.is_dir():
        sys.exit(f"EROARE: directorul nu exista: {args.isec_dir}")

    unmap_dir = resolve_unmap_dir(args.work_dir)
    print(f"Citesc .unmap din {unmap_dir} ...", file=sys.stderr)
    unmap_records = load_unmap(unmap_dir) if unmap_dir else []
    if not unmap_records:
        print("  ATENTIE: niciun fisier .unmap gasit. Artefactele de liftover NU vor\n"
              "  fi marcate - cele ~13 'liftover_ref_swap' vor aparea fals ca 'real'.\n"
              "  Adu unmap_backup.tar.gz (vezi HANDOFF sectiunea 12) sau da --work-dir.",
              file=sys.stderr)
    print(f"  {len(unmap_records)} variante de laborator pierdute la liftover",
          file=sys.stderr)
    unmap_hg38 = lift_unmap(unmap_records, args.chain, args.crossmap)
    print(f"  {len(unmap_hg38)} pozitii recuperate in coordonate hg38", file=sys.stderr)

    rows = collect(args.isec_dir, unmap_hg38, args)
    all_callers = {r[1] for r in rows}

    print(f"\nPraguri: GQ>={args.min_gq:g}  DP>={args.min_dp:g}  "
          f"VAF het {args.het_min}-{args.het_max} / hom >={args.hom_min}",
          file=sys.stderr)
    n_before = len(rows)
    if not args.keep_all:
        rows = [r for r in rows if r[13] == "pass"]
        print(f"  {n_before - len(rows)} apeluri eliminate, {len(rows)} pastrate "
              f"(--keep-all le pastreaza pe toate)", file=sys.stderr)
    # DUPA filtrare: "both" = ambele metode o apeleaza SI ambele trec pragurile.
    rows, header = collapse(rows, all_callers, args.keep_all)

    outdir = args.outdir or (args.isec_dir.parent / "pipeline_only")
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / "pipeline_only_variants.tsv"
    with out_path.open("w", newline="") as out:
        # lineterminator: dialectul implicit 'excel' scrie \r\n, iar CR-ul ajunge in
        # ultima coloana si strica orice awk/grep pe ea.
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)

    real = [r for r in rows if r[7] == "real"]
    print(f"\nVariante depistate de pipeline si absente la laborator "
          f"(reason=real): {len(real)} variante", file=sys.stderr)
    groups = {}
    for r in real:
        groups.setdefault(r[6], []).append(r)
    for found_in, variants in sorted(groups.items(), key=lambda x: -len(x[1])):
        n_samp = len({r[0] for r in variants})
        n_pos = len({(r[1], r[2], r[3], r[4]) for r in variants})
        print(f"  {found_in:<12}: {len(variants):>4} variante | "
              f"{n_samp} esantioane | {n_pos} pozitii distincte", file=sys.stderr)

    col = header.index("found_in") + 1
    print(f"\n  -> {out_path}", file=sys.stderr)
    print(f"\n  Confirmate de ambele metode: "
          f"awk -F'\\t' '${col}==\"both\"' {out_path.name}", file=sys.stderr)


if __name__ == "__main__":
    main()
