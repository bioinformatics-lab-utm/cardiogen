# Parametri folosiți la Variant Calling

Referință: `/reference/hg38/hg38.fa`
Regiuni (panel): `TruSight_Cardio_TargetedRegions_v1.0.hg38.bed`

---

## GATK HaplotypeCaller (v4.5.0.0)

| Parametru | Valoare |
|---|---|
| Regiuni țintă (`-L`) | panel TruSight Cardio |
| Padding (`--interval-padding`) | 100 bp |
| Threads (`--native-pair-hmm-threads`) | 32 |
| Memorie Java (`-Xmx`) | 100 GB |
| Ploidie | 2 (default) |
| Prag apelare (`-stand-call-conf`) | 30 (default) |
| Argumente extra | niciunul |
| Filtrare variante | niciuna (FILTER = `.`) |

---

## DeepVariant (v1.6.1)

| Parametru | Valoare |
|---|---|
| Model (`--model_type`) | WES |
| Regiuni țintă (`--regions`) | panel TruSight Cardio (fără padding) |
| Shards (`--num_shards`) | 32 |
| Output | VCF + gVCF |
| Argumente extra | niciunul |
| Filtrare variante | model intern → FILTER = `PASS` / `RefCall` |
