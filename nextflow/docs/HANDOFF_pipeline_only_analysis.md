# Handoff — analiza „variante depistate de pipeline, absente la laborator"

> Document pentru continuarea proiectului pe alt HPC. Rezumă o sesiune de analiză
> asupra comparației pipeline-vs-laborator. Codul și rezultatele sunt în repo;
> aici sunt **contextul, deciziile și raționamentul** care nu se văd din cod.
> Ultimul commit la scrierea acestui doc: `643a3d4` (branch `nextflow`).

---

## 1. Întrebarea de la care s-a pornit

Comparația pipeline vs. referința de laborator raporta **~22300 de variante
„pipeline-only"** (fișierele `0000.vcf` din `bcftools isec`, agregate în
`results/06_comparison/compare_isec_callers/variants_0000_deepvariant_vs_gatk.tsv`).
Utilizatorul (Dubciuc / ddubciuc@gmail.com) a întrebat, pe drept: *am cumva o
grămadă de False Positives?*

**Răspunsul, pe scurt: NU.** Cele 22300 nu erau variante suspecte — erau în
proporție de ~95% non-variante și zgomot. Numărătoarea era greșită, nu
calling-ul. Vezi secțiunea 3.

## 2. Context tehnic esențial (ușor de înțeles greșit)

- **Datele sunt PANEL țintit**, nu exom și nu genom: TruSight Cardio,
  3251 de regiuni, ~174 gene. BED: `reference/bed_hg38/TruSight_Cardio_TargetedRegions_v1.0.hg38.bed`.
- **Referința de laborator e hg19**, liftată la hg38 cu CrossMap
  (`modules/liftover.nf`). VCF-uri Illumina/Starling (filtre LowGQX, LowDP,
  LowVariantFreq, SB, R8 etc.).
- **Un singur aligner (bwamem).** DeepVariant și GATK citesc **ACELAȘI BAM**.
  Consecință critică: `found_in=both` NU e confirmare independentă — e aceeași
  aliniere judecată de două ori. Nu protejează de erori de mapare.
- **VCF-urile GATK sunt NEFILTRATE**: `FILTER="."` peste tot. Fără VQSR, fără
  hard filters. Orice comparație GATK e brută până se aplică filtrarea.
- **hap.py din pipeline (`modules/happy.nf`) compară cu LABORATORUL**, nu cu
  GIAB (toată config GIAB e comentată în `nextflow.config:39-56`). Truth =
  VCF laborator liftat, regiune = BED cardio, `--pass-only`.
- Utilizatorul a rulat **separat** un benchmark GIAB (HG002), rezultate raportate
  de el: DeepVariant+bwamem SNP Recall 0.9701 / Precision **0.9989** / F1 0.9843;
  INDEL 0.9730 / 0.9971 / 0.9849. GATK ușor mai slab (SNP prec 0.9883, INDEL 0.9824).
  **Concluzie: pipeline-ul e solid.** Diferențele față de laborator sunt despre
  laborator sau despre liftover, nu despre calling.
  (Notă: precizia GIAB e măsurată în regiuni high-confidence, care exclud zonele
  grele — deci nu se transferă automat la pozițiile care diferă de laborator.)

## 3. De ce 22300 era o cifră falsă

Descompunerea celor ~22300 (DeepVariant, `0000.vcf` agregat):

| Ce conțineau                                          | ~Cât   | Variante reale? |
|-------------------------------------------------------|--------|-----------------|
| Înregistrări `RefCall` (GT 0/0 sau ./., QUAL 0)       | 21256  | **NU** — DeepVariant spune explicit „aici NU e variantă" |
| Apeluri reale dar sub praguri de calitate             | ~1300  | Nu, zgomot      |
| Artefacte de liftover (vezi 4)                        | ~13    | Nu              |
| Candidați reali                                       | 345    | Da (cu rezerve) |

Cauza rădăcină: `bcftools isec` compară doar `CHROM:POS:REF:ALT` și numără
înregistrările `RefCall` ca și cum ar fi apeluri. DeepVariant emite multe
`RefCall`; GATK nu — de asta DeepVariant „pierdea" 21066 vs GATK 319, ceea ce
părea alarmant dar era doar un artefact de format.

## 4. Liftover-ul — verificat, NU e vinovatul principal

CrossMap scrie variantele neconvertite în `<output>.unmap`, dar `modules/liftover.nf`
**nu le publică** — se pierd tăcut în `work/`. Le-am recuperat de acolo.

- **3764 variante de laborator pierdute la liftover** (226 probe):
  3042 `Fail(REF==ALT)` + 722 `Fail(Unmap)`.
- `Fail(REF==ALT)`: hg38 și-a schimbat referința la acea poziție. 2650 din ele
  sunt GT 1/1 → în hg38 proba = referința → pipeline-ul corect nu apelează nimic.
- Doar **13** din 715 pipeline-only DeepVariant sunt artefacte reale de liftover
  (marcate `reason=liftover_ref_swap` în TSV). Restul e curat.
- Cele 722 `Fail(Unmap)` = doar 38 poziții unice, repetate; **niciuna nu cade în
  panelul cardio** → irelevante.

## 5. Metodologia de filtrare (implementată în cod)

Script: **`bin/list_pipeline_only.py`** (rulează zero-config, ~360 linii).
Output: **`results/06_comparison/pipeline_only/pipeline_only_variants.tsv`**.

Ce face, în ordine:
1. Citește `0000.vcf` din `results/06_comparison/bcftools_isec/<caller>/.../`.
2. **Exclude non-variantele**: păstrează doar `FILTER∈{PASS,.}` și GT variant
   (nu 0/0, ./., etc.). Asta elimină RefCall-urile.
3. **Adnotează artefactele de liftover**: citește `.unmap` din `work/`, lifează
   pozițiile hg19→hg38 (CrossMap din PATH), marchează coloana `reason`
   (`real` / `liftover_ref_swap` / `liftover_unmapped`).
4. **Filtrare de calitate** (praguri = filtrele laboratorului, cu GQ minim 20):
   - `GQ >= 20`, `DP >= 20`
   - VAF coerent cu genotipul: het (0/1) ∈ [0.25, 0.75], hom (1/1) >= 0.90
   - **VAF pentru GATK se calculează din AD** (GATK nu emite câmpul VAF).
     Decisiv: fără asta, filtrul pe VAF ar tăia doar DeepVariant. A redus
     „doar gatk" de la 389 la 42.
   - Reglabil: `--min-gq`, `--min-dp`, `--het-min/max`, `--hom-min`, `--keep-all`.
5. **Colapsează la un rând per variantă** (nu unul per metodă). Metricile
   fiecărei metode în coloane proprii: `deepvariant_gq/dp/vaf`, `gatk_gq/dp/vaf`.
   Verificat: cele două metode nu diferă NICIODATĂ pe GT când cad de acord pe
   variantă (286/286), deci `gt` e o singură coloană.
6. **`found_in` se calculează DUPĂ filtrare**: `both` = ambele metode o apelează
   ȘI ambele trec pragurile. Dacă una pică, devine numele celeilalte.

Coloane TSV finale (14):
`sample, chrom, pos, ref, alt, gt, found_in, reason, deepvariant_gq, deepvariant_dp, deepvariant_vaf, gatk_gq, gatk_dp, gatk_vaf`

## 6. Rezultatele curente (praguri implicite)

TSV: **345 rânduri** (variante distincte per probă). Defalcare:

| found_in    | reason | Nr  | Note |
|-------------|--------|-----|------|
| both        | real   | 274 | **Setul de lucru.** Confirmate de ambele metode. |
| both        | liftover_ref_swap | 12 | De exclus (artefact conversie) |
| gatk        | real   | 42  | Cele mai slabe — GATK nefiltrat pe INFO |
| deepvariant | real   | 17  | Fără confirmare independentă |

Setul de încredere: `awk -F'\t' '$7=="both" && $8=="real"' pipeline_only_variants.tsv` → **274**.

**Structura celor 274**: stau pe doar **22 de poziții distincte**. 4 recurente
acoperă ~240 din ele:

| Poziție hg38     | Genă    | Nr probe | DP median | În VCF-uri laborator (la ALTE probe) |
|------------------|---------|----------|-----------|--------------------------------------|
| chr19:44905910 C>G | APOE  | 149      | 375       | 39/226 → **variantă reală, raportare inconsistentă** |
| chr8:143213308 T>G | GPIHBP1 | 42     | 116       | **0/226 → SUSPECT, necesită IGV** |
| chr2:73448097 TCTC>T | ALMS1 | 35     | 89        | 51/226 → reală |
| chr2:43828388 T>C | ABCG5  | 14       | 26        | 33/226 → reală (dar sub-acoperită) |

Restul de 18 poziții: 1-6 probe fiecare. Include gene relevante clinic:
LDLR (chr19:11131368, 6 probe), MYBPC3 (chr11:47343125+47343126, 2 inserții
lipite la aceeași probă → posibil un singur eveniment sau artefact aliniere),
TNNI3 (chr19:55156661, 1 probă, DP=1186 → candidat rar puternic), PCSK9, APOB,
ABCG8, TPM1, ACTA1, LAMA4, AKAP9, LTBP2, ZHX3, LDLRAP1.

**Acoperire pe probe**: 201/226 probe au ≥1 variantă `both+real`. DAR cifra e
înșelătoare — **175 au DOAR polimorfisme comune** (top 4), 24 au comune + ceva
rar, doar **2 au exclusiv ceva rar**. **Probe cu ceva de investigat clinic = 26.**

## 7. „Evoluția" rezultatelor (de explicat clar utilizatorului)

Nu s-au schimbat variantele apelate — s-a schimbat cât de corect sunt numărate.
Pipeline-ul a apelat identic tot timpul. Fiecare pas = o întrebare mai bună pe
aceleași date, nu o corecție de calling:

```
22300 → toate înregistrările din 0000.vcf (95% RefCall/zgomot)
  345 → apeluri reale care trec pragurile de calitate
  274 → + confirmate de ambele metode, fără artefacte liftover
   22 → poziții distincte (restul = aceeași variantă repetată)
    ? → câte supraviețuiesc IGV + Sanger  ← singurele cu adevărat "sigure"
```

## 8. IMPORTANT: ce înseamnă „sigur" (a nu se supravinde)

Cele 274 sunt **candidați bine filtrați, NU variante confirmate.** Rezerve, în ordine:
1. **Ambele metode = același BAM/aligner.** `both` nu prinde erorile de aliniere
   (clasa dominantă la indel-uri și regiuni omoloage).
2. **GPIHBP1**: 42 probe la tine, 0 la laborator. Nerezolvat.
3. **Nimic verificat vizual sau adnotat.** „Comun vs rar" e inferență din
   recurență, nu fapt (recurența e și semnătura unui artefact sistematic).
4. Standardul clinic rămâne **confirmare Sanger**. Niciun pipeline nu-l înlocuiește.

## 9. Pași următori recomandați (în ordine)

1. **Adnotare gnomAD/dbSNP** (VEP sau snpEff) — minute. Rezolvă definitiv
   „comun vs rar". Confirmă ipoteza că top-4 sunt polimorfisme comune.
2. **IGV** pe GPIHBP1 (chr8:143213308) și cele 2 MYBPC3 (chr11:47343125-6) —
   prinde artefactele de aliniere.
3. **Verificare batch**: probele cu prefix `G-*` par să împartă LDLR + LAMA4 →
   posibil batch/run effect, nu biologie. De verificat maparea probă↔run.
   (Vezi și proba 001GN, run 220929_KDC72: doar 42/333 variante lab cad în BED →
   posibil ALT PANEL. De exclus din comparație dacă e cazul.)
4. **Sanger** pe candidații rari rămași.

## 10. Fixuri de pipeline încă NEAPLICATE (propuse, neconfirmate de user)

- `modules/liftover.nf`: publică `*.unmap` în `results/` (acum se pierd în `work/`).
- `modules/bcftools_isec.nf`: filtrează înainte de isec
  (`bcftools view -f PASS -i 'GT="alt"'`) ca `summary.txt` să nu mai raporteze
  cifre umflate de RefCall.
- GATK: aplică hard filters (`QD<2, FS>60, MQ<40, SOR>3, MQRankSum<-12.5,
  ReadPosRankSum<-8`) — se așteaptă ca cele 42 „doar gatk" să scadă mult.
  ATENȚIE: aceste filtre folosesc câmpuri INFO, NU se pot aplica din TSV — doar pe VCF.

## 11. Bug reparat pe parcurs

`bin/list_isec_diffs.py`, `bin/compare_isec_callers.py`, `bin/list_pipeline_only.py`
scriau TSV cu `\r\n` (dialectul `csv` implicit `excel`). CR-ul ajungea în ultima
coloană → orice `awk`/`grep` pe ea eșua (inclusiv exemplele din propria lor
documentație). Reparat cu `lineterminator="\n"` în toate trei.

## 12. Note de mediu

- CrossMap NU e în `.venv`; e la `/home/nick/miniconda3/bin/CrossMap`. Pentru a
  rula `list_pipeline_only.py` pe HPC-ul nou, asigură CrossMap în PATH (altfel
  adnotarea liftover e sărită cu avertisment, dar restul merge).
- `bedtools` folosit pentru verificări ad-hoc de acoperire/intersecție.
- Fișierele `.unmap` trăiau în `work/` (106 GB, NEmutat între sisteme). Au fost
  extrase (2.5 MB) în `results/06_comparison/unmap_backup/` și arhivate ca
  `unmap_backup.tar.gz` (289 KB). Pe HPC-ul nou:
  `tar xzf unmap_backup.tar.gz -C results/06_comparison/`. Scriptul le găsește
  automat acolo (caută unmap_backup, apoi work/). Dacă lipsesc, AVERTIZEAZĂ și
  cele ~13 `liftover_ref_swap` apar fals ca `real` — nu ignora avertismentul.
- Din `work/` (106 GB) NU trebuie mutat nimic altceva: e cache Nextflow legat de
  mașina asta (18101 symlink-uri cu căi absolute care se rup la mutare) și nu
  poate fi refolosit pentru `-resume` pe alt sistem. Pentru continuarea ANALIZEI
  ajunge `results/06_comparison/` (161 MB) + cod (prin git) + unmap_backup.
