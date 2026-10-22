# MDS Follow-up CITE-seq: Sample Identity Audit and Resolution

Documentation of a sample mislabelling and sequencing depth problem in a CITE-seq dataset, and how it was diagnosed and resolved before
any downstream processing.

**Status:** resolved. All 8 GEX and 8 protein libraries identified and above
sequencing depth spec.

---

## 1. Dataset

CITE-seq (5' GEX plus cell surface protein) on 8 MDS follow-up patient samples.

| Item | Detail |
|---|---|
| Antibody panel | TotalSeq-C cocktail V2, Cat 399910, Lot B446606 |
| 10x protocol | CG000734 Rev A, Chromium GEM-X Single Cell 5' v3 with Feature Barcode |
| Samples | MDS001, MDS010, MDS016, MDS023, MDS065, MDS167, MDS201, MDS212 |
| Libraries | 8 GEX (SI-TT index set) and 8 protein (SI-TN index set) |
| Depth spec | 200 Mreads per GEX library, 50 Mreads per protein library |

Data arrived across multiple submissions, delivered over several weeks
via separate Nextcloud share links.

---

## 2. The problem

Two independent problems were entangled, which is what made this confusing.

### 2.1 Sample names did not match library contents

The sequencing facility generated FASTQ filenames from the original
submission sheet rather than from the indexes actually present in each pool.
As a result some files carry a patient name that does not correspond to the
library inside them.

The library preparation technician flagged this but his corrections were incomplete:

- First he listed two renames, both for GEX libraries, one of which was
  attributed to the wrong sequencing run.
- Secondly, he stated the mismatch applied only to three **protein**
  libraries, which the data contradicts.

### 2.2 Sequencing depth below spec

A facility-side error meant 6 of 8 GEX libraries received fewer than the
requested 200 Mreads in the main run. 

### 2.3 Why this mattered

Merging by filename would have mixed reads from different patients into a
single library. No downstream QC step detects this. The error would have
propagated silently into cell type annotation, differential expression and
any clinical interpretation.

---

## 3. Diagnostic principle

**The filename is metadata. The index in the read header is the data.**

Every read in an Illumina FASTQ carries its own index sequence in the header:

```
@LH00583:195:23KFMVLT3:5:1190:43396:10733 1:N:0:CGTCAAGGGC+GAGTGACCTA
                                                  ^^^^^^^^^^ ^^^^^^^^^^
                                                      i7         i5
```

The submission sheet maps each index to a patient sample. So sample identity
can be recovered from the sequence data itself, independent of any filename,
any sample sheet, and any email.

This is the entire basis of the audit. Three properties were measured per
file:

1. **Observed i7 index** gives the intended sample identity.
2. **Index purity** distinguishes a real library (near 100 percent one index)
   from demultiplexing background (low count, mixed indexes).
3. **Read count** checks the depth spec.

---

## 4. Index reference

Taken from the submission sheet. This table is the ground truth for the audit.

### GEX libraries (SI-TT set)

| Sample | Index code | i7 | i5 |
|---|---|---|---|
| MDS001 | SI-TT-H9 | AGAACTTAGA | AAAGGACTCG |
| MDS010 | SI-TT-H10 | TTATCTAGGG | TAGAGCCTTT |
| MDS016 | SI-TT-H11 | ACAATCGATC | CATTCCGTCA |
| MDS023 | SI-TT-H12 | TGATGATTCA | CGACTCCTAC |
| MDS065 | SI-TT-B10 | GCCCGATGGA | CTAGACGATT |
| MDS167 | SI-TT-C1 | TGCGCGGTTT | TTTATCCTTG |
| MDS201 | SI-TT-B11 | TCTTACTTGC | CTAGAGGTCA |
| MDS212 | SI-TT-B12 | CGTCAAGGGC | GAGTGACCTA |

### Protein libraries (SI-TN set)

| Sample | Index code | i7 | i5 |
|---|---|---|---|
| MDS001 | SI-TN-C6 | TCAAGTAAAG | AACAGCACTG |
| MDS010 | SI-TN-C7 | GCTGCTCCCA | GACCCTGTAG |
| MDS016 | SI-TN-C8 | ATCTGTAGTT | CATTGGGCCT |
| MDS023 | SI-TN-C9 | ATGGCGCAAA | GCTTGCATCG |
| MDS065 | SI-TN-C10 | TAATGGGCAA | GCCTATGCGC |
| MDS167 | SI-TN-C11 | TAACTGTAGT | TATCGGTCAT |
| MDS201 | SI-TN-C12 | AGCCTTCTCT | ACGCGGACAG |
| MDS212 | SI-TN-D1 | AGCGCCTTGC | CACTTCGTAC |

---

## 5. Step by step resolution

### Step 1: Verify transfer integrity

Large transfers over Nextcloud can truncate silently. Checksums were verified
before anything else, so that later anomalies could not be blamed on download
failure.

```bash
find . -name "checksums.md5" | while read f; do
  (cd "$(dirname "$f")" && md5sum -c checksums.md5)
done
```

All files passed. Any anomaly found later was therefore present at the source.

### Step 2: Detect a non-library delivery

The first share link delivered a directory that passed its checksums but was
only 30 KB in total.

```bash
zcat MDS201-FU-GEX_S383_L005_R1_001.fastq.gz | wc -l
# 68  ->  17 read pairs
```

A GEM-X 5' v3 GEX library at roughly 20 to 50k read pairs per cell and ~5000
cells should be 100 to 250 Mreads, several GB per lane. 17 read pairs is
background.

The index explained why:

```bash
zcat MDS201-FU-GEX_S383_L005_R1_001.fastq.gz \
  | awk 'NR%4==1{split($2,a,":"); print a[4]}' | sort | uniq -c
#   17 CGTCAAGGGC+GAGTGACCTA
```

The file was named MDS201 but carried MDS212's index (SI-TT-B12). The
submission sheet says this pool contained MDS201 (SI-TT-B11). The run had
been demultiplexed against an index that was not in the pool, so the output
is index-hopped background and the real reads went to `Undetermined`.

**Lesson applied for the rest of the audit:** checksum validity proves the
transfer, not the content. Always check read count and index purity.

### Step 3: Audit every file by index

Applied to every run directory. Reads only the first 10k records per file, so
it completes in minutes even on a 281 GB delivery.

```bash
for f in *_R1_001.fastq.gz; do
  read -r cnt obs < <(zcat "$f" | head -40000 \
    | awk 'NR%4==1{split($2,a,":"); split(a[4],b,"+"); print b[1]}' \
    | sort | uniq -c | sort -rn | head -1)
  echo -e "$f\t$obs\t$cnt/10000"
done | tee index_audit.tsv
```

Output columns: filename, most frequent i7, purity out of 10000 reads.

### Step 4: Cross-check read counts against the delivered QC

The facility supplies MultiQC output with every run. Reading the counts from
there avoids decompressing whole files and allows independent verification of
the depth figures quoted by the technician.

```bash
cut -f1,5 multiqc_report_data/multiqc_fastqc.txt
```

Summed across lanes, these matched the technician's table exactly, which
confirmed his table was labelled by current (incorrect) filename.

### Step 5: Reconcile against the submission sheet

This is where the picture closed. The submission sheet's listing of which
samples were pooled in each run was cross-referenced against the libraries
found to be below spec. The match was exact, which identified a missing
delivery by deduction rather than by asking.

### Step 6: Build a corrected file structure without renaming in place

Originals are left untouched. See section 8.

---

## 6. Findings per run

| Submission | Flowcell | Lanes | Content | Verdict |
|---|---|---|---|---|
| FS63996037_NGS000665 | B23JYKCLT3 | L003, L004 | 8 GEX + 8 protein | All real. 3 GEX misnamed |
| FS63996038_NGS000713 | B23KCMTLT3 | L006 | 5 real GEX + 1 background | 1 GEX misnamed, 1 file discarded |
| FS63996039_NGS000717 | B23KFMVLT3 | L005 | 17 read pairs | Discarded, background |
| FS63996040_NGS000745 | n/a | n/a | Never delivered | Redundant, not required |
| FS63996048_NGS000747 | B23KCVVLT3 | L001 to L003 | MDS212 GEX top-up | Real, correctly named |

### Run 665

All 16 libraries genuine: 4.5 to 9 GB per lane, 100 percent index pure.

All 8 protein libraries correctly named.

Three GEX libraries form a three-way cyclic rotation. The other five are
correct.

### Run 713

Five genuine libraries, matching the submission sheet's listing for this pool
exactly (SI-TT H9, H10, H11, B10, B11).

`MDS201-FU-GEX_S199_L006` is background: 415 read pairs, 33 KB, 41 percent
CGTCAAGGGC. That index (MDS212) was not in this pool. Same failure mode as
run 717. Discarded.

`MDS167-FU-GEX_S198_L006` carries TCTTACTTGC and is MDS201, the same rotation
seen in run 665. Merging this by filename would have contaminated MDS167 with
MDS201 reads.

### Run 747

Single library, MDS212 GEX, correctly named, 434.3 Mreads across three lanes.

---

## 7. Corrected mapping

### GEX libraries requiring rename

| Run | Current filename | Observed i7 | Index code | True sample |
|---|---|---|---|---|
| 665 | `MDS167-FU-GEX_S192` | TCTTACTTGC | SI-TT-B11 | **MDS201** |
| 665 | `MDS201-FU-GEX_S193` | CGTCAAGGGC | SI-TT-B12 | **MDS212** |
| 665 | `MDS212-FU-GEX_S194` | TGCGCGGTTT | SI-TT-C1 | **MDS167** |
| 713 | `MDS167-FU-GEX_S198` | TCTTACTTGC | SI-TT-B11 | **MDS201** |

### Correct as delivered

- Run 665 GEX: MDS001, MDS010, MDS016, MDS023, MDS065
- Run 665 protein: all 8
- Run 713 GEX: MDS001, MDS010, MDS016, MDS065
- Run 747 GEX: MDS212

### Discarded

- `FS63996039_NGS000717` entire delivery (17 read pairs)
- `FS63996038_NGS000713/MDS201-FU-GEX_S199_L006` (415 read pairs, 41 percent pure)

### Final GEX depth per patient

| Sample | Source files | 665 | 713 | 747 | Total Mreads |
|---|---|---|---|---|---|
| MDS001 | `MDS001_S187`, `MDS001_S194` | 159.0 | 188.8 | | 347.8 |
| MDS010 | `MDS010_S188`, `MDS010_S195` | 129.6 | 152.7 | | 282.3 |
| MDS016 | `MDS016_S189`, `MDS016_S196` | 162.7 | 176.3 | | 339.0 |
| MDS023 | `MDS023_S190` | 256.7 | | | 256.7 |
| MDS065 | `MDS065_S191`, `MDS065_S197` | 141.6 | 166.9 | | 308.5 |
| MDS167 | `MDS212_S194` | 218.9 | | | 218.9 |
| MDS201 | `MDS167_S192`, `MDS167_S198` | 186.7 | 209.9 | | 396.6 |
| MDS212 | `MDS201_S193`, `MDS212_S246` | 147.3 | | 434.3 | 581.6 |

All libraries exceed the 200 Mreads spec. Protein libraries range 92.4 to
112.4 Mreads against a 50 Mreads spec.

**Depth imbalance note:** GEX depth spans 218.9 to 581.6 Mreads, a 2.7-fold
range, with MDS212 the outlier. Sequencing saturation will differ across
samples, which affects gene detection rates and therefore cross-patient
comparison of cell states. Handle explicitly: either downsample MDS212 for
cross-patient analyses, or record saturation as a covariate at the
normalisation and pseudobulk stage.

---

## 8. Building the corrected structure

Symlinks, not renames. The original delivery stays byte-identical and the
mapping remains auditable.

Two collisions make a flat directory unsafe:

- `S194` is MDS212 in run 665 and MDS001 in run 713
- MDS212 has an `L003` in both run 665 and run 747

Both are resolved by one directory per sample with lanes renumbered
sequentially.

```bash
R665=/path/to/FS63996037_NGS000665/B23JYKCLT3
R713=/path/to/FS63996038_NGS000713/B23KCMTLT3
R747=/path/to/FS63996048_NGS000747/B23KCVVLT3
OUT=/path/to/fastq_clean

cat > gex_map.tsv <<EOF
MDS001-FU-GEX	$R665	MDS001-FU-GEX_S187
MDS001-FU-GEX	$R713	MDS001-FU-GEX_S194
MDS010-FU-GEX	$R665	MDS010-FU-GEX_S188
MDS010-FU-GEX	$R713	MDS010-FU-GEX_S195
MDS016-FU-GEX	$R665	MDS016-FU-GEX_S189
MDS016-FU-GEX	$R713	MDS016-FU-GEX_S196
MDS023-FU-GEX	$R665	MDS023-FU-GEX_S190
MDS065-FU-GEX	$R665	MDS065-FU-GEX_S191
MDS065-FU-GEX	$R713	MDS065-FU-GEX_S197
MDS167-FU-GEX	$R665	MDS212-FU-GEX_S194
MDS201-FU-GEX	$R665	MDS167-FU-GEX_S192
MDS201-FU-GEX	$R713	MDS167-FU-GEX_S198
MDS212-FU-GEX	$R665	MDS201-FU-GEX_S193
MDS212-FU-GEX	$R747	MDS212-FU-GEX_S246
EOF

while IFS=$'\t' read -r sample dir prefix; do
  mkdir -p "$OUT/$sample"
  n=$(ls "$OUT/$sample"/*_R1_001.fastq.gz 2>/dev/null | wc -l)
  for r1 in "$dir/${prefix}"_L*_R1_001.fastq.gz; do
    n=$((n+1)); lane=$(printf "L%03d" $n)
    ln -s "$r1" "$OUT/$sample/${sample}_S1_${lane}_R1_001.fastq.gz"
    ln -s "${r1/_R1_/_R2_}" "$OUT/$sample/${sample}_S1_${lane}_R2_001.fastq.gz"
  done
done < gex_map.tsv
```

Protein libraries are correctly named in run 665 and need only a direct
symlink.

### Mandatory verification before CellRanger

```bash
for d in $OUT/*/; do
  echo -n "$(basename $d) "
  zcat $(ls $d/*_R1_001.fastq.gz | head -1) | head -4000 \
    | awk 'NR%4==1{split($2,a,":"); split(a[4],b,"+"); print b[1]}' \
    | sort -u | tr '\n' ' '
  echo
done
```

Each directory must report exactly one index, and it must be the index
assigned to that patient in section 4. This is the last point at which a
mapping error is still visible.

---


## 10. Reusable lessons

1. **Checksums validate transfer, not content.** A 30 KB FASTQ passed
   `md5sum -c` cleanly. Always check read count and index purity as well.

2. **Treat the index as ground truth.** Filenames, sample sheets and emails
   are all secondary. The index is in the data.

3. **Index purity separates a library from background.** A near-100 percent
   single index with high read count is a library. A low count with mixed
   indexes is demultiplexing output from an index that was not in the pool.

4. **A consistent rotation is a sample-sheet error, not random corruption.**
   The same rotation appeared in two independent runs, which pointed to the
   demultiplexing configuration rather than to sample handling.

5. **Audit every file, including the ones reported as correct.** The
   technician's list of renames was incomplete. One of three swaps in run 665
   was never reported, and it was the one that changed which libraries counted
   as being at spec.

6. **Reconcile inventories arithmetically.** Matching the set of below-spec
   libraries against the submission sheet's pool listings identified a missing
   delivery by deduction, which made the request to the facility specific and
   actionable rather than a general complaint.

7. **Symlink, never rename in place.** Preserves the original delivery and
   keeps the mapping reversible and auditable.

8. **Verify immediately before the irreversible step.** Once CellRanger runs,
   a mapping error is baked into the count matrix and invisible.

---

## 11. Next steps

1. Build the corrected symlink structure (section 8).
2. Run the mandatory index verification (section 8).
3. Run `cellranger multi`, pairing each GEX library with its correspondingly
   named protein library.
4. Review per-library median reads per cell and sequencing saturation.
5. Run SNP concordance against the DN dataset.
6. Decide the depth normalisation strategy for cross-patient comparison.
