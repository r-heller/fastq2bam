# Audit and refinement — 30 October 2026

Purpose: a reproducible teaching path from a single-end FASTQ through QC,
alignment and indexed BAM inspection. The acceptance target is a fresh matched
software demonstration with known expected alignments, safe output handling,
useful error messages and synchronized English/German instructions.

## Findings and completed repairs

| Severity | Evidence and user impact | Completed action |
|---|---|---|
| High | The old demo described ShortRead ERR127302 reads as E. coli and selected an E. coli reference. ENA identifies ERR127302 as Homo sapiens, RNA-Seq, PAIRED. That combination cannot support the claimed biological example. | Replaced the default with explicitly synthetic matching reference/reads and exact position truth; preserved historical `data/` untouched. Real inputs must specify both FASTQ and reference. |
| High | Raw external command strings did not reliably preserve argument boundaries or stop before downstream use of failed/empty output. | Checked `system2` calls with quoted arguments, required outputs, complete HISAT2 indices, samtools quickcheck, and retained intermediates on failure. |
| Medium | A non-existent output path returns NA from `Sys.readlink`; default `nzchar(NA)` incorrectly caused every fresh run to be refused. | Explicit NA-aware symlink test; fresh, existing and dangling-link regressions. |
| Medium | Actual Rsubread execution rejected the original synthetic FASTA's long line, wrote only a log and returned without an R exception. The pipeline then announced a completed index. | Wrap generated FASTA to 80 bases, require `.files`, `.reads`, `.00.b.array`, `.00.b.tab` before alignment; a log alone cannot pass. |
| Medium | Implicit RNA/DNA selection and absent first-use prerequisites obscured what ran. | Explicit DNA/RNA and aligner options, Phred+33 setting, one thread, dependency instructions, local input fingerprints and matching guides. |

Source verification: [ENA read-run metadata](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERR127302&result=read_run&fields=run_accession,scientific_name,library_strategy,library_layout&format=tsv), retrieved 30 October 2026 and retained as `ENA-ERR127302.tsv` in the local audit evidence.
The [ShortRead overview](https://bioconductor.org/packages/devel/bioc/vignettes/ShortRead/inst/doc/Overview.html) identifies the example accession. These are provenance checks, not new analyses of the original reads.

## Observed validation

macOS 27.0.1 ARM64; R 4.6.1; Bioconductor 3.23; Rsubread 2.26.0,
Rsamtools 2.28.0, ShortRead 1.70.0 and GenomicAlignments 1.48.0.
Dependencies were installed in an isolated maintenance library.

- Command regressions passed: inert metacharacters/spaces, nonzero exit handling,
  empty/missing output, complete/partial HISAT2 and Rsubread indices, and fresh/
  existing/dangling output destinations.
- Synthetic fixture regressions passed: 100 unique read identities, known exact
  sequence positions, 75-base qualities, wrapped FASTA and no overwrite.
- Complete default pipeline passed with Rsubread and Rsamtools: 100/100 uniquely
  mapped reads, sorted BAM and BAI, QC HTML/PDF, coverage and provenance. All
  read identities, reference names, starts and widths matched expected truth.
- A final direct reread also verified positive strands for all 100 reads after
  adding that requirement to the pipeline. The final script parsed successfully.
- Five final QC/alignment PDFs were rendered and visually inspected: per-cycle
  quality, GC, coverage, MAPQ and aligned read length. Constant synthetic quality,
  read length and MAPQ appropriately produce degenerate distributions.
- Rsubread emits a Phred-offset heuristic warning for constant ASCII-63 qualities.
  The fixture explicitly encodes Phred+33 Q30; ShortRead measures Q30, and align
  receives `phredOffset=33`. The warning remains visible; it is not suppressed.

Evidence is under
`/Users/rh/Documents/GitHub/.repository-audit/2026-10-30-rheller-refinement/fastq2bam/`.
The successful run is ignored `results/run-13e456b01fe08/`; its provenance receipt
states the exact 100-read validation. Earlier failed attempts and logs are retained,
not presented as successful results. The original `data/` files were not changed.

## Roadmap and readiness recorded on 30 October

Essential repairs and the local single-end synthetic first-use path are complete.
The configured Linux CI has not run remotely. HISAT2, samtools and FastQC branches
were not exercised end to end on this machine; their command boundaries have
regression coverage. Windows/Linux execution and real biological inputs remain
unqualified. No publication or deployment occurred.

Before using real data, review organism/assembly, quality encoding, read type and
alignment settings and validate suitable controls. Paired-end support, streaming
large FASTQ files, statistical downstream analysis and a broader aligner comparison
are future milestones, outside this teaching release. README.md gives the exact
current prerequisites and runnable commands. Next distribution action: review the
candidate and run its Linux workflow on an authorized branch.

## Native integration update — 14 November 2026

Actual FastQC 0.13 rejected the old `-o` option. The pipeline now uses
`--outdir`. Actual HISAT2 also failed when its internally reconstructed shell
command encountered an apostrophe in a run path. Alignment now uses fixed
relative names through temporary symbolic links, without copying the large
inputs again. The working directory and links are cleaned up on success and
failure; original inputs and partial alignment outputs remain intact.

Native DNA and RNA-mode synthetic runs on macOS ARM64 completed QC, alignment,
sorting, indexing and the exact 100-read truth comparison with FastQC 0.13 and
samtools. HISAT2 was built from upstream commit
`0d244324f98de541bce04d45c75e83bc3522f7f4` using its ARM support.
The run paths contained spaces and an apostrophe. An independent 40-assertion
filesystem probe checked small/large indexes, plain/gzip read aliases, retained
input hashes, cleanup and injected alignment failures. The native CI step
requires all tools and runs both modes; its hosted result is a separate gate.

Visual inspection also found a silent ShortRead heatmap error: equal bin counts
produce nonunique colour breaks. For that case only, the report now plots the
complete original count table as points. Every cycle, score and count is
preserved; the synthetic quality values are unchanged. Unit tests cover equal
counts with multiple scores and samples, ordinary counts and an empty table.
The corrected panel was rendered and inspected, and CI rejects the observed
error text in generated PDFs.

FastQC's actual synthetic-run summary retains one FAIL (per-base sequence
content) and two WARN modules (per-sequence GC and overrepresented sequences).
Other modules pass. Successful execution and exact alignment are separate from
those diagnostic quality judgements; no all-modules-pass claim is made.

The earlier failures are retained in the audit evidence. These contiguous reads
do not qualify splice-junction detection, biological accuracy, paired-end or
large-input streaming. The historical scope and limitations above remain useful;
the earlier statement that native tools were untested is superseded by these
observed runs.
