# fastq2bam

An R teaching pipeline: inspect single-end FASTQ reads, run HISAT2 or Rsubread,
sort and index BAM, and inspect read quality and coverage.

## Install and first use

Use R >= 4.1 and the matching Bioconductor release. The verified local environment
was R 4.6.1 with Bioconductor 3.23. Install dependencies explicitly:

```r
install.packages(c("BiocManager", "data.table", "ggplot2", "R.utils"))
BiocManager::install(c("ShortRead", "Rsamtools", "GenomicAlignments",
                      "Biostrings", "BiocParallel", "Rsubread"),
                    ask = FALSE, update = FALSE)
```

From the repository root:

```sh
Rscript --vanilla tests/test_commands.R
Rscript --vanilla tests/test_demo.R
Rscript --vanilla tests/test_quality_report.R
Rscript --vanilla 00_fastq_to_bam_pipeline.R
```

The default demonstration generates a matching 20 kb synthetic reference and
100 exact 75-base single-end reads (seed 42; Phred+33, quality 30). It downloads
no sequence data. It checks every aligned read name, reference, start position
and alignment width against saved truth. The fresh `results/run-*/` directory
contains inputs, input hashes/provenance, expected positions, QC HTML/PDF,
alignment intermediates and the final BAM/BAI. This verifies software behavior;
the fixture is not a biological sample or an alignment-accuracy study.

HISAT2 plus `hisat2-build` can replace Rsubread; samtools and FastQC are optional.
`FASTQ2BAM_ALIGNER=auto` prefers HISAT2, then Rsubread. Select `hisat2` or `rsubread`
explicitly when reproducing an analysis. All compute uses one thread.

The synthetic alignment workflow exercises Rsubread and native HISAT2 with
samtools and FastQC. HISAT2 runs in both DNA and RNA modes with spaces and an
apostrophe in the output path. These contiguous reads check the command and BAM
round trip; they do not test splice-junction detection. FastQC uses its supported
`--outdir` option, including with FastQC 0.13.

For equal quality-bin counts, ShortRead's heatmap can render an error panel.
The report uses the complete count table as points in that case, keeping every
cycle, quality score and count. CI also checks PDFs for this rendered error.
Review FastQC's module results separately from its command exit status. The
small synthetic fixture produces content warnings and a base-content failure;
successful processing does not mean every quality module passes.

HISAT2 alignment uses temporary symbolic links with fixed relative names because
its wrapper reconstructs shell commands. The links avoid copying large inputs
again and are removed on success or failure; inputs and alignment outputs remain
in the run directory. The native path requires symbolic-link support.

## Your own data

Supply both a FASTQ (plain or gzip) and an uncompressed reference FASTA:

```sh
FASTQ2BAM_FASTQ=/absolute/reads.fastq.gz \
FASTQ2BAM_REFERENCE=/absolute/reference.fasta \
FASTQ2BAM_TYPE=dna \
FASTQ2BAM_ALIGNER=rsubread \
FASTQ2BAM_RUN_DIR=/absolute/new-output-directory \
Rscript --vanilla 00_fastq_to_bam_pipeline.R
```

Use `FASTQ2BAM_TYPE=rna` for RNA alignment settings; the default is DNA, including
HISAT2's no-spliced-alignment setting. This tutorial processes one single-end
FASTQ with Phred+33 qualities and loads it into memory. It does not implement
paired-end analysis, production-scale streaming, quantification, variant calling
or clinical interpretation. Verify the organism, assembly, read type, quality
encoding and aligner settings before real analysis. Format FASTA lines below
1,000 bytes for Rsubread (80 bases per line is suitable).

Inputs are copied into a fresh run directory, with MD5 fingerprints for local
integrity tracking. Existing destinations and symlinks are refused. A partial
failure retains inputs, intermediate files and aligner logs; choose a fresh
output directory after correcting it. Missing/empty outputs and incomplete
indices stop the pipeline. No existing source data is replaced.

Missing R packages stop with an installation message. The legacy installer is
available only with `FASTQ2BAM_INSTALL_PACKAGES=1`; explicit installation above
is the reproducible path. Read the [German guide](ANLEITUNG_FASTQ_zu_BAM.md) and
[maintenance evidence](MAINTENANCE.md) for the tested paths and limitations.
