# FASTQ zu BAM — Schritt-für-Schritt-Anleitung in R

## Überblick

Diese Anleitung beschreibt den kompletten Workflow, um **FASTQ-Dateien** (rohe Sequenzierungsdaten) in **BAM-Dateien** (alignierte Reads) zu konvertieren — alles gesteuert aus R heraus.

### Was passiert bei FASTQ → BAM?

```
FASTQ (rohe Reads)
  │
  ├── 1. Qualitätskontrolle (QC)
  │      → Reads inspizieren, Qualitätsscores prüfen
  │
  ├── 2. Trimming (optional)
  │      → Adapter und Low-Quality-Basen entfernen
  │
  ├── 3. Alignment / Mapping
  │      → Reads gegen ein Referenzgenom alignieren
  │      → Ergebnis: SAM-Datei
  │
  ├── 4. SAM → BAM Konvertierung
  │      → Binärformat, komprimiert, schneller
  │
  ├── 5. Sortieren & Indexieren
  │      → BAM nach genomischer Position sortieren
  │      → Index (.bai) für schnellen Zugriff erstellen
  │
  └── 6. BAM in R einlesen & inspizieren
         → Rsamtools / GenomicAlignments
```

---

## Voraussetzungen

### Systemtools (müssen installiert sein)

| Tool | Zweck | Installation |
|------|-------|-------------|
| **samtools** | SAM/BAM-Manipulation | `sudo apt install samtools` (Linux) / `brew install samtools` (Mac) |
| **HISAT2** oder **STAR** | Aligner | `sudo apt install hisat2` / Download von [HISAT2-Website](http://daehwankimlab.github.io/hisat2/) |
| **FastQC** | QC-Reports | `sudo apt install fastqc` / Download von [Babraham](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
| **fastp** (optional) | Trimming | `sudo apt install fastp` / `conda install fastp` |

### Prüfen ob Tools verfügbar sind

Vor dem Start im Terminal prüfen:

```bash
samtools --version
hisat2 --version
fastqc --version
```

### R-Pakete

```r
# CRAN
install.packages(c("data.table", "ggplot2"))

# Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ShortRead", "Rsamtools", "GenomicAlignments", "Biostrings"))
```

---

## Schritt 1: Projektstruktur anlegen

Am Anfang jedes Projekts eine saubere Ordnerstruktur erstellen:

```r
# Projektordner
project_dir <- "fastq_to_bam_pipeline"
dir.create(project_dir, showWarnings = FALSE)

# Unterordner
dirs <- c("data/fastq", "data/reference", "data/bam",
          "results/qc", "results/alignment", "scripts")
for (d in file.path(project_dir, dirs)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
```

---

## Schritt 2: Beispiel-FASTQ herunterladen

Wir verwenden eine kleine Beispiel-FASTQ-Datei aus dem **ShortRead**-Paket (ca. 1 MB, E. coli Reads) — perfekt zum Testen.

```r
# Option A: ShortRead-Beispieldaten (bereits installiert)
library(ShortRead)
fastq_dir <- system.file("extdata", "E-MTAB-1147", package = "ShortRead")
fastq_files <- list.files(fastq_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE)
cat("Verfügbare FASTQ-Dateien:\n")
print(fastq_files)

# Option B: Von ENA herunterladen (kleiner SARS-CoV-2 Datensatz, ~5 MB)
# download.file(
#   url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR121/077/SRR12132977/SRR12132977.fastq.gz",
#   destfile = file.path(project_dir, "data/fastq/SRR12132977.fastq.gz"),
#   mode = "wb"
# )
```

---

## Schritt 3: FASTQ einlesen & inspizieren

```r
library(ShortRead)

# FASTQ einlesen
reads <- readFastq(fastq_files[1])

# Grundlegende Informationen
cat("Anzahl Reads:", length(reads), "\n")
cat("Read-Länge:", unique(width(reads)), "bp\n")

# Erste 3 Reads anschauen
sread(reads)[1:3]         # DNA-Sequenzen
quality(reads)[1:3]       # Qualitätsscores (ASCII)
id(reads)[1:3]            # Read-IDs

# Qualitätsscores als Zahlen
qual_matrix <- as(quality(reads), "matrix")
cat("Mittlere Qualität pro Read:", mean(rowMeans(qual_matrix, na.rm = TRUE)), "\n")
```

---

## Schritt 4: Qualitätskontrolle (QC)

### 4a: Mit ShortRead in R

```r
# QC-Report generieren
qa_result <- qa(fastq_files[1], type = "fastq")
report(qa_result, dest = file.path(project_dir, "results/qc"))
cat("QC-Report erstellt unter:", file.path(project_dir, "results/qc"), "\n")

# Eigene QC-Plots
library(ggplot2)

# Qualitätsverteilung pro Zyklus
qual_df <- data.frame(
  cycle = rep(seq_len(ncol(qual_matrix)), each = min(1000, nrow(qual_matrix))),
  quality = as.vector(qual_matrix[seq_len(min(1000, nrow(qual_matrix))), ])
)
qual_df <- qual_df[!is.na(qual_df$quality), ]

p_qual <- ggplot(qual_df, aes(x = factor(cycle), y = quality)) +
  geom_boxplot(outlier.size = 0.3, fill = "#4393C3") +
  scale_x_discrete(breaks = seq(1, max(qual_df$cycle), by = 5)) +
  labs(title = "Basenqualität pro Zyklus",
       x = "Zyklus (Position im Read)", y = "Phred-Score") +
  theme_minimal()

ggsave(file.path(project_dir, "results/qc/quality_per_cycle.pdf"), p_qual,
       width = 10, height = 5)
```

### 4b: Mit FastQC (extern)

```r
# FastQC über system() aufrufen
fastqc_path <- Sys.which("fastqc")
if (nchar(fastqc_path) > 0) {
  cmd <- paste("fastqc", fastq_files[1],
               "-o", file.path(project_dir, "results/qc"))
  system(cmd)
  cat("FastQC-Report erstellt.\n")
} else {
  cat("FastQC nicht gefunden. Bitte installieren.\n")
}
```

---

## Schritt 5: Referenzgenom vorbereiten

Für das Alignment brauchen wir ein Referenzgenom + Index.

```r
# Beispiel: E. coli K-12 Referenzgenom herunterladen (klein, ~5 MB)
ref_url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
ref_gz  <- file.path(project_dir, "data/reference/ecoli_k12.fna.gz")
ref_fa  <- file.path(project_dir, "data/reference/ecoli_k12.fna")

# Herunterladen
download.file(ref_url, ref_gz, mode = "wb")

# Entpacken
R.utils::gunzip(ref_gz, destname = ref_fa, remove = FALSE, overwrite = TRUE)
cat("Referenzgenom:", ref_fa, "\n")

# HISAT2-Index bauen
hisat2_path <- Sys.which("hisat2-build")
if (nchar(hisat2_path) > 0) {
  index_prefix <- file.path(project_dir, "data/reference/ecoli_k12_idx")
  cmd <- paste("hisat2-build", ref_fa, index_prefix)
  system(cmd)
  cat("HISAT2-Index erstellt.\n")
} else {
  cat("HISAT2 nicht gefunden. Alternative: Rsubread::buildindex()\n")
}
```

### Alternative: Index mit Rsubread (rein in R)

```r
if (requireNamespace("Rsubread", quietly = TRUE)) {
  Rsubread::buildindex(
    basename = file.path(project_dir, "data/reference/ecoli_k12_rsubread"),
    reference = ref_fa
  )
  cat("Rsubread-Index erstellt.\n")
}
```

---

## Schritt 6: Alignment (FASTQ → SAM/BAM)

### 6a: Mit HISAT2 (über system())

```r
# Pfade definieren
index_prefix <- file.path(project_dir, "data/reference/ecoli_k12_idx")
sam_file     <- file.path(project_dir, "results/alignment/aligned.sam")
bam_file     <- file.path(project_dir, "data/bam/aligned.bam")

# HISAT2 aufrufen
hisat2_cmd <- paste(
  "hisat2",
  "-x", index_prefix,                    # Index-Prefix
  "-U", fastq_files[1],                  # Single-End Reads
  "-S", sam_file,                         # Ausgabe: SAM
  "--threads", 4,                         # Threads
  "--summary-file", file.path(project_dir, "results/alignment/hisat2_summary.txt")
)

cat("Starte Alignment...\n")
exit_code <- system(hisat2_cmd)

if (exit_code == 0) {
  cat("Alignment erfolgreich. SAM-Datei:", sam_file, "\n")
} else {
  cat("Fehler beim Alignment. Exit-Code:", exit_code, "\n")
}
```

### 6b: Alternative mit Rsubread (komplett in R)

```r
if (requireNamespace("Rsubread", quietly = TRUE)) {
  index_base <- file.path(project_dir, "data/reference/ecoli_k12_rsubread")
  bam_rsubread <- file.path(project_dir, "data/bam/aligned_rsubread.bam")

  Rsubread::align(
    index = index_base,
    readfile1 = fastq_files[1],
    output_file = bam_rsubread,
    nthreads = 4
  )
  cat("Rsubread-Alignment fertig:", bam_rsubread, "\n")
}
```

---

## Schritt 7: SAM → BAM konvertieren, sortieren, indexieren

```r
# SAM → BAM
bam_unsorted <- file.path(project_dir, "data/bam/aligned_unsorted.bam")
cmd_convert  <- paste("samtools view -bS", sam_file, ">", bam_unsorted)
system(cmd_convert)

# BAM sortieren
bam_sorted <- file.path(project_dir, "data/bam/aligned_sorted.bam")
cmd_sort   <- paste("samtools sort", bam_unsorted, "-o", bam_sorted, "-@ 4")
system(cmd_sort)

# BAM indexieren
cmd_index <- paste("samtools index", bam_sorted)
system(cmd_index)

cat("Fertig!\n")
cat("BAM-Datei:  ", bam_sorted, "\n")
cat("Index-Datei:", paste0(bam_sorted, ".bai"), "\n")

# Statistiken
cmd_flagstat <- paste("samtools flagstat", bam_sorted)
system(cmd_flagstat)

# Alternative: samtools via Rsamtools
if (requireNamespace("Rsamtools", quietly = TRUE)) {
  Rsamtools::sortBam(bam_unsorted, destination = gsub("\\.bam$", "", bam_sorted))
  Rsamtools::indexBam(bam_sorted)
}
```

---

## Schritt 8: BAM in R einlesen & inspizieren

```r
library(Rsamtools)

# BAM-Datei öffnen
bam <- BamFile(bam_sorted)
bam

# Header-Informationen
seqinfo(bam)

# Reads zählen
counts <- countBam(bam)
cat("Gesamtzahl Reads:", counts$records, "\n")
cat("Nukleotide gesamt:", counts$nucleotides, "\n")

# Erste Reads einlesen
param <- ScanBamParam(
  what = c("qname", "flag", "rname", "pos", "mapq", "cigar", "seq", "qual"),
  which = GRanges()  # Alles einlesen
)

# Chunk-weise einlesen (speichereffizient)
yieldSize(bam) <- 1000
open(bam)
chunk <- scanBam(bam, param = ScanBamParam(
  what = c("qname", "rname", "pos", "mapq", "cigar")
))[[1]]
close(bam)

# Als data.frame
reads_df <- data.frame(
  read_id   = chunk$qname[1:10],
  chr       = chunk$rname[1:10],
  position  = chunk$pos[1:10],
  mapq      = chunk$mapq[1:10],
  cigar     = chunk$cigar[1:10],
  stringsAsFactors = FALSE
)
print(reads_df)
```

---

## Schritt 9: BAM visualisieren

```r
library(GenomicAlignments)
library(ggplot2)

# Alle Alignments einlesen
ga <- readGAlignments(bam_sorted)
cat("Alignierte Reads:", length(ga), "\n")

# Coverage berechnen
cov <- coverage(ga)

# Coverage-Plot für das erste Chromosom
chr_name <- names(cov)[1]
cov_vec  <- as.numeric(cov[[chr_name]])

# Nur einen Ausschnitt plotten (erste 10.000 bp)
plot_range <- seq_len(min(10000, length(cov_vec)))

cov_df <- data.frame(
  position = plot_range,
  coverage = cov_vec[plot_range]
)

p_cov <- ggplot(cov_df, aes(x = position, y = coverage)) +
  geom_area(fill = "#2166AC", alpha = 0.7) +
  labs(title = paste("Coverage —", chr_name),
       x = "Genomische Position (bp)", y = "Read-Tiefe") +
  theme_minimal()

ggsave(file.path(project_dir, "results/alignment/coverage_plot.pdf"),
       p_cov, width = 12, height = 4)
cat("Coverage-Plot gespeichert.\n")

# Mapping-Qualitätsverteilung
mapq_df <- data.frame(mapq = mcols(ga)$mapq)

p_mapq <- ggplot(mapq_df, aes(x = mapq)) +
  geom_histogram(bins = 50, fill = "#B2182B", color = "white") +
  labs(title = "Verteilung der Mapping-Qualität (MAPQ)",
       x = "MAPQ-Score", y = "Anzahl Reads") +
  theme_minimal()

ggsave(file.path(project_dir, "results/alignment/mapq_distribution.pdf"),
       p_mapq, width = 8, height = 5)
```

---

## Zusammenfassung

| Schritt | Was passiert | Tool/Paket |
|---------|-------------|------------|
| 1 | Projektstruktur anlegen | base R |
| 2 | FASTQ herunterladen | `download.file()` / ShortRead |
| 3 | FASTQ einlesen & inspizieren | ShortRead |
| 4 | Qualitätskontrolle | ShortRead / FastQC |
| 5 | Referenzgenom + Index | HISAT2 / Rsubread |
| 6 | Alignment (FASTQ → SAM) | HISAT2 / Rsubread |
| 7 | SAM → BAM → Sort → Index | samtools / Rsamtools |
| 8 | BAM einlesen | Rsamtools |
| 9 | Coverage & QC-Plots | GenomicAlignments / ggplot2 |

---

## Häufige Probleme

**"samtools: command not found"**
→ samtools ist nicht installiert oder nicht im PATH. Prüfe mit `Sys.which("samtools")` in R.

**"HISAT2 index not found"**
→ Der Index muss vor dem Alignment gebaut werden (Schritt 5). Prüfe ob die `.ht2`-Dateien existieren.

**"Error in readFastq: file not found"**
→ Pfad prüfen. `list.files()` zum Debuggen verwenden.

**BAM-Datei ist leer (0 Reads)**
→ Reads passen nicht zum Referenzgenom. Bei Testdaten darauf achten, dass Organismus übereinstimmt.

**Speicherprobleme bei großen Dateien**
→ `yieldSize` bei BamFile setzen und chunk-weise verarbeiten. Bei FASTQ: `FastqStreamer()` verwenden.

---

## Weiterführende Ressourcen

- [Bioconductor ShortRead Vignette](https://bioconductor.org/packages/ShortRead/)
- [Rsamtools Dokumentation](https://bioconductor.org/packages/Rsamtools/)
- [GenomicAlignments Vignette](https://bioconductor.org/packages/GenomicAlignments/)
- [HISAT2 Manual](http://daehwankimlab.github.io/hisat2/manual/)
- [samtools Documentation](http://www.htslib.org/doc/samtools.html)
