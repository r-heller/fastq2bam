###############################################################################
#
#  FASTQ → BAM Pipeline in R
#  ==========================
#
#  Komplette Pipeline: FASTQ einlesen → QC → Alignment → BAM → Inspektion
#
#  Voraussetzungen:
#    - R >= 4.1.0
#    - Systemtools: samtools, hisat2 (oder Rsubread als Alternative)
#    - Bioconductor: ShortRead, Rsamtools, GenomicAlignments
#
#  Autor: Raban Heller
#  Datum: 2026-04-09
#
###############################################################################


# ============================================================================
# 0. SETUP — Pakete & Konfiguration
# ============================================================================

# --- Hilfsfunktion: Paket laden, ggf. installieren --------------------------
.ensure_pkg <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

# --- Systemtool prüfen ------------------------------------------------------
.check_tool <- function(tool) {
  path <- Sys.which(tool)
  if (nchar(path) == 0) {
    message("[WARNUNG] '", tool, "' nicht gefunden. ",
            "Bitte installieren oder PATH pruefen.")
    return(FALSE)
  }
  message("[OK] ", tool, ": ", path)
  return(TRUE)
}

# --- CRAN-Pakete ------------------------------------------------------------
.ensure_pkg("data.table")
.ensure_pkg("ggplot2")
.ensure_pkg("R.utils")

# --- Bioconductor-Pakete ----------------------------------------------------
.ensure_pkg("ShortRead",          bioc = TRUE)
.ensure_pkg("Rsamtools",          bioc = TRUE)
.ensure_pkg("GenomicAlignments",  bioc = TRUE)
.ensure_pkg("Biostrings",         bioc = TRUE)

# --- Systemtools prüfen -----------------------------------------------------
message("\n=== Systemtools pruefen ===")
has_samtools <- .check_tool("samtools")
has_hisat2   <- .check_tool("hisat2")
has_fastqc   <- .check_tool("fastqc")
has_rsubread <- requireNamespace("Rsubread", quietly = TRUE)

if (!has_hisat2 && !has_rsubread) {
  message("[INFO] Weder HISAT2 noch Rsubread verfuegbar.")
  message("       Installiere Rsubread als Fallback...")
  .ensure_pkg("Rsubread", bioc = TRUE)
  has_rsubread <- requireNamespace("Rsubread", quietly = TRUE)
}

message("")


# ============================================================================
# 1. KONFIGURATION
# ============================================================================

# --- Projektpfade -----------------------------------------------------------
PROJECT_DIR <- getwd()  # Aktuelles Arbeitsverzeichnis = Projektordner

PATHS <- list(
  fastq     = file.path(PROJECT_DIR, "data", "fastq"),
  reference = file.path(PROJECT_DIR, "data", "reference"),
  bam       = file.path(PROJECT_DIR, "data", "bam"),
  qc        = file.path(PROJECT_DIR, "results", "qc"),
  alignment = file.path(PROJECT_DIR, "results", "alignment"),
  scripts   = file.path(PROJECT_DIR, "scripts")
)

# Ordner erstellen
invisible(lapply(PATHS, dir.create, recursive = TRUE, showWarnings = FALSE))

# --- Parameter --------------------------------------------------------------
THREADS  <- 4L      # Anzahl CPU-Threads
ALIGNER  <- "auto"  # "hisat2", "rsubread", oder "auto"

message("=== Projekt-Konfiguration ===")
message("Projektverzeichnis: ", PROJECT_DIR)
message("Threads:            ", THREADS)
message("Aligner:            ", ALIGNER)
message("")


# ============================================================================
# 2. BEISPIEL-FASTQ HERUNTERLADEN
# ============================================================================

message("=== Schritt 2: FASTQ-Daten bereitstellen ===")

# Wir verwenden die ShortRead-Beispieldaten (immer verfuegbar nach Installation)
fastq_source <- system.file("extdata", "E-MTAB-1147", package = "ShortRead")
fastq_files  <- list.files(fastq_source, pattern = "\\.fastq\\.gz$",
                           full.names = TRUE)

if (length(fastq_files) == 0) {
  stop("Keine FASTQ-Dateien im ShortRead-Paket gefunden.")
}

# Dateien ins Projektverzeichnis kopieren
for (f in fastq_files) {
  dest <- file.path(PATHS$fastq, basename(f))
  if (!file.exists(dest)) file.copy(f, dest)
}

# Arbeitskopien verwenden
fastq_files <- list.files(PATHS$fastq, pattern = "\\.fastq\\.gz$",
                          full.names = TRUE)

message("FASTQ-Dateien (", length(fastq_files), "):")
invisible(lapply(fastq_files, function(f) {
  sz <- file.info(f)$size
  message("  ", basename(f), " (", round(sz / 1024, 1), " KB)")
}))
message("")


# ============================================================================
# 3. FASTQ EINLESEN & INSPIZIEREN
# ============================================================================

message("=== Schritt 3: FASTQ inspizieren ===")

library(ShortRead)

# Erste Datei einlesen
reads <- readFastq(fastq_files[1])

message("Datei:         ", basename(fastq_files[1]))
message("Anzahl Reads:  ", length(reads))
message("Read-Laenge:   ", paste(range(width(reads)), collapse = "–"), " bp")

# Qualitaetsscores
qual_matrix <- as(quality(reads), "matrix")
mean_qual   <- mean(rowMeans(qual_matrix, na.rm = TRUE), na.rm = TRUE)
message("Mittl. Qualitaet: ", round(mean_qual, 1), " (Phred)")

# Erste 3 Reads zeigen
message("\nErste 3 Sequenzen:")
print(sread(reads)[1:3])

message("")


# ============================================================================
# 4. QUALITAETSKONTROLLE (QC)
# ============================================================================

message("=== Schritt 4: Qualitaetskontrolle ===")

library(ggplot2)

# --- 4a: ShortRead QC-Report ------------------------------------------------
# SerialParam vermeidet Serialisierungsfehler auf Windows
BiocParallel::register(BiocParallel::SerialParam())
qa_result <- qa(fastq_files[1], type = "fastq")
report(qa_result, dest = PATHS$qc)
message("ShortRead QC-Report: ", PATHS$qc)

# --- 4b: Qualitaet pro Zyklus plotten ---------------------------------------
n_sample <- min(2000L, nrow(qual_matrix))
sampled  <- qual_matrix[sample(nrow(qual_matrix), n_sample), ]

qual_df <- data.frame(
  cycle   = rep(seq_len(ncol(sampled)), each = n_sample),
  quality = as.vector(sampled)
)
qual_df <- qual_df[!is.na(qual_df$quality), ]

p_qual <- ggplot(qual_df, aes(x = factor(cycle), y = quality)) +
  geom_boxplot(outlier.size = 0.2, fill = "#4393C3", alpha = 0.8) +
  scale_x_discrete(breaks = seq(1, max(qual_df$cycle), by = 5)) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red", linewidth = 0.5) +
  annotate("text", x = 1, y = 21, label = "Q20", color = "red",
           hjust = 0, size = 3) +
  labs(title = "Basenqualitaet pro Zyklus",
       subtitle = basename(fastq_files[1]),
       x = "Zyklus (Position im Read)",
       y = "Phred-Score") +
  theme_minimal(base_size = 11)

ggsave(file.path(PATHS$qc, "quality_per_cycle.pdf"), p_qual,
       width = 10, height = 5)
message("Plot gespeichert: quality_per_cycle.pdf")

# --- 4c: GC-Content ---------------------------------------------------------
gc_content <- letterFrequency(sread(reads), "GC", as.prob = TRUE)[, 1]
gc_df <- data.frame(gc = gc_content)

p_gc <- ggplot(gc_df, aes(x = gc)) +
  geom_histogram(bins = 50, fill = "#1B7837", color = "white", alpha = 0.8) +
  labs(title = "GC-Content Verteilung",
       x = "GC-Anteil", y = "Anzahl Reads") +
  theme_minimal(base_size = 11)

ggsave(file.path(PATHS$qc, "gc_content.pdf"), p_gc, width = 8, height = 5)
message("Plot gespeichert: gc_content.pdf")

# --- 4d: FastQC (wenn verfuegbar) -------------------------------------------
if (has_fastqc) {
  cmd <- paste("fastqc", fastq_files[1], "-o", PATHS$qc, "-t", THREADS)
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  message("FastQC-Report erstellt.")
}

message("")


# ============================================================================
# 5. REFERENZGENOM VORBEREITEN
# ============================================================================

message("=== Schritt 5: Referenzgenom herunterladen & Index bauen ===")

# E. coli K-12 MG1655 (klein, ~5 MB — ideal zum Testen)
ref_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/",
  "GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
)
ref_gz <- file.path(PATHS$reference, "ecoli_k12.fna.gz")
ref_fa <- file.path(PATHS$reference, "ecoli_k12.fna")

if (!file.exists(ref_fa)) {
  message("Lade E. coli K-12 Referenzgenom herunter...")
  download.file(ref_url, ref_gz, mode = "wb", quiet = TRUE)
  R.utils::gunzip(ref_gz, destname = ref_fa, remove = FALSE, overwrite = TRUE)
  message("Referenzgenom: ", ref_fa)
}

# --- Index bauen ------------------------------------------------------------
index_prefix <- file.path(PATHS$reference, "ecoli_k12_idx")

use_rsubread <- FALSE

if (has_hisat2 && (ALIGNER %in% c("hisat2", "auto"))) {
  # HISAT2-Index
  idx_files <- list.files(PATHS$reference, pattern = "ecoli_k12_idx.*\\.ht2$")
  if (length(idx_files) == 0) {
    message("Baue HISAT2-Index...")
    cmd <- paste("hisat2-build -q", ref_fa, index_prefix)
    system(cmd, ignore.stdout = TRUE)
    message("HISAT2-Index fertig.")
  } else {
    message("HISAT2-Index bereits vorhanden.")
  }
} else if (has_rsubread) {
  use_rsubread <- TRUE
  rsubread_idx <- file.path(PATHS$reference, "ecoli_k12_rsubread")
  idx_files <- list.files(PATHS$reference, pattern = "ecoli_k12_rsubread")
  if (length(idx_files) == 0) {
    message("Baue Rsubread-Index...")
    Rsubread::buildindex(basename = rsubread_idx, reference = ref_fa)
    message("Rsubread-Index fertig.")
  }
} else {
  stop("Kein Aligner verfuegbar. Bitte HISAT2 oder Rsubread installieren.")
}

message("")


# ============================================================================
# 6. ALIGNMENT (FASTQ → SAM/BAM)
# ============================================================================

message("=== Schritt 6: Alignment ===")

sam_file <- file.path(PATHS$alignment, "aligned.sam")
bam_raw  <- file.path(PATHS$bam, "aligned_unsorted.bam")

if (!use_rsubread) {
  # --- HISAT2 -----------------------------------------------------------------
  message("Aligner: HISAT2")
  summary_file <- file.path(PATHS$alignment, "hisat2_summary.txt")

  hisat2_cmd <- paste(
    "hisat2",
    "-x", index_prefix,
    "-U", fastq_files[1],
    "-S", sam_file,
    "--threads", THREADS,
    "--summary-file", summary_file,
    "--no-unal"                         # Unalignierte Reads nicht ausgeben
  )

  t0 <- Sys.time()
  exit_code <- system(hisat2_cmd)
  dt <- difftime(Sys.time(), t0, units = "secs")

  if (exit_code == 0) {
    message("Alignment fertig (", round(dt, 1), " Sek.)")

    # Summary anzeigen
    if (file.exists(summary_file)) {
      message("\nHISAT2-Summary:")
      writeLines(readLines(summary_file))
    }
  } else {
    stop("HISAT2-Alignment fehlgeschlagen. Exit-Code: ", exit_code)
  }

} else {
  # --- Rsubread ---------------------------------------------------------------
  message("Aligner: Rsubread")
  bam_rsubread <- file.path(PATHS$bam, "aligned_rsubread")

  t0 <- Sys.time()
  Rsubread::align(
    index      = file.path(PATHS$reference, "ecoli_k12_rsubread"),
    readfile1  = fastq_files[1],
    output_file = paste0(bam_rsubread, ".BAM"),
    nthreads   = THREADS
  )
  dt <- difftime(Sys.time(), t0, units = "secs")

  message("Alignment fertig (", round(dt, 1), " Sek.)")

  # Rsubread erzeugt direkt BAM
  bam_raw <- paste0(bam_rsubread, ".BAM")
  sam_file <- NULL  # Kein SAM bei Rsubread
}

message("")


# ============================================================================
# 7. SAM → BAM → SORTIEREN → INDEXIEREN
# ============================================================================

message("=== Schritt 7: SAM -> BAM konvertieren, sortieren, indexieren ===")

bam_sorted <- file.path(PATHS$bam, "aligned_sorted.bam")
bam_index  <- paste0(bam_sorted, ".bai")

if (!is.null(sam_file) && file.exists(sam_file)) {
  # SAM → BAM
  if (has_samtools) {
    message("Konvertiere SAM -> BAM...")
    system(paste("samtools view -bS -@", THREADS, sam_file, ">", bam_raw))

    message("Sortiere BAM...")
    system(paste("samtools sort -@", THREADS, bam_raw, "-o", bam_sorted))

    message("Indexiere BAM...")
    system(paste("samtools index", bam_sorted))

    # Aufraumen: SAM und unsortierte BAM loeschen
    unlink(sam_file)
    unlink(bam_raw)

    message("SAM/BAM-Verarbeitung fertig (samtools).")

  } else if (requireNamespace("Rsamtools", quietly = TRUE)) {
    message("Verwende Rsamtools fuer SAM -> BAM...")
    # Manuell: samtools nicht vorhanden, aber Rsamtools kann sortieren/indexieren
    # Rsamtools braucht aber eine BAM-Datei als Input
    message("[WARNUNG] Ohne samtools kann SAM nicht direkt konvertiert werden.")
    message("          Bitte samtools installieren.")
  }

} else if (file.exists(bam_raw)) {
  # Rsubread hat direkt BAM erzeugt → sortieren & indexieren
  if (has_samtools) {
    system(paste("samtools sort -@", THREADS, bam_raw, "-o", bam_sorted))
    system(paste("samtools index", bam_sorted))
  } else {
    library(Rsamtools)
    sortBam(bam_raw, destination = sub("\\.bam$", "", bam_sorted))
    indexBam(bam_sorted)
  }
  message("BAM sortiert & indexiert.")
}

# --- Statistiken ------------------------------------------------------------
if (has_samtools && file.exists(bam_sorted)) {
  message("\nBAM-Statistiken (samtools flagstat):")
  system(paste("samtools flagstat", bam_sorted))
}

message("")


# ============================================================================
# 8. BAM IN R EINLESEN & INSPIZIEREN
# ============================================================================

message("=== Schritt 8: BAM in R einlesen ===")

library(Rsamtools)

if (file.exists(bam_sorted)) {
  bam <- BamFile(bam_sorted)

  # Grundinfos
  si <- seqinfo(bam)
  message("Referenzsequenzen im BAM:")
  print(si)

  # Reads zaehlen
  counts <- countBam(bam)
  message("\nGesamtzahl Reads:  ", counts$records)
  message("Nukleotide gesamt: ", counts$nucleotides)

  # Chunk-weise einlesen
  yieldSize(bam) <- 100
  open(bam)
  chunk <- scanBam(bam, param = ScanBamParam(
    what = c("qname", "rname", "pos", "mapq", "cigar")
  ))[[1]]
  close(bam)

  # Als data.frame
  n_show <- min(10, length(chunk$qname))
  reads_df <- data.frame(
    read_id  = chunk$qname[seq_len(n_show)],
    chr      = as.character(chunk$rname[seq_len(n_show)]),
    position = chunk$pos[seq_len(n_show)],
    mapq     = chunk$mapq[seq_len(n_show)],
    cigar    = chunk$cigar[seq_len(n_show)],
    stringsAsFactors = FALSE
  )

  message("\nErste ", n_show, " alignierte Reads:")
  print(reads_df)

} else {
  message("[WARNUNG] Keine sortierte BAM-Datei gefunden.")
}

message("")


# ============================================================================
# 9. COVERAGE & VISUALISIERUNG
# ============================================================================

message("=== Schritt 9: Coverage berechnen & visualisieren ===")

if (file.exists(bam_sorted)) {
  library(GenomicAlignments)
  library(ggplot2)

  # Alignments einlesen (mapq muss explizit angefordert werden)
  ga <- readGAlignments(bam_sorted, param = ScanBamParam(what = "mapq"))
  message("Alignierte Reads: ", length(ga))

  # Coverage berechnen
  cov <- coverage(ga)
  chr_name <- names(cov)[1]
  cov_vec  <- as.numeric(cov[[chr_name]])

  message("Chromosom:       ", chr_name)
  message("Genomlaenge:     ", length(cov_vec), " bp")
  message("Max. Coverage:   ", max(cov_vec), "x")
  message("Mittl. Coverage: ", round(mean(cov_vec), 2), "x")

  # --- Coverage-Plot ----------------------------------------------------------
  # Ausschnitt: Erste 50 kb oder gesamtes Genom wenn kleiner
  plot_end <- min(50000L, length(cov_vec))
  cov_df <- data.frame(
    position = seq_len(plot_end),
    coverage = cov_vec[seq_len(plot_end)]
  )

  p_cov <- ggplot(cov_df, aes(x = position, y = coverage)) +
    geom_area(fill = "#2166AC", alpha = 0.7) +
    labs(title = paste("Read-Coverage —", chr_name),
         subtitle = paste0("Positionen 1–", format(plot_end, big.mark = "."),
                           " | Mittl. Coverage: ",
                           round(mean(cov_vec), 1), "x"),
         x = "Genomische Position (bp)",
         y = "Read-Tiefe") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(PATHS$alignment, "coverage_plot.pdf"), p_cov,
         width = 12, height = 4)
  message("Coverage-Plot gespeichert.")

  # --- MAPQ-Verteilung --------------------------------------------------------
  mapq_vals <- mcols(ga)$mapq
  mapq_df <- data.frame(mapq = mapq_vals[!is.na(mapq_vals)])

  p_mapq <- ggplot(mapq_df, aes(x = mapq)) +
    geom_histogram(bins = 50, fill = "#B2182B", color = "white", alpha = 0.8) +
    labs(title = "Verteilung der Mapping-Qualitaet (MAPQ)",
         x = "MAPQ-Score",
         y = "Anzahl Reads") +
    theme_minimal(base_size = 11)

  ggsave(file.path(PATHS$alignment, "mapq_distribution.pdf"), p_mapq,
         width = 8, height = 5)
  message("MAPQ-Plot gespeichert.")

  # --- Read-Laengenverteilung -------------------------------------------------
  rlen_df <- data.frame(width = width(ga))

  p_rlen <- ggplot(rlen_df, aes(x = width)) +
    geom_histogram(bins = 50, fill = "#762A83", color = "white", alpha = 0.8) +
    labs(title = "Read-Laengenverteilung (aligniert)",
         x = "Alignierte Laenge (bp)",
         y = "Anzahl Reads") +
    theme_minimal(base_size = 11)

  ggsave(file.path(PATHS$alignment, "read_length_distribution.pdf"), p_rlen,
         width = 8, height = 5)
  message("Read-Laengen-Plot gespeichert.")
}

message("")


# ============================================================================
# FERTIG
# ============================================================================

message("=============================================")
message("  Pipeline abgeschlossen!")
message("=============================================")
message("")
message("Ergebnisse:")
message("  QC-Reports:    ", PATHS$qc)
message("  BAM-Datei:     ", bam_sorted)
message("  BAM-Index:     ", bam_index)
message("  Plots:         ", PATHS$alignment)
message("")
message("Naechste Schritte:")
message("  - Reads zaehlen:  Rsubread::featureCounts() / GenomicAlignments::summarizeOverlaps()")
message("  - DE-Analyse:     DESeq2 / edgeR / limma")
message("  - Visualisierung: IGV, bambamR")
message("")

# Session-Info
message("=== Session Info ===")
sessionInfo()
