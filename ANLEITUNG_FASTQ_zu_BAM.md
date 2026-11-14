# FASTQ zu BAM mit R

Diese Anleitung gehört zu `00_fastq_to_bam_pipeline.R`. Sie zeigt einen
nachvollziehbaren technischen Ablauf mit einzelnen Single-End-Reads. Für
biologische Schlussfolgerungen müssen Probe, Referenz und Parameter fachlich
zusammenpassen.

## Voraussetzungen

R ab 4.1 und eine dazu passende Bioconductor-Version. Lokal geprüft wurden R
4.6.1, ShortRead 1.70.0, Rsamtools 2.28.0 und Rsubread 2.26.0 auf macOS ARM64.
Die vollständigen tatsächlich geladenen Versionen stehen im Laufprotokoll.

```r
install.packages(c("BiocManager", "data.table", "ggplot2", "R.utils"))
BiocManager::install(c("ShortRead", "Rsamtools", "GenomicAlignments",
                      "Biostrings", "BiocParallel", "Rsubread"),
                    ask = FALSE, update = FALSE)
```

HISAT2 mit `hisat2-build` ist eine Alternative zu Rsubread. samtools und FastQC
sind optional. Fehlende Pakete werden standardmäßig gemeldet, nicht stillschweigend
installiert. Alle Rechenschritte verwenden einen Thread.

Der Workflow prüft neben Rsubread auch HISAT2 mit samtools und FastQC im DNA-
und RNA-Modus. Die Ausgabeordner enthalten dabei Leerzeichen und ein Apostroph.
Diese zusammenhängenden synthetischen Reads prüfen keine Spleißstellen.
FastQC verwendet die auch in Version 0.13 unterstützte Option `--outdir`.

## Demonstration ausführen

Im Repository-Verzeichnis:

```sh
Rscript --vanilla tests/test_commands.R
Rscript --vanilla tests/test_demo.R
Rscript --vanilla tests/test_quality_report.R
Rscript --vanilla 00_fastq_to_bam_pipeline.R
```

Ohne eigene Eingaben erzeugt der Lauf eine synthetische 20-kb-Referenz und 100
exakte Reads mit je 75 Basen, Seed 42 und Phred+33-Qualität 30. Sequenzdaten werden
nicht heruntergeladen. Für jeden Read sind Name, Referenz und Position bekannt;
der Lauf gleicht diese Wahrheit mit dem fertigen BAM ab.

Die Schritte sind:

1. Einen neuen Ausgabeordner unter `results/` anlegen; Eingaben und Hashwerte
   sowie den synthetischen Herkunftsvermerk speichern.
2. FASTQ mit ShortRead lesen, Sequenzlängen und Qualitäten prüfen, QC-Bericht und
   PDF-Grafiken erzeugen.
3. Referenzindex erstellen und die benötigten Indexdateien kontrollieren. Ein
   bloßes Logfile gilt nicht als fertiger Index.
4. Mit HISAT2 oder Rsubread ausrichten. DNA ist die Voreinstellung; für RNA
   `FASTQ2BAM_TYPE=rna` ausdrücklich setzen.
5. BAM sortieren und indexieren, Alignment und Coverage einlesen und darstellen.
6. Beim synthetischen Beispiel alle 100 Positionen und Read-Identitäten prüfen.

Der Lauf schreibt die konkreten Ausgabepfade und `sessionInfo()` ins Protokoll.
Das Beispiel belegt die Dateiverarbeitung, nicht die Genauigkeit eines Aligners
bei Fehlern, Spleißstellen, Varianten oder biologischen Proben.

Bei gleichen Häufigkeiten aller Qualitätsklassen kann die ShortRead-Heatmap
statt Daten einen Grafikfehler ausgeben. Dann stellt der Bericht die vollständige
Zähltabelle als Punkte dar; Zyklus, Qualität und Häufigkeit bleiben unverändert.
Der Workflow prüft die PDFs zusätzlich auf diese gerenderte Fehlermeldung.
FastQC-Modulergebnisse sind getrennt vom Exit-Code zu lesen: Das kleine
synthetische Beispiel enthält Inhaltswarnungen und einen fehlgeschlagenen
Basenzusammensetzungstest. Erfolgreiche Verarbeitung bedeutet keinen vollständigen
Qualitätsnachweis.

## Eigene Daten verwenden

FASTQ und Referenz müssen gemeinsam angegeben werden. Die FASTA-Referenz bleibt
unkomprimiert; Zeilen mit 80 Basen sind für beide Indexer geeignet. Der Ablauf
setzt Single-End und Phred+33 voraus und lädt die Reads in den Arbeitsspeicher.

```sh
FASTQ2BAM_FASTQ=/absoluter/pfad/reads.fastq.gz \
FASTQ2BAM_REFERENCE=/absoluter/pfad/referenz.fasta \
FASTQ2BAM_TYPE=dna \
FASTQ2BAM_ALIGNER=rsubread \
FASTQ2BAM_RUN_DIR=/absoluter/pfad/neuer-lauf \
Rscript --vanilla 00_fastq_to_bam_pipeline.R
```

Referenzorganismus, Assembly, Datentyp und Parameter vorab prüfen. Paired-End,
Quantifizierung, Variantenanalyse und produktionsgeeignete Verarbeitung großer
Dateien sind hier nicht implementiert. Die historischen Dateien unter `data/`
werden weder als biologisch passende Beispielkombination behauptet noch geändert.

## Fehler nachvollziehen

Vorhandene Ausgabeordner und Symlinks werden abgelehnt. Bei Fehlern bleiben
Eingaben, SAM/BAM-Zwischenstände und Logs erhalten; nach Korrektur einen neuen
Laufordner wählen. Externe Befehle müssen Exit-Code 0 liefern; ihre Argumente
werden getrennt und korrekt quotiert. Fehlende oder leere erwartete Dateien
stoppen den Ablauf. Rsubread kann einen Indexfehler nur ins Log schreiben;
deshalb wird zusätzlich die vollständige Indexausgabe geprüft.

Der HISAT2-Wrapper setzt intern erneut Shell-Befehle zusammen. Deshalb verwendet
das Alignment feste relative Namen über temporäre symbolische Verknüpfungen;
große Eingaben werden dabei nicht erneut kopiert. Die Verknüpfungen werden auch
bei Fehlern entfernt, während Eingaben und Alignment-Ausgaben erhalten bleiben.
Dieser Werkzeugpfad benötigt Unterstützung für symbolische Verknüpfungen.

Die tatsächlich beobachteten Ergebnisse und noch ungeprüften Werkzeugpfade stehen
in [MAINTENANCE.md](MAINTENANCE.md).
