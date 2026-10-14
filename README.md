# Genome assembly and annotation of two freshwater sponges

Scripts used to assemble, annotate, curate, and validate the genomes of *Spongilla lacustris* (`GZXS202303`) and *Ephydatia fluviatilis* (`NJ151201`). Species-specific files and directories are named `S_lacustris` and `E_fluviatilis`.

## Input data

- Oxford Nanopore long reads and Illumina short reads for both species.
- Intermediate genome assemblies and sequencing alignments generated during assembly and polishing.
- Reference databases and protein datasets for contamination screening and genome annotation.

## Repository structure

- `workflows/01_assembly_polishing/`: assembly post-processing, polishing, and `purge_dups`/NextPolish wrappers.
- `workflows/02_contamination_curation/`: contamination screening, contig curation, and final assembly/annotation filtering.
- `workflows/03_repeat_annotation/`: repeat annotation summaries and final-assembly filtering.
- `workflows/04_structural_annotation/`: BRAKER configuration, transcript extraction, and GFF3/CDS/protein consistency checks.
- `workflows/05_functional_annotation/`: eggNOG-mapper, InterProScan, Swiss-Prot/DIAMOND, and GO/KEGG/COG annotation processing.
- `workflows/06_quality_control/`: BUSCO, Merqury, QUAST, and assembly-quality summaries.
- `workflows/07_figures/`: source-data preparation and manuscript figure scripts.
