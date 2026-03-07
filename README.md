# human-antibiotic-dysbiosis-neutrophils
Code and analysis workflows for the manuscript on antibiotic-induced gut dysbiosis and neutrophil response in humans.

<br>

## Overview

This repository contains the computational workflows used to process and analyze sequencing datasets generated in the study:

### “Short-term antibiotic treatment alters human blood neutrophil activation and gene regulation through microbiota disruption.”

The project investigates how acute perturbation of the gut microbiota influences human neutrophil biology by integrating microbiome, transcriptomic, and chromatin accessibility profiling.

<br>

## Scope of this repository

This repository contains only the computational pipelines related to sequencing data processing and downstream analyses, including:

- 16S rRNA sequencing to characterize stool microbiota composition

- RNA-seq to profile transcriptional changes in circulating neutrophils

- ATAC-seq to assess chromatin accessibility and regulatory landscape remodeling in circulating neutrophils

Other experimental components of the study are not included in this repository.

All analyses are organized by data type and implemented using reproducible computational workflows.

<br>

## Experimental Design

Healthy male adult volunteers were treated with a short course of broad-spectrum antibiotics to induce transient microbiota disruption.

Samples collected:

| Sample type |	Timepoints | Analysis
| ----------- | -----------| ---------
| Stool | Pre / Post antibiotics |	16S rRNA microbiota profiling
| Peripheral blood neutrophils |	Pre / Post antibiotics (+/- LPS) |	RNA-seq
| Peripheral  blood neutrophils |	Pre / Post antibiotics |	ATAC-seq

<br>

## Data Types

## 16S rRNA Sequencing

Purpose: 
Characterize microbiota disruption and taxonomic shifts induced by antibiotic exposure.

Main analyses:
- microbiota diversity
- taxonomic composition
- differential abundance
- functional profile prediction

Pipeline:

- Raw read processing: `nf-core/ampliseq` (Nextflow)
- Downstream analysis: `R` 

Details available in:
``` 16S_rRNA/README.md ```

## RNA-seq

Purpose: Identify transcriptional programs in blood neutrophils following microbiota disruption.

Main analyses:
- differential gene expression
- pathway enrichment
- transcriptional signatures of inflammatory activation
- integration with chromatin accessibility data

Pipeline:

- Raw read processing: `nf-core/rnaseq` (Nextflow)
- Downstream analysis: `R`
  
Details available in:
``` RNA-seq/README.md ```

## ATAC-seq

Purpose: Characterize chromatin accessibility remodelling in blood neutrophils following microbiota disruption.

Main analyses:
- peak calling
- differential accessibility
- motif enrichment
- integration with transcriptomic changes

Pipeline:

- Raw read processing: `nf-core/atacseq` (Nextflow)
- Downstream analysis: `R`
  
Details available in:
``` ATAC-seq/README.md ```

<br>

## Data Availability

Raw sequencing data are available at:

NCBI Sequence Read Archive
Accession: PRJNA1402444 and PRJNA1402306

Processed count matrices and intermediate results required to reproduce the analyses are included in this repository where possible.

<br>
<br>

## Repository versions
### Manuscript status: under review
v1.0-submission – version used for initial manuscript submission
