# Eusocial *Synalpheus* Transcription factor Enrichment Analysis using CrusTF database

## Transcription Factor Family Analysis in *S. elizabethae*

This project aims to identify transcription factor (TF) families associated with differentially expressed transcripts in S. elizabethae by performing sequence matching and BLAST analysis against the CrusTF database. We then compare the representation of these TF families between all expressed transcripts and those that are differentially expressed between queen and worker transcriptomes. The goal is to explore whether specific TF families are enriched or uniquely regulated in a eusocial context.

Prior to this analysis, we conducted an RNA-seq workflow using Galaxy, which included de novo transcriptome assembly and differential expression analysis to identify genes with significant expression differences between queen and worker castes.

CrusTF database: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-4305-2

Galaxy Tutorial: https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/full-de-novo/tutorial.html

---
## Project Summary

- Extracted sequences from the full transcriptome and a list of all DEGs
- Matched those sequences to a curated TF family database using `blastx`
- Identified which TF families each transcript belonged to
- Compared TF family distribution between:
  - **Transcriptome** 
  - **DEGs** 

This approach helps highlight whether certain TF families are more active (or suppressed) between queen and worker.

---

## Tools Used
- [Seqtk](https://github.com/lh3/seqtk): For extracting specific transcript sequences
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi): To match transcripts to the TF family database using `blastx`
- Python 3 with [Biopython](https://biopython.org/): For parsing and matching sequences
- OSC (Ohio Supercomputer Center) HPC system: For batch processing with SLURM

---

## Workflow
1. **Extract DEG Sequences**

Use `seqtk` to extract sequences for DEGs from the full transcriptome:
```bash
seqtk subseq S.elizabethae_transcriptome_assembly.fasta S.elizabethae_DEGs_sequences.txt > S.elizabethae_DEGs_sequences.fasta
```

2. **Concatenate Protein files and build a BLAST database**:
    ```bash
    cat Protein_Species/*.fasta > all_pro.fasta
    makeblastdb -in all_pro.fasta -dbtype prot -title "pro_db"
    ```
3. **Run `blastx`**:
    ```bash
    blastx -query S.elizabethae_transcriptome_assembly.fasta -db pro_db  -out blastx_results_S.elizabethae_transcriptome.txt  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -evalue 1e-5 -max_target_seqs 1 -num_threads 8
    ```
    Repeat this step for DEGs file

4.  **Match TFs to Families using Python script**:
    ```bash
    python3 S.elizabethae_transcriptome_TF_family_annotated.py
    ```
    Repeat this step for DEGs file

5. **Run the full analysis on HPC (SLURM)**:
    ```bash
    sbatch S.elizabethae_transcriptome_TF_family_annotated.sh
    ```
    Repeat this step for DEGs file

6. **Compare TF Family Distributions**
    This comparison step is being developed in R and will be updated shortly.

