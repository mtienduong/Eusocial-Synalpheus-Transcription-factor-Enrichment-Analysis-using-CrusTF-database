# Synalpheus Transcription factor Enrichment Analysis using CrusTF database

## Transcription Factor Family Analysis in *S. elizabethae*

This project identifies transcription factor (TF) families associated with differentially expressed transcripts in S. elizabethae using sequence matching and BLAST analysis based on CrusTF database, then compare their representation between all expressed transcripts (using the worker transcriptome) and differentially expressed genes in queen transcriptome. The goal is to explore which TF families may be enriched or specifically regulated in eusocial species.

CrusTF database: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-4305-2

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
    blastx -query S.elizabethae_transcriptome_assembly.fasta -db pro_db \
  -out blastx_results_S.elizabethae_transcriptome.txt \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -evalue 1e-5 -max_target_seqs 1 -num_threads 8
    ```
    Repeat this step for DEGs file

4.  **Match TFs to Families using Python script**:
    ```bash
    python3 S.elizabethae_transcriptome_TF_family_matches.py
    ```
    Repeat this step for DEGs file

5. **Run the full analysis on HPC (SLURM)**:
    ```bash
    sbatch run_tf_family_matcher.sh
    ```
    Repeat this step for DEGs file

6. **Compare TF Family Distributions**
      This comparison step is done in R. Currently not included in the respository.

