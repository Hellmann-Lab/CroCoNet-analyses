#!/bin/bash


#### hg38 ####

# add plasmid FASTA to hg38 FASTA
cp GRCh38.primary_assembly.genome.fa GRCh38.primary_assembly.genome.withdCas9.fa
cat pAAVS1_TetOn_zKRAB_dCas9_P2A_mCherry_insert.fa >> GRCh38.primary_assembly.genome.withdCas9.fa
grep ">" GRCh38.primary_assembly.genome.withdCas9.fa

# remove version suffix from transcript, gene, and exon IDs to match IDs in the liftoff annotation
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat gencode.v32.primary_assembly.annotation.gtf \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > gencode.v32.primary_assembly.annotation.modified.gtf

# Construct the gene ID allowlist
BIOTYPE_PATTERN=\
"(protein_coding|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
PAR_PATTERN="tag \"PAR\""
cat gencode.v32.primary_assembly.annotation.modified.gtf \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "gene_allowlist"

grep -E "^#" gencode.v32.primary_assembly.annotation.modified.gtf > gencode.v32.primary_assembly.annotation.prefiltered.gtf

# Filter to the gene allowlist
grep -Ff "gene_allowlist" gencode.v32.primary_assembly.annotation.modified.gtf \
    >> gencode.v32.primary_assembly.annotation.prefiltered.gtf

# filter GTF
cellranger mkgtf gencode.v32.primary_assembly.annotation.prefiltered.gtf gencode.v32.primary_assembly.annotation.filtered.gtf --attribute=gene_biotype:protein_coding

# add plasmid GTF to hg38 GTF
cp gencode.v32.primary_assembly.annotation.filtered.gtf gencode.v32.primary_assembly.annotation.filtered.withdCas9.gtf
cat pAAVS1_TetOn_zKRAB_dCas9_P2A_mCherry_insert.gtf >> gencode.v32.primary_assembly.annotation.filtered.withdCas9.gtf
tail gencode.v32.primary_assembly.annotation.filtered.withdCas9.gtf

# make reference
cellranger mkref --genome=GRCh38_withdCas9 \
--fasta=GRCh38.primary_assembly.genome.withdCas9.fa \
--genes=gencode.v32.primary_assembly.annotation.filtered.withdCas9.gtf

#### macFas6 ####

# add plasmid FASTA to hg38 FASTA
cp macFas6.fa macFas6.withdCas9.fa
cat pC_AAVS1_TetOn_zKRAB_dCas9_P2A_mCherry_insert.fa >> macFas6.withdCas9.fa
grep ">" macFas6.withdCas9.fa

# filter GTF
cellranger mkgtf macFas6_liftoff.gtf macFas6_liftoff.filtered.gtf --attribute=gene_biotype:protein_coding

# add plasmid GTF to hg38 GTF
cp macFas6_liftoff.filtered.gtf macFas6_liftoff.filtered.withdCas9.gtf
cat pC_AAVS1_TetOn_zKRAB_dCas9_P2A_mCherry_insert.gtf >> macFas6_liftoff.filtered.withdCas9.gtf
tail macFas6_liftoff.filtered.withdCas9.gtf

# make reference
cellranger mkref --genome=macFas6_withdCas9 \
--fasta=macFas6.withdCas9.fa \
--genes=macFas6_liftoff.filtered.withdCas9.gtf