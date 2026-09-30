#!/usr/bin/env bash
set -e

echo "=== STEP 2: Building Transcript-to-Gene Map & Filtering FASTA ==="

GTF="../ref/Mus_musculus.GRCm38.100.gtf"
FASTA="../ref/Mus_musculus.GRCm38.cdna.all.fa"
TX2GENE="../ref/tx2gene.tsv"
FILTERED_FASTA="../ref/Filtered_FASTA.fa"

# 1. Parse GTF file to extract transcript_id.version -> gene_id.version mapping
echo "-> Parsing GTF to extract transcript and gene IDs..."
awk '$3 == "transcript" {
    tx_id = ""; gene_id = ""; tx_ver = ""; gene_ver = "";
    
    if (match($0, /transcript_id "([^"]+)"/, a)) tx_id = a[1];
    if (match($0, /transcript_version "([^"]+)"/, a)) tx_ver = a[1];
    if (match($0, /gene_id "([^"]+)"/, a)) gene_id = a[1];
    if (match($0, /gene_version "([^"]+)"/, a)) gene_ver = a[1];
    
    if (tx_id != "" && gene_id != "") {
        full_tx = (tx_ver != "") ? tx_id "." tx_ver : tx_id;
        full_gene = (gene_ver != "") ? gene_id "." gene_ver : gene_id;
        print full_tx "\t" full_gene;
    }
}' ${GTF} | sort -u > ${TX2GENE}

echo "-> Created tx2gene map with $(wc -l < ${TX2GENE}) entries."

# 2. Extract transcript IDs list
cut -f1 ${TX2GENE} > ../ref/valid_transcripts.txt

# 3. Filter FASTA using seqkit
echo "-> Filtering cDNA FASTA to match annotated GTF transcripts..."
seqkit grep -f ../ref/valid_transcripts.txt ${FASTA} > ${FILTERED_FASTA}

echo "-> Filtered FASTA contains $(grep -c "^>" ${FILTERED_FASTA}) sequence entries."
echo "=== Reference preparation complete! ==="
