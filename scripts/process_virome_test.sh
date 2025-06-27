#!/bin/bash

# ----------------------
# PARAMÈTRES UTILISATEUR
# ----------------------
R1="./PRJEB46943/ERR8081229_1.fastq.gz"
R2="./PRJEB46943/ERR8081229_2.fastq.gz"
SAMPLE="ERR8081229"
THREADS=4
ADAPTER_SEQ="CTGTCTCTTATACACATCT"

# ----------------------
# ÉTAPE 1 : TRIMMING AVEC FASTP
# ----------------------
echo "Étape 1 : trimming avec fastp..."

fastp \
  --in1 "$R1" \
  --in2 "$R2" \
  --out1 "${SAMPLE}_clean_R1.fastq.gz" \
  --out2 "${SAMPLE}_clean_R2.fastq.gz" \
  --adapter_sequence "$ADAPTER_SEQ" \
  --qualified_quality_phred 13 \
  --unqualified_percent_limit 10 \
  --length_required 32 \
  --trim_poly_g \
  --trim_poly_x \
  --detect_adapter_for_pe \
  --correction \
  --thread $THREADS \
  --html "${SAMPLE}_fastp.html" \
  --json "${SAMPLE}_fastp.json" \

# ----------------------
# ÉTAPE 2 : DEREPLICATION AVEC VSEARCH
# ----------------------
echo "Étape 2 : déduplication avec vsearch..."

zcat "${SAMPLE}_clean_R1.fastq.gz" | \
vsearch --derep_prefix - \
        --output "${SAMPLE}_derep.fasta" \
        --quiet

        
echo "Pipeline terminé !"
echo "Résultat : ${SAMPLE}_derep.fasta"

# ----------------------
# ÉTAPE 3 : ASSEMBLAGE AVEC SPADES
# ----------------------
echo "Étape 3 : assemblage avec spades..."

# Fichier en entrée : reads dédupliqués
INPUT="${SAMPLE}_derep.fasta"
OUTDIR="spades_derep_${SAMPLE}"

# Vérification
if [[ ! -f "$INPUT" ]]; then
    echo "Erreur : fichier d'entrée introuvable ($INPUT)"
    exit 1
fi

# Assemblage avec SPAdes en mode méta à partir de reads non appariés
echo "Assemblage avec SPAdes (entrée : reads dédupliqués non appariés)..."
spades.py \
  --meta \
  --only-assembler \
  -s "$INPUT" \
  -o "$OUTDIR" \
  -t 8 \
  -m 32

# Vérification du succès
if [[ -f "$OUTDIR/contigs.fasta" ]]; then
    echo "Assemblage terminé avec succès."
    echo "Résultat : $OUTDIR/contigs.fasta"
else
    echo "Erreur : assemblage échoué ou contigs manquants."
fi
