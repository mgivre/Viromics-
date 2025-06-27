#!/usr/bin/env python3
"""
Script de traitement des vOTUs - Pipeline complet
Usage: python votu_pipeline.py <input_tsv> [reference_fasta]
"""

import sys
import os
import subprocess
import argparse
from pathlib import Path
import logging
from datetime import datetime

def setup_logging(output_dir):
    """Configuration du logging"""
    log_file = output_dir / "pipeline.log"
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

def run_command(cmd, description, logger):
    """Exécute une commande système avec gestion d'erreur"""
    logger.info(f"Exécution: {description}")
    logger.info(f"Commande: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    
    try:
        result = subprocess.run(cmd, shell=isinstance(cmd, str), 
                              capture_output=True, text=True, check=True)
        if result.stdout:
            logger.info(f"Sortie: {result.stdout.strip()}")
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"Erreur lors de l'exécution de: {description}")
        logger.error(f"Code d'erreur: {e.returncode}")
        logger.error(f"Stderr: {e.stderr}")
        raise

def count_sequences(fasta_file):
    """Compte le nombre de séquences dans un fichier FASTA"""
    if not os.path.exists(fasta_file):
        return 0
    
    with open(fasta_file, 'r') as f:
        return sum(1 for line in f if line.startswith('>'))

def extract_viral_sequences(input_tsv, output_fasta, logger):
    """Extrait les séquences virales du fichier TSV"""
    logger.info("=== Extraction des séquences virales ===")
    
    viral_count = 0
    total_count = 0
    
    with open(input_tsv, 'r') as infile, open(output_fasta, 'w') as outfile:
        # Ignorer l'en-tête
        next(infile)
        
        for line in infile:
            total_count += 1
            fields = line.strip().split('\t')
            
            if len(fields) >= 3 and fields[2] == '1':  # Colonne 3 == 1 pour viral
                viral_count += 1
                seq_id = fields[0]
                sequence = fields[-1]  # Dernière colonne
                outfile.write(f">{seq_id}\n{sequence}\n")
    
    logger.info(f"Séquences totales dans TSV: {total_count}")
    logger.info(f"Séquences virales extraites: {viral_count}")
    if total_count > 0:
        percentage = (viral_count / total_count) * 100
        logger.info(f"Pourcentage viral: {percentage:.2f}%")
    
    return viral_count

def combine_sequences(viral_fasta, reference_fasta, output_fasta, logger):
    """Combine les séquences virales avec les séquences de référence"""
    logger.info("=== Combinaison avec les séquences de référence ===")
    
    # Copier les séquences virales
    with open(viral_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        outfile.write(infile.read())
    
    # Ajouter les séquences de référence
    with open(reference_fasta, 'r') as infile, open(output_fasta, 'a') as outfile:
        outfile.write(infile.read())
    
    viral_count = count_sequences(viral_fasta)
    ref_count = count_sequences(reference_fasta)
    total_count = count_sequences(output_fasta)
    
    logger.info(f"Séquences virales: {viral_count}")
    logger.info(f"Séquences de référence: {ref_count}")
    logger.info(f"Total combiné: {total_count}")
    
    return total_count

def run_prodigal(input_fasta, output_dir, logger):
    """Exécute Prodigal pour la prédiction de gènes"""
    logger.info("=== Gene calling avec Prodigal (mode meta) ===")
    
    output_prefix = output_dir / "vOTUs"
    gff_file = f"{output_prefix}.gff"
    faa_file = f"{output_prefix}.faa"
    fna_file = f"{output_prefix}.fna"
    
    cmd = [
        "prodigal",
        "-i", str(input_fasta),
        "-o", gff_file,
        "-a", faa_file,
        "-d", fna_file,
        "-p", "meta",
        "-f", "gff"
    ]
    
    run_command(cmd, "Prédiction de gènes avec Prodigal", logger)
    
    gene_count = count_sequences(faa_file)
    logger.info(f"Gènes prédits: {gene_count}")
    
    return faa_file, gene_count

def run_mmseqs_alignment(faa_file, output_dir, logger):
    """Exécute l'alignement all-vs-all avec MMseqs2"""
    logger.info("=== Alignement all-vs-all avec MMseqs2 ===")
    
    db_path = output_dir / "mmseqs_db"
    result_path = output_dir / "mmseqs_result"
    tmp_dir = output_dir / "tmp"
    alignment_output = output_dir / "vOTUs_alignment.tsv"
    
    # Créer la base de données
    logger.info("Création de la base de données MMseqs2...")
    run_command([
        "mmseqs", "createdb", str(faa_file), str(db_path)
    ], "Création base de données MMseqs2", logger)
    
    # Recherche all-vs-all
    logger.info("Recherche all-vs-all...")
    run_command([
        "mmseqs", "search", str(db_path), str(db_path), str(result_path), str(tmp_dir),
        "--threads", str(os.cpu_count()),
        "-e", "1e-5",
        "--max-seqs", "10000"
    ], "Recherche all-vs-all MMseqs2", logger)
    
    # Conversion en format tabulaire
    logger.info("Conversion en format tabulaire...")
    run_command([
        "mmseqs", "convertalis", str(db_path), str(db_path), str(result_path), str(alignment_output),
        "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
    ], "Conversion format tabulaire", logger)
    
    # Nettoyage
    if tmp_dir.exists():
        run_command(f"rm -rf {tmp_dir}", "Nettoyage fichiers temporaires", logger)
    
    # Compter les alignements
    if alignment_output.exists():
        with open(alignment_output, 'r') as f:
            alignment_count = sum(1 for _ in f)
        logger.info(f"Alignements trouvés: {alignment_count}")
    else:
        raise FileNotFoundError("Fichier d'alignement non généré")
    
    return str(alignment_output), alignment_count

def main():
    parser = argparse.ArgumentParser(description="Pipeline de traitement des vOTUs")
    parser.add_argument("input_tsv", help="Fichier TSV contenant les OTUs")
    parser.add_argument("--reference", "-r", 
                       default="mmseqs_output/OTUs_mmseqs.fna",
                       help="Fichier de séquences de référence")
    parser.add_argument("--output", "-o", 
                       default="votu_analysis",
                       help="Répertoire de sortie")
    
    args = parser.parse_args()
    
    # Vérifications
    if not os.path.exists(args.input_tsv):
        print(f"ERREUR: Fichier TSV introuvable: {args.input_tsv}")
        sys.exit(1)
    
    if not os.path.exists(args.reference):
        print(f"ERREUR: Fichier de référence introuvable: {args.reference}")
        sys.exit(1)
    
    # Configuration
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)
    
    logger = setup_logging(output_dir)
    
    try:
        logger.info("=== Début du pipeline vOTU ===")
        logger.info(f"Fichier TSV: {args.input_tsv}")
        logger.info(f"Fichier de référence: {args.reference}")
        logger.info(f"Répertoire de sortie: {output_dir}")
        
        # Étape 1: Extraction des vOTUs
        viral_fasta = output_dir / "viral_sequences.fna"
        viral_count = extract_viral_sequences(args.input_tsv, viral_fasta, logger)
        
        if viral_count == 0:
            logger.error("Aucune séquence virale trouvée")
            sys.exit(1)
        
        # Étape 2: Combinaison avec références
        combined_fasta = output_dir / "combined_sequences.fna"
        total_seqs = combine_sequences(viral_fasta, args.reference, combined_fasta, logger)
        
        # Étape 3: Gene calling
        faa_file, gene_count = run_prodigal(combined_fasta, output_dir, logger)
        
        # Étape 4: Alignement
        alignment_file, alignment_count = run_mmseqs_alignment(faa_file, output_dir, logger)
        
        # Résumé
        logger.info("=== Pipeline terminé avec succès ===")
        logger.info("Fichiers générés:")
        logger.info(f"  - Séquences virales: {viral_fasta}")
        logger.info(f"  - Séquences combinées: {combined_fasta}")
        logger.info(f"  - Protéines prédites: {faa_file}")
        logger.info(f"  - Alignements: {alignment_file}")
        
    except Exception as e:
        logger.error(f"Erreur dans le pipeline: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()