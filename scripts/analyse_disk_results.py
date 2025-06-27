#!/usr/bin/env python3
import os
import pandas as pd
import logging
from collections import defaultdict
import glob
import re

# Configuration du logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s'
)
logger = logging.getLogger()

# Configuration des chemins
tsv_path = "benchmark.tsv"
output_dir = "output_analysis"
status_file = os.path.join(output_dir, "processing_status.tsv")

def load_expected_sequences():
    """Charge les séquences qui devraient être traitées"""
    try:
        logger.info(f"Lecture du fichier {tsv_path}")
        df = pd.read_csv(tsv_path, sep='\t', header=None, low_memory=False)
        seq_ids = df[0].tolist()  # Première colonne = IDs
        logger.info(f"Nombre de séquences attendues: {len(seq_ids)}")
        return seq_ids
    except Exception as e:
        logger.error(f"Erreur lors de la lecture du fichier TSV: {e}")
        return []

def analyze_genomad_on_disk(seq_id):
    """Analyse les résultats geNomad directement sur disque"""
    # Format du dossier: seq{id}_genomad
    genomad_output = os.path.join(output_dir, f"seq{seq_id}_genomad")
    
    if not os.path.exists(genomad_output):
        return "not_started"
    
    # Chercher tous les fichiers de sortie possibles
    all_files = glob.glob(os.path.join(genomad_output, "*"))
    
    if not all_files:
        return "directory_empty"
    
    # Fichiers de résultats geNomad typiques
    plasmid_summary = os.path.join(genomad_output, f"seq{seq_id}_plasmid_summary.tsv")
    virus_summary = os.path.join(genomad_output, f"seq{seq_id}_virus_summary.tsv")
    
    # Autres fichiers possibles (geNomad peut utiliser différents noms)
    possible_files = [
        os.path.join(genomad_output, f"seq{seq_id}_summary", f"seq{seq_id}_plasmid_summary.tsv"),
        os.path.join(genomad_output, f"seq{seq_id}_summary", f"seq{seq_id}_virus_summary.tsv"),
        os.path.join(genomad_output, "plasmid_summary.tsv"),
        os.path.join(genomad_output, "virus_summary.tsv")
    ]
    
    # Chercher aussi dans les sous-dossiers
    for root, dirs, files in os.walk(genomad_output):
        for file in files:
            if "plasmid_summary" in file or "virus_summary" in file:
                possible_files.append(os.path.join(root, file))
    
    # Vérifier la présence des fichiers summary
    found_summaries = []
    for file_path in [plasmid_summary, virus_summary] + possible_files:
        if os.path.exists(file_path):
            found_summaries.append(file_path)
    
    if not found_summaries:
        # Peut-être en cours de traitement - chercher d'autres indices
        log_files = [f for f in all_files if f.endswith('.log')]
        tmp_files = [f for f in all_files if 'tmp' in f.lower() or f.endswith('.tmp')]
        
        if log_files or tmp_files:
            return "in_progress_maybe"
        else:
            return "incomplete_no_summary"
    
    # Analyser le contenu des fichiers summary trouvés
    results_found = False
    try:
        for file_path in found_summaries:
            if os.path.getsize(file_path) > 0:  # Fichier non vide
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    if len(lines) > 1:  # Plus qu'un header
                        results_found = True
                        break
        
        if results_found:
            return "completed_with_hits"
        else:
            return "completed_no_hits"
            
    except Exception as e:
        return f"error_reading: {str(e)[:50]}"

def analyze_checkv_on_disk(seq_id):
    """Analyse les résultats CheckV directement sur disque"""
    # Format du dossier: seq{id}_checkv
    checkv_output = os.path.join(output_dir, f"seq{seq_id}_checkv")
    
    if not os.path.exists(checkv_output):
        return "not_started"
    
    # Chercher tous les fichiers
    all_files = glob.glob(os.path.join(checkv_output, "*"))
    
    if not all_files:
        return "directory_empty"
    
    # Fichier de qualité CheckV
    quality_file = os.path.join(checkv_output, "quality_summary.tsv")
    
    # Autres emplacements possibles
    possible_quality_files = [
        os.path.join(checkv_output, f"seq{seq_id}_quality_summary.tsv"),
        os.path.join(checkv_output, "checkv_quality_summary.tsv"),
    ]
    
    # Chercher aussi dans les sous-dossiers
    for root, dirs, files in os.walk(checkv_output):
        for file in files:
            if "quality_summary" in file:
                possible_quality_files.append(os.path.join(root, file))
    
    # Trouver le fichier de qualité
    quality_file_found = None
    for qf in [quality_file] + possible_quality_files:
        if os.path.exists(qf):
            quality_file_found = qf
            break
    
    if not quality_file_found:
        # Peut-être en cours - chercher d'autres indices
        log_files = [f for f in all_files if f.endswith('.log')]
        tmp_files = [f for f in all_files if 'tmp' in f.lower() or f.endswith('.tmp')]
        
        if log_files or tmp_files:
            return "in_progress_maybe"
        else:
            return "incomplete_no_quality"
    
    try:
        if os.path.getsize(quality_file_found) == 0:
            return "quality_file_empty"
        
        # Vérifier le contenu
        with open(quality_file_found, 'r') as f:
            lines = f.readlines()
            if len(lines) <= 1:  # Seulement header ou vide
                return "completed_no_results"
            else:
                return "completed_with_results"
                
    except Exception as e:
        return f"error_reading: {str(e)[:50]}"

def scan_existing_directories():
    """Scanne les dossiers existants pour identifier les séquences déjà traitées"""
    logger.info("Scan des dossiers existants...")
    
    existing_sequences = set()
    
    if not os.path.exists(output_dir):
        logger.warning(f"Le dossier {output_dir} n'existe pas")
        return existing_sequences

    for item in os.listdir(output_dir):
        item_path = os.path.join(output_dir, item)
        if os.path.isdir(item_path):
            match = re.match(r"seq(\d+)_(genomad|checkv)", item)
            if match:
                seq_id = match.group(1)
                existing_sequences.add(seq_id)
    
    logger.info(f"Trouvé {len(existing_sequences)} séquences avec des dossiers de traitement")
    return existing_sequences

def create_status_from_disk():
    """Crée le fichier de statut basé sur l'analyse du disque"""
    logger.info("=== CRÉATION DU FICHIER DE STATUT DEPUIS LE DISQUE ===")
    
    expected_sequences = load_expected_sequences()
    if not expected_sequences:
        logger.error("Impossible de charger les séquences attendues")
        return
    
    # Scanner aussi les dossiers existants
    existing_sequences = scan_existing_directories()
    
    # Combiner les séquences attendues et existantes
    all_sequences = set(map(str, expected_sequences)) | existing_sequences
    all_sequences = sorted(list(all_sequences))  # Trier pour un traitement ordonné
    
    logger.info(f"Total de séquences à analyser: {len(all_sequences)}")
    logger.info(f"  - Séquences attendues: {len(expected_sequences)}")
    logger.info(f"  - Séquences avec dossiers existants: {len(existing_sequences)}")
    
    # Analyser chaque séquence
    results = []
    
    for i, seq_id in enumerate(all_sequences):
        if i % 100 == 0:
            logger.info(f"Progression: {i}/{len(all_sequences)} ({i/len(all_sequences)*100:.1f}%)")
        
        genomad_status = analyze_genomad_on_disk(seq_id)
        checkv_status = analyze_checkv_on_disk(seq_id)
        
        # Mapper vers les statuts standards
        genomad_mapped = map_genomad_status(genomad_status)
        checkv_mapped = map_checkv_status(checkv_status)
        
        results.append({
            'sequence_id': seq_id,
            'geNomad_status': genomad_mapped,
            'checkV_status': checkv_mapped,
            'genomad_disk_detail': genomad_status,
            'checkv_disk_detail': checkv_status,
            'in_expected': str(seq_id) in map(str, expected_sequences),
            'has_genomad_dir': os.path.exists(os.path.join(output_dir, f"seq{seq_id}_genomad")),
            'has_checkv_dir': os.path.exists(os.path.join(output_dir, f"seq{seq_id}_checkv"))
        })
    
    # Créer le DataFrame et sauvegarder
    df = pd.DataFrame(results)
    
    # Sauvegarder le fichier de statut standard
    standard_df = df[['sequence_id', 'geNomad_status', 'checkV_status']]
    standard_df.to_csv(status_file, sep='\t', index=False)
    logger.info(f"Fichier de statut créé: {status_file}")
    
    # Sauvegarder aussi un fichier détaillé
    detailed_file = os.path.join(output_dir, "processing_status_detailed.tsv")
    df.to_csv(detailed_file, sep='\t', index=False)
    logger.info(f"Fichier détaillé créé: {detailed_file}")
    
    return df

def map_genomad_status(disk_status):
    """Mappe les statuts disque vers les statuts standards geNomad"""
    mapping = {
        'not_started': 'pending',
        'directory_empty': 'pending',
        'in_progress_maybe': 'running',
        'incomplete_no_summary': 'incomplete',
        'completed_with_hits': 'completed',
        'completed_no_hits': 'completed_no_hits'
    }
    
    if disk_status.startswith('error_reading'):
        return 'error'
    
    return mapping.get(disk_status, 'error')

def map_checkv_status(disk_status):
    """Mappe les statuts disque vers les statuts standards CheckV"""
    mapping = {
        'not_started': 'pending',
        'directory_empty': 'pending',
        'in_progress_maybe': 'running',
        'incomplete_no_quality': 'incomplete',
        'quality_file_empty': 'incomplete',
        'completed_no_results': 'completed',
        'completed_with_results': 'completed'
    }
    
    if disk_status.startswith('error_reading'):
        return 'error'
    
    return mapping.get(disk_status, 'error')

def analyze_and_report():
    """Analyse complète et génération de rapport"""
    logger.info("=== ANALYSE COMPLÈTE DES RÉSULTATS SUR DISQUE ===")
    
    # Créer le fichier de statut depuis le disque
    df = create_status_from_disk()
    if df is None:
        return
    
    # Calculer les statistiques
    total = len(df)
    
    # Statistiques geNomad
    genomad_stats = df['geNomad_status'].value_counts()
    genomad_completed = genomad_stats.get('completed', 0) + genomad_stats.get('completed_no_hits', 0)
    
    # Statistiques CheckV
    checkv_stats = df['checkV_status'].value_counts()
    checkv_completed = checkv_stats.get('completed', 0)
    
    # Séquences complètement terminées
    fully_completed = len(df[
        (df['geNomad_status'].isin(['completed', 'completed_no_hits'])) &
        (df['checkV_status'] == 'completed')
    ])
    
    # Séquences pas encore commencées
    not_started = len(df[
        (df['geNomad_status'] == 'pending') &
        (df['checkV_status'] == 'pending')
    ])
    
    # Séquences en cours
    in_progress = len(df[
        (df['geNomad_status'] == 'running') |
        (df['checkV_status'] == 'running')
    ])
    
    # Séquences avec des dossiers mais pas dans le fichier attendu
    unexpected_sequences = len(df[df['in_expected'] == False])
    
    # Affichage des résultats
    logger.info(f"\n=== RÉSULTATS DE L'ANALYSE ===")
    logger.info(f"Total des séquences analysées: {total}")
    logger.info(f"Séquences attendues dans benchmark.tsv: {len(df[df['in_expected'] == True])}")
    logger.info(f"Séquences avec dossiers mais non attendues: {unexpected_sequences}")
    logger.info(f"Complètement terminées: {fully_completed} ({fully_completed/total*100:.1f}%)")
    logger.info(f"Pas encore commencées: {not_started} ({not_started/total*100:.1f}%)")
    logger.info(f"En cours de traitement: {in_progress} ({in_progress/total*100:.1f}%)")
    
    logger.info(f"\n=== DÉTAILS GENOMAD ===")
    for status, count in genomad_stats.items():
        logger.info(f"  {status}: {count} ({count/total*100:.1f}%)")
    
    logger.info(f"\n=== DÉTAILS CHECKV ===")
    for status, count in checkv_stats.items():
        logger.info(f"  {status}: {count} ({count/total*100:.1f}%)")
    
    # Statuts disque détaillés
    logger.info(f"\n=== DÉTAILS TECHNIQUES (DISQUE) ===")
    logger.info("geNomad sur disque:")
    disk_genomad = df['genomad_disk_detail'].value_counts()
    for status, count in disk_genomad.items():
        logger.info(f"  {status}: {count}")
    
    logger.info("CheckV sur disque:")
    disk_checkv = df['checkv_disk_detail'].value_counts()
    for status, count in disk_checkv.items():
        logger.info(f"  {status}: {count}")
    
    # Recommandations
    logger.info(f"\n=== RECOMMANDATIONS ===")
    if fully_completed == len(df[df['in_expected'] == True]):
        logger.info("✅ Toutes les séquences attendues sont traitées!")
        logger.info("   Vous pouvez procéder à l'analyse des résultats.")
    elif fully_completed / len(df[df['in_expected'] == True]) > 0.95:
        logger.info("⚠️  Presque terminé!")
        remaining = len(df[df['in_expected'] == True]) - fully_completed
        logger.info(f"   Il reste {remaining} séquences attendues à traiter.")
        logger.info("   Vous pouvez relancer le script corrigé.")
    elif not_started > fully_completed:
        logger.info("🔄 Traitement en cours de démarrage")
        logger.info("   La majorité des séquences n'ont pas encore été traitées.")
        logger.info("   Relancez le script corrigé avec un nombre de processus adapté.")
    else:
        logger.info("🔄 Traitement en cours")
        logger.info("   Relancez le script corrigé pour terminer les séquences restantes.")
    
    if unexpected_sequences > 0:
        logger.info(f"\n⚠️  Trouvé {unexpected_sequences} séquences avec des dossiers mais non listées dans benchmark.tsv")
        logger.info("   Cela peut indiquer des restes de traitements précédents.")
    
    # Créer un rapport résumé
    report_file = os.path.join(output_dir, "disk_analysis_report.txt")
    with open(report_file, 'w') as f:
        f.write("=== RAPPORT D'ANALYSE DU DISQUE ===\n")
        f.write(f"Date: {pd.Timestamp.now()}\n\n")
        f.write(f"Total séquences analysées: {total}\n")
        f.write(f"Séquences attendues: {len(df[df['in_expected'] == True])}\n")
        f.write(f"Complètement terminées: {fully_completed} ({fully_completed/total*100:.1f}%)\n")
        f.write(f"Pas encore commencées: {not_started} ({not_started/total*100:.1f}%)\n")
        f.write(f"En cours: {in_progress} ({in_progress/total*100:.1f}%)\n")
        f.write(f"Séquences inattendues: {unexpected_sequences}\n\n")
        
        expected_completed = len(df[(df['in_expected'] == True) & 
                                  (df['geNomad_status'].isin(['completed', 'completed_no_hits'])) &
                                  (df['checkV_status'] == 'completed')])
        expected_total = len(df[df['in_expected'] == True])
        
        if expected_completed == expected_total:
            f.write("STATUT: ✅ TRAITEMENT COMPLET DES SÉQUENCES ATTENDUES\n")
        else:
            f.write(f"STATUT: 🔄 {expected_total - expected_completed} séquences attendues restantes\n")
    
    logger.info(f"Rapport sauvegardé: {report_file}")
    
    return {
        'total': total,
        'fully_completed': fully_completed,
        'not_started': not_started,
        'in_progress': in_progress,
        'unexpected_sequences': unexpected_sequences
    }

if __name__ == "__main__":
    try:
        stats = analyze_and_report()
        
        if stats:
            completion_rate = stats['fully_completed'] / stats['total'] * 100
            print(f"\n🎯 RÉSUMÉ FINAL: {completion_rate:.1f}% des séquences sont complètement traitées")
            
            if stats['unexpected_sequences'] > 0:
                print(f"⚠️  {stats['unexpected_sequences']} séquences inattendues trouvées")
            
    except KeyboardInterrupt:
        logger.info("Analyse interrompue par l'utilisateur")
    except Exception as e:
        logger.error(f"Erreur lors de l'analyse: {e}")
        import traceback
        traceback.print_exc()