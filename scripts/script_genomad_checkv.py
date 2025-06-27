import os
import subprocess
import pandas as pd
import logging
from multiprocessing import Pool
from datetime import datetime

# --- Configuration du logging ---
log_dir = "logs"
os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, f"analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger()

# --- Chemins ---
tsv_path = "benchmark.tsv"
output_dir = "output_analysis"
genomad_db = "/srv/scratch/yazidima/Maha/genomad_new_db/genomad_db"
checkv_db = "/srv/scratch/givreex/checkv-db-v1.5"
status_file = os.path.join(output_dir, "processing_status.tsv")  # Fichier pour suivre le statut des échantillons

# --- Préparation des répertoires ---
os.makedirs(output_dir, exist_ok=True)
fasta_dir = os.path.join(output_dir, "fasta")
os.makedirs(fasta_dir, exist_ok=True)

# --- Lecture et initialisation du fichier de statut ---
def initialize_status_file():
    if not os.path.exists(status_file):
        with open(status_file, 'w') as f:
            f.write("sequence_id\tgeNomad_status\tcheckV_status\n")
        return {}
    else:
        try:
            status_df = pd.read_csv(status_file, sep='\t')
            # Créer un dictionnaire pour un accès rapide au statut
            status_dict = {}
            for _, row in status_df.iterrows():
                status_dict[row['sequence_id']] = {
                    'geNomad_status': row['geNomad_status'],
                    'checkV_status': row['checkV_status']
                }
            return status_dict
        except Exception as e:
            logger.error(f"Erreur lors de la lecture du fichier de statut: {e}")
            return {}

# --- Mise à jour du fichier de statut ---
def update_status(seq_id, genomad_status=None, checkv_status=None):
    try:
        if os.path.exists(status_file):
            status_df = pd.read_csv(status_file, sep='\t')
        else:
            status_df = pd.DataFrame(columns=["sequence_id", "geNomad_status", "checkV_status"])
        
        # Vérifier si la séquence existe déjà dans le fichier de statut
        if seq_id in status_df['sequence_id'].values:
            idx = status_df.index[status_df['sequence_id'] == seq_id].tolist()[0]
            if genomad_status is not None:
                status_df.at[idx, 'geNomad_status'] = genomad_status
            if checkv_status is not None:
                status_df.at[idx, 'checkV_status'] = checkv_status
        else:
            new_row = {'sequence_id': seq_id, 
                      'geNomad_status': genomad_status if genomad_status is not None else 'pending',
                      'checkV_status': checkv_status if checkv_status is not None else 'pending'}
            status_df = pd.concat([status_df, pd.DataFrame([new_row])], ignore_index=True)
        
        # Sauvegarder le fichier de statut
        status_df.to_csv(status_file, sep='\t', index=False)
    except Exception as e:
        logger.error(f"Erreur lors de la mise à jour du fichier de statut pour {seq_id}: {e}")

# --- Lecture du TSV ---
def load_sequences():
    try:
        logger.info(f"Lecture du fichier {tsv_path}")
        df = pd.read_csv(tsv_path, sep='\t', header=None, low_memory=False)
        seq_col = df.columns[-1]
        sequences = df[[0, seq_col]].values.tolist()
        logger.info(f"Nombre de séquences chargées: {len(sequences)}")
        return sequences
    except Exception as e:
        logger.error(f"Erreur lors de la lecture du fichier TSV: {e}")
        return []

def check_genomad_results(genomad_output, seq_id):
    plasmid_summary = os.path.join(genomad_output, f"{seq_id}_plasmid_summary.tsv")
    virus_summary = os.path.join(genomad_output, f"{seq_id}_virus_summary.tsv")

    summaries = [plasmid_summary, virus_summary]
    found_files = [f for f in summaries if os.path.exists(f)]

    if not found_files:
        return "not_found"

    try:
        # Vérifie si au moins un fichier a plus d'une ligne (donc des hits)
        for file in found_files:
            with open(file, "r") as f:
                lines = f.readlines()
                if len(lines) > 1:
                    return "results_found"
        return "no_results"
    except Exception as e:
        logger.error(f"Erreur lors de la lecture des fichiers summary pour {seq_id}: {e}")
        return "error"
    
def check_checkv_results(checkv_output):
    """Vérifie si les résultats de CheckV sont complets et valides"""
    quality_file = os.path.join(checkv_output, "quality_summary.tsv")
    
    if not os.path.exists(quality_file):
        return False
    
    # Vérifier si les fichiers ne sont pas vides
    try:
        if os.path.exists(quality_file) and os.path.getsize(quality_file) > 0:
            return True
        return False
    except Exception as e:
        logger.error(f"Erreur lors de la vérification des résultats CheckV: {e}")
        return False

def process_sequence(entry):
    seq_id, sequence = entry

    fasta_path = os.path.join(fasta_dir, f"{seq_id}.fasta")
    genomad_output = os.path.join(output_dir, f"{seq_id}_genomad")
    checkv_output = os.path.join(output_dir, f"{seq_id}_checkv")

    # Créer FASTA si nécessaire
    if not os.path.exists(fasta_path):
        with open(fasta_path, "w") as f:
            f.write(f">{seq_id}\n{sequence}\n")

    genomad_status = "error"
    checkv_status = "error"

    # ---- geNomad ----
    genomad_result_status = check_genomad_results(genomad_output, seq_id)

    if genomad_result_status == "results_found":
        logger.info(f"[✔] Résultats geNomad valides avec détections pour {seq_id}")
        genomad_status = 'completed'

    elif genomad_result_status == "no_results":
        logger.info(f"[✓] geNomad terminé sans détection pour {seq_id}")
        genomad_status = 'completed_no_hits'

    elif genomad_result_status == "not_found":
        logger.info(f"[⚙] Lancement de geNomad pour {seq_id}")
        os.makedirs(genomad_output, exist_ok=True)
        try:
            result = subprocess.run(
                ["genomad", "end-to-end", fasta_path, genomad_output, genomad_db],
                check=False,
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                logger.error(f"[✗] Erreur geNomad pour {seq_id}")
                logger.error(f"↳ Code: {result.returncode}")
                logger.error(f"↳ stdout: {result.stdout}")
                logger.error(f"↳ stderr: {result.stderr}")
                genomad_status = 'failed'
            else:
                logger.debug(f"[geNomad stdout pour {seq_id}]:\n{result.stdout}")
                logger.debug(f"[geNomad stderr pour {seq_id}]:\n{result.stderr}")

                new_result = check_genomad_results(genomad_output, seq_id)
                if new_result == "results_found":
                    genomad_status = 'completed'
                elif new_result == "no_results":
                    genomad_status = 'completed_no_hits'
                else:
                    genomad_status = 'incomplete'
        except Exception as e:
            logger.error(f"[✗] Exception geNomad pour {seq_id}: {e}")
            genomad_status = 'error'
    else:
        genomad_status = 'error'

    # ---- CheckV ----
    try:
        if check_checkv_results(checkv_output):
            logger.info(f"[✔] Résultats CheckV valides pour {seq_id}")
            checkv_status = 'completed'
        else:
            logger.info(f"[⚙] Lancement CheckV pour {seq_id}")
            os.makedirs(checkv_output, exist_ok=True)
            result = subprocess.run(
                ["checkv", "end_to_end", fasta_path, checkv_output, "-d", checkv_db],
                check=False,
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                logger.error(f"[✗] Erreur CheckV pour {seq_id}")
                logger.error(f"↳ Code: {result.returncode}")
                logger.error(f"↳ stdout: {result.stdout}")
                logger.error(f"↳ stderr: {result.stderr}")
                checkv_status = 'failed'
            else:
                if check_checkv_results(checkv_output):
                    logger.info(f"[✓] CheckV terminé avec succès pour {seq_id}")
                    checkv_status = 'completed'
                else:
                    checkv_status = 'incomplete'
    except Exception as e:
        logger.error(f"[✗] Exception CheckV pour {seq_id}: {e}")
        checkv_status = 'error'

    return {
        'sequence_id': seq_id,
        'genomad_status': genomad_status,
        'checkv_status': checkv_status
    }

# --- MAIN ---
if __name__ == "__main__":
    sequences = load_sequences()
    if not sequences:
        logger.error("Aucune séquence à traiter. Arrêt du script.")
        exit(1)

    initialize_status_file()

    results = []
    with Pool(10) as pool:
        results = pool.map(process_sequence, sequences)

    # Mise à jour centralisée des statuts
    for result in results:
        update_status(
            result['sequence_id'],
            genomad_status=result['genomad_status'],
            checkv_status=result['checkv_status']
        )

    logger.info("Traitement terminé!")
