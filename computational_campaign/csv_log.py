from datetime import datetime
import csv
import pandas as pd


# Creates a CSV file for the computational campaign
def init_csv(base_dir):
    timestamp  = datetime.now().strftime('%d-%m-%Y_%H-%M-%S')
    csv_path   = base_dir / 'results' / f'{timestamp}.csv'

    # Ensure the 'results' folder exists
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    
    csv_file   = open(csv_path, 'w', newline='')
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow([
        'esecuzione_pk', 'istanza_pk', 'istanza',
        'positive', 'negative', 'neutre',
        'specie', 'strip', 'fori', 'cluster',
        'file_id', 'z', 'time', 'status',
        'sinergie', 'conflitti', 'thread', 'sn',
        'time_limit', 'solver', 'maximal'
    ])
    return csv_file, csv_writer, csv_path

# Writes row on the CSV containing the instance execution benchmarks
def write_row(csv_writer, esecuzione_pk, istanza_pk, instance,
              positive, negative, neutre, H, K, DIM_STRIP, cluster, file_id,
              z_val, wall_time, status_code, sinergie, conflitti,
              num_workers, constraint_mode, time_limit):
    sn_label = 'Hard' if constraint_mode == 'hard' else 'Soft'
    csv_writer.writerow([
        esecuzione_pk,
        istanza_pk,
        instance,
        positive,
        negative,
        neutre,
        H,
        K,
        DIM_STRIP,
        cluster,
        file_id,
        z_val,
        wall_time,
        status_code,
        sinergie,
        conflitti,
        num_workers,
        sn_label,
        time_limit,
        'OR-Tools',
        0
    ])

# Finalize computational campaign CSV and saves it
def finalize_csv(csv_file, csv_path):
    csv_file.close()
    df = pd.read_csv(csv_path)
    df['maximal'] = df.groupby('istanza_pk')['z'].transform('max')
    df.to_csv(csv_path, index=False)
    print(f"CSV file saved in: {csv_path}")

    timestamp = df['time_first_feasible'].name
    table.to_csv(f'results/aggregated_{timestamp}.csv', index=False)
