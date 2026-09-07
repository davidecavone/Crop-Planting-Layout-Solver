from datetime import datetime
import csv
import pandas as pd

from csv_log import *

# Creates a CSV file for the computational campaign
def init_csv(base_dir):
    timestamp  = datetime.now().strftime('%d-%m-%Y_%H-%M-%S')
    csv_path   = base_dir / 'results' / f'{timestamp}.csv'
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

# Aggregate computational campaign CSV data in a table
def analyze(csv_path):
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"Error: fil '{csv_path}' not found.")
        sys.exit(1)

    df['60s']  = df['time_first_feasible'] <= 60
    df['120s'] = df['time_first_feasible'] <= 120
    df['240s'] = df['time_first_feasible'] <= 240
    df['480s'] = df['time_first_feasible'] <= 480

    table = df.groupby(['H', 'K', 'N']).agg(
        count=('instance', 'count'),
        s60=('60s', 'sum'),
        s120=('120s', 'sum'),
        s240=('240s', 'sum'),
        s480=('480s', 'sum')
    ).reset_index()

    timestamp = df['time_first_feasible'].name
    table.to_csv(f'results/aggregated_{timestamp}.csv', index=False)
