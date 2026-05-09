import pandas as pd
import sys

def analyze(csv_path):
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"Errore: file '{csv_path}' non trovato.")
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