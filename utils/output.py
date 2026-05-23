import matplotlib.pyplot as plt
import matplotlib.patches as patches
from datetime import datetime
import csv

def print_solution(instance, solver, presence, start, end, size, HSI):
    print(f"\nIstanza: {instance}")
    print(f"Objective: {solver.objective_value}")
    for (h, s, i) in HSI:
        if solver.value(presence[h, s, i]):
            print(f"  z[{h},{s},{i}]: start={solver.value(start[h,s,i])}, "
                  f"end={solver.value(end[h,s,i])}, size={solver.value(size[h,s,i])}")
    print("-------------------------------------------------------------------------")

def save_solution_image(instance, solver, presence, start, size, HSI, H, K, DIM_STRIP):
    fig, ax = plt.subplots(figsize=(max(12, DIM_STRIP * 0.15), 6))
    cmap   = plt.colormaps.get_cmap('tab10')
    colors = [cmap(h / H) for h in range(H)]
    for (h, s, i) in HSI:
        if solver.value(presence[h, s, i]):
            x = solver.value(start[h, s, i])
            w = solver.value(size[h, s, i])
            ax.add_patch(patches.Rectangle((x, s), w, 1, facecolor=colors[h], edgecolor='none'))
    ax.tick_params(axis='x', rotation=90)
    ax.set_xlim(0, DIM_STRIP)
    ax.set_ylim(0, K)
    ax.set_xticks([x + 0.5 for x in range(DIM_STRIP)])
    ax.set_xticklabels(range(1, DIM_STRIP + 1))
    ax.set_yticks([s + 0.5 for s in range(K)])
    ax.set_yticklabels(range(1, K + 1))
    ax.set_xticks(range(DIM_STRIP + 1), minor=True)
    ax.set_yticks(range(K + 1), minor=True)
    ax.grid(which='minor', color='black', linewidth=0.5)
    ax.set_title(f'Funzione obiettivo: {solver.objective_value}')
    handles = [patches.Patch(color=colors[h], label=f'Specie {h+1}') for h in range(H)]
    ax.legend(handles=handles)
    plt.tight_layout()
    plt.savefig(f"plots/output_{instance}.png", dpi=150, bbox_inches='tight')
    plt.close()

def init_csv(base_dir):
    timestamp  = datetime.now().strftime('%d-%m-%Y_%H-%M-%S')
    csv_path   = base_dir / 'results' / f'{timestamp}.csv'
    csv_file   = open(csv_path, 'w', newline='')
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow([
        'esecuzione_pk', 'istanza_pk', 'istanza',
        'positive', 'negative', 'neutre',
        'specie', 'strip', 'fori', 'cluster',
        'file_id', 'z', 'ctime', 'time', 'status',
        'sinergie', 'conflitti', 'thread', 'sn',
        'time_limit', 'solver', 'maximal'
    ])
    return csv_file, csv_writer, csv_path

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
        0,             # ctime (N/A per OR-Tools)
        wall_time,
        status_code,
        sinergie,
        conflitti,
        num_workers,
        sn_label,
        time_limit,
        'OR-Tools',
        0              # maximal: aggiornato nel post-processing
    ])

def finalize_csv(csv_file, csv_path):
    csv_file.close()
    import pandas as pd
    df = pd.read_csv(csv_path)
    df['maximal'] = df.groupby('istanza_pk')['z'].transform('max')
    df.to_csv(csv_path, index=False)
    print(f"CSV finale salvato in: {csv_path}")