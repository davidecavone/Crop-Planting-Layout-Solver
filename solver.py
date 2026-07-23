from ortools.sat.python import cp_model as cp
from utils.parsing import *
from utils.model import *
from utils.output import *
from multiprocessing import Pool
import traceback

# ---------------------------------------------------------------------------
# Configurazione campagna
# ---------------------------------------------------------------------------

CAMPAIGN_TIME_LIMITS = [60, 120, 240, 480]

CONFIGURATIONS = [
    ('hard', 1),
    ('soft', 1),
    ('hard', 2),
    ('soft', 2),
]

# Quanti task girare in parallelo.
# Sul Ryzen 9 5950X (32 thread): ogni task usa max 2 thread di CP-SAT,
# quindi con 8 worker usi ~12-16 thread. Puoi alzare a 12 se vuoi.
POOL_SIZE = 8

# ---------------------------------------------------------------------------
# Worker: gestisce una singola (istanza, constraint_mode, num_workers)
# e percorre tutta la sequenza di TL al suo interno
# ---------------------------------------------------------------------------

def run_task(args):
    """
    Esegue la sequenza completa di time limit per una combinazione
    (istanza, constraint_mode, num_workers).
    Ritorna una lista di dict (uno per TL), pronti da scrivere sul CSV.
    """
    istanza_pk, instance, constraint_mode, num_workers, allelopathy_threshold, export_plots = args

    results = []

    # Parsing una sola volta per tutti i TL
    try:
        K, M, H, a, o, c_min, c_max, d, positive, negative, neutre, file_id = parse_dat_file(instance)
        cluster = int(instance.split('_')[4])
    except Exception as e:
        print(f"  [ERRORE parsing] {instance}: {e}")
        return results

    optimal_found  = False
    cached_result  = None

    for time_limit in CAMPAIGN_TIME_LIMITS:

        # --- Replica se ottimale già trovato ---
        if optimal_found and cached_result is not None:
            results.append({**cached_result, 'time_limit': time_limit, 'replicated': True})
            continue

        # --- Risoluzione reale ---
        try:
            solver, status, HSI, presence, start, end, size, DIM_STRIP, P = build_and_solve(
                K, M, H, a, o, c_min, c_max, d,
                constraint_mode, allelopathy_threshold,
                num_workers, time_limit
            )
        except Exception as e:
            print(f"  [ERRORE solver] {instance} | tl={time_limit}s | {constraint_mode} | w={num_workers}: {e}")
            traceback.print_exc()
            continue

        if status in (cp.OPTIMAL, cp.FEASIBLE):
            z_val     = int(solver.objective_value)
            sinergie  = max(z_val, 0)
            conflitti = max(-z_val, 0)
        else:
            z_val     = 0
            sinergie  = 0
            conflitti = 0

        status_code = map_status(status)
        wall_time   = solver.wall_time

        print(f"  {instance} | tl={time_limit}s | w={num_workers} | {constraint_mode} | "
              f"status={solver.status_name(status)} | time={wall_time:.3f}s")

        result = {
            'istanza_pk':     istanza_pk,
            'instance':       instance,
            'positive':       positive,
            'negative':       negative,
            'neutre':         neutre,
            'H':              H,
            'K':              K,
            'DIM_STRIP':      DIM_STRIP,
            'cluster':        cluster,
            'file_id':        file_id,
            'z_val':          z_val,
            'wall_time':      wall_time,
            'status_code':    status_code,
            'sinergie':       sinergie,
            'conflitti':      conflitti,
            'num_workers':    num_workers,
            'constraint_mode': constraint_mode,
            'time_limit':     time_limit,
            'replicated':     False,
        }
        results.append(result)

        # Cache per i TL successivi
        if status_code == 1:
            optimal_found = True
            cached_result = result.copy()

        # Plot (ogni processo ha il suo matplotlib, nessun conflitto)
        if export_plots and status in (cp.OPTIMAL, cp.FEASIBLE):
            try:
                save_solution_image(instance, solver, presence, start, size, HSI, H, K, DIM_STRIP)
            except Exception as e:
                print(f"  [ERRORE plot] {instance}: {e}")

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    base_dir = Path(__file__).parent
    allelopathy_threshold, export_results, export_plots = parse_config_file(base_dir / 'config.ini')
    instances = parse_instances_list(base_dir / "instances.txt")

    # Costruisce la lista di task: una per ogni (istanza x configurazione)
    tasks = [
        (istanza_pk, instance, constraint_mode, num_workers, allelopathy_threshold, export_plots)
        for istanza_pk, instance in enumerate(instances, start=1)
        for constraint_mode, num_workers in CONFIGURATIONS
    ]

    print(f"Istanze: {len(instances)} | Configurazioni: {len(CONFIGURATIONS)} | "
          f"Task totali: {len(tasks)} | Pool size: {POOL_SIZE}")

    if export_results:
        csv_file, csv_writer, csv_path = init_csv(base_dir)

    esecuzione_pk     = 0
    task_completati   = 0

    # imap_unordered: i risultati arrivano man mano che i worker finiscono
    with Pool(processes=POOL_SIZE) as pool:
        for task_results in pool.imap_unordered(run_task, tasks, chunksize=1):
            task_completati += 1
            print(f"\n[Task {task_completati}/{len(tasks)} completato]")

            for result in task_results:
                esecuzione_pk += 1

                if export_results:
                    write_row(
                        csv_writer,
                        esecuzione_pk,
                        result['istanza_pk'],
                        result['instance'],
                        result['positive'],
                        result['negative'],
                        result['neutre'],
                        result['H'],
                        result['K'],
                        result['DIM_STRIP'],
                        result['cluster'],
                        result['file_id'],
                        result['z_val'],
                        result['wall_time'],
                        result['status_code'],
                        result['sinergie'],
                        result['conflitti'],
                        result['num_workers'],
                        result['constraint_mode'],
                        result['time_limit'],
                    )
                    csv_file.flush()

    if export_results:
        finalize_csv(csv_file, csv_path)


if __name__ == '__main__':
    main()
