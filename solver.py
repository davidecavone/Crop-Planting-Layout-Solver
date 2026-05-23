from ortools.sat.python import cp_model as cp
from utils.parsing import *
from utils.model import *
from utils.output import *

# Tempi limite della campagna
CAMPAIGN_TIME_LIMITS = [60, 120, 240, 480]

# Chiama le funzioni di parsing e legge il file di configurazione e la lista delle istanze
base_dir = Path(__file__).parent
allelopathy_threshold, export_results, export_plots = parse_config_file(base_dir/'config.ini')
instances = parse_instances_list(base_dir/"instances.txt")[:36]

# Apre il CSV una volta sola
if export_results:
    csv_file, csv_writer, csv_path = init_csv(base_dir)

esecuzione_pk = 0

CONFIGURATIONS = [
    ('hard', 1),
    ('soft', 1),
    ('hard', 2),
    ('soft', 2),
]

# Dizionario per tenere traccia dei risultati ottimali già trovati
# chiave: (istanza_pk, constraint_mode, num_workers)
optimal_cache = {}

for time_limit in CAMPAIGN_TIME_LIMITS:
    for constraint_mode, num_workers in CONFIGURATIONS:
        for istanza_pk, instance in enumerate(instances, start=1):

            esecuzione_pk += 1
            print(f"\n[{esecuzione_pk}] {instance} | tl={time_limit}s | workers={num_workers} | mode={constraint_mode}")

            cache_key = (istanza_pk, constraint_mode, num_workers)

            # Se era già ottimale con un time limit precedente, replica
            if cache_key in optimal_cache:
                print(f"  -> Replica (ottimale precedente)")
                if export_results:
                    cached = optimal_cache[cache_key]
                    write_row(csv_writer, esecuzione_pk, istanza_pk, instance,
                            cached['positive'], cached['negative'], cached['neutre'],
                            cached['H'], cached['K'], cached['DIM_STRIP'], cached['cluster'], cached['file_id'],
                            cached['z_val'], cached['wall_time'],
                            cached['status_code'], cached['sinergie'], cached['conflitti'],
                            num_workers, constraint_mode, time_limit)
                    csv_file.flush()
                continue

            # Esecuzione reale
            K, M, H, a, o, c_min, c_max, d, positive, negative, neutre, file_id = parse_dat_file(instance)
            cluster = int(instance.split('_')[4])

            solver, status, HSI, presence, start, end, size, DIM_STRIP, P = build_and_solve(
                K, M, H, a, o, c_min, c_max, d,
                constraint_mode, allelopathy_threshold,
                num_workers, time_limit
            )

            print(f"  Status: {solver.status_name(status)} | Time: {solver.wall_time:.3f}s")

            if status in (cp.OPTIMAL, cp.FEASIBLE):
                z_val     = int(solver.objective_value)
                sinergie  = max(z_val, 0)
                conflitti = max(-z_val, 0)
            else:
                z_val     = 0
                sinergie  = 0
                conflitti = 0

            status_code = map_status(status)
            wall_time   = round(solver.wall_time, 3)

            # Se ottimale, salva in cache per i time limit successivi
            if status_code == 1:
                optimal_cache[cache_key] = {
                    'positive':    positive,
                    'negative':    negative,
                    'neutre':      neutre,
                    'cluster':     cluster,
                    'H':           H,
                    'K':           K,
                    'DIM_STRIP':   DIM_STRIP,
                    'file_id':     file_id,
                    'z_val':       z_val,
                    'wall_time':   wall_time,
                    'status_code': status_code,
                    'sinergie':    sinergie,
                    'conflitti':   conflitti,
                }

            if export_results:

                write_row(csv_writer, esecuzione_pk, istanza_pk, instance,
                        positive, negative, neutre, H, K, DIM_STRIP, cluster, file_id,
                        z_val, wall_time, status_code, sinergie, conflitti,
                        num_workers, constraint_mode, time_limit)
                csv_file.flush()

            if status in (cp.OPTIMAL, cp.FEASIBLE):
                print_solution(instance, solver, presence, start, end, size, HSI)
                if export_plots:
                    save_solution_image(instance, solver, presence, start, size, HSI, H, K, DIM_STRIP)

# Post processing
if export_results:
    finalize_csv(csv_file, csv_path)