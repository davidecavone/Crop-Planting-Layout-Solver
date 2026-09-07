import argparse
from multiprocessing import Pool
from pathlib import Path
import sys
import traceback

# Add project root (Crop-Planting-Layout-Solver) to sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# Sibling module (same folder: computational_campaign/)
from csv_log import *

# Project-level modules (from utils/)
from ortools.sat.python import cp_model as cp
from utils.model import *
from utils.output import *
from utils.parsing import *

# Computational campaign configuration
from config import CAMPAIGN_TIME_LIMITS, CONFIGURATIONS, POOL_SIZE

# Solve the instance for every combination of constraint mode, number of workers, and time limit
def run_task(args):

    istanza_pk, instance, constraint_mode, num_workers, allelopathy_threshold, export_plots = args

    results = []

    # Parsing una sola volta per tutti i TL
    try:
        K, M, H, a, o, c_min, c_max, d, positive, negative, neutre, file_id = parse_dat_file(instance)
        cluster = int(instance.split('_')[4])
    except Exception as e:
        print(f"  [ERROR parsing] {instance}: {e}")
        return results

    optimal_found  = False
    cached_result  = None

    for time_limit in CAMPAIGN_TIME_LIMITS:

        # If optimal solution is found, copy it for higher time limits
        if optimal_found and cached_result is not None:
            results.append({**cached_result, 'time_limit': time_limit, 'replicated': True})
            continue

        # Solver logic
        try:
            # Build the model
            solver, status, HSI, presence, start, end, size, DIM_STRIP, P = build_and_solve(
                K, M, H, a, o, c_min, c_max, d,
                constraint_mode, allelopathy_threshold,
                num_workers, time_limit
            )
        except Exception as e:
            print(f"  [ERROR solver] {instance} | tl={time_limit}s | {constraint_mode} | w={num_workers}: {e}")
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
        
        # Debug print
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

        # Cache for higher time limits
        if status_code == 1:
            optimal_found = True
            cached_result = result.copy()

        # Export plots
        if export_plots and status in (cp.OPTIMAL, cp.FEASIBLE):
            try:
                save_solution_image(instance, solver, presence, start, size, HSI, H, K, DIM_STRIP)
            except Exception as e:
                print(f"  [ERROR plot] {instance}: {e}")

    return results

def main():
    base_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Run computational campaign.")
    
    # CLI arguments
    # Instance list file path
    parser.add_argument(
        "instances_file",
        type=Path,
        help="Path to the .txt file containing the list of instance file names"
    )
    # Allelopathy threshold
    parser.add_argument(
        "--allelopathy-threshold",
        type=int,
        default=-100,
        help="Below this allelopathy threshold two species are considered incompatibles (default: -100)"
    )
    # Export results (boolean flag defaulting to True)
    parser.add_argument(
        "--export-results",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Export campaign results in a CSV file (default: True)"
    )
    # Export plots (boolean flag defaulting to False)
    parser.add_argument(
        "--export-plots",
        action="store_true",
        default=False,
        help="Exports found solutions as PNG images (default: False)"
    )

    args = parser.parse_args()

    instances = parse_instances_list(args.instances_file)
    allelopathy_threshold = args.allelopathy_threshold
    export_results = args.export_results
    export_plots = args.export_plots

    # Makes a task list containing every configuration possible for eeach instance
    tasks = [
        (istanza_pk, instance, constraint_mode, num_workers, allelopathy_threshold, export_plots)
        
        # pairs each element with a 1-based counter, example: (11, I_2_6_33_1_1.dat)
        for istanza_pk, instance in enumerate(instances, start=1)
        for constraint_mode, num_workers in CONFIGURATIONS
    ]

    print(f"Number of instances: {len(instances)} | Number of testing configurations types: {len(CONFIGURATIONS)} | "
          f"Total tasks: {len(tasks)} | Pool size: {POOL_SIZE}")

    if export_results:
        csv_file, csv_writer, csv_path = init_csv(base_dir)

    esecuzione_pk     = 0
    task_completati   = 0

    with Pool(processes=POOL_SIZE) as pool:
        for task_results in pool.imap_unordered(run_task, tasks, chunksize=1):
            task_completati += 1
            print(f"\n[{task_completati}/{len(tasks)} tasks completed.]")

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
