import argparse
import sys
import traceback
from ortools.sat.python import cp_model as cp

from utils.parsing import *
from utils.model import *
from utils.output import *


def solve_single_instance(instance_path, time_limit, constraint_mode, num_workers, allelopathy_threshold, export_plots):

    # Debug print
    print(f"\nSolving instance : {instance_path}")
    print(f"Solver parameters  : Time Limit={time_limit}s | Mode={constraint_mode} | Workers={num_workers} | Allelopathy Threshold={allelopathy_threshold}")
    
    # Calls the parse_dat_file function from the parsing module
    try:
        K, M, H, a, o, c_min, c_max, d, positive, negative, neutre, file_id = parse_dat_file(instance_path)
    except Exception as e:
        print(f"\n[Parsing Error] {instance_path} cannot be read.")
        return

    # Calls the CP-SAT solver and the build_and_solve function from the model module
    print("\nStarting CP-SAT Solver...")
    try:
        solver, status, HSI, presence, start, end, size, DIM_STRIP, P = build_and_solve(
            K, M, H, a, o, c_min, c_max, d,
            constraint_mode, allelopathy_threshold,
            num_workers, time_limit
        )
    except Exception as e:
        print(f"\n[Critical Error]: ", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)

    # In case of optimal or feasible solution, saves the found values
    if status in (cp.OPTIMAL, cp.FEASIBLE):
        z_val      = int(solver.objective_value)
        sinergie   = max(z_val, 0)
        conflitti  = max(-z_val, 0)
        status_str = solver.status_name(status)
    else:
        z_val      = sinergie = conflitti = 0
        status_str = solver.status_name(status)

    # Saves the time needed to find the solution
    wall_time = solver.wall_time

    # Prints solution
    print("\n" + "="*45)
    print(" RESULTS ".center(45))
    print("="*45)
    print(f" Solution State : {status_str}")
    print(f" Time (s)       : {wall_time:.3f}")
    print(f" Objective (Z)   : {z_val}")
    print(f" Synergies        : {sinergie}")
    print(f" Conflicts       : {conflitti}")
    print("="*45)

    # Exports plots using the output module
    if export_plots and status in (cp.OPTIMAL, cp.FEASIBLE):
        print("\nGenerating plot...")
        try:
            save_solution_image(instance_path, solver, presence, start, size, HSI, H, K, DIM_STRIP)
            print("Plot saved successfully.")
        except Exception as e:
            print(f"[ERROR. Plot saving failed.]: {e}")


def main():
    # Takes parameters from the CLI
    parser = argparse.ArgumentParser(
        description="Crop Planting Layout solver."
    )
    # Instance file name
    parser.add_argument(
        "instance",
        type=str,
        help="Instance path"
    )
    # Constraint mode
    parser.add_argument(
        "--mode",
        type=str,
        default="hard",
        choices=["hard", "soft"],
        help="Adjacency constraint: 'hard' or 'soft' (default: hard)"
    )
    parser.add_argument(
        "--allelopathy-threshold",
        type=int,
        default=-100,
        help="Below this allelopathy threshold two species are considered incompatibles (default: -100)"
    )
    # CP-SAT solver time limit
    parser.add_argument(
        "--time-limit",
        type=int,
        default=60,
        help="Solver Time Limit (default: 60)"
    )
    # Number of workers
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of parallel threads used by the CP-SAT solver (default: 4)"
    )
    
    # Output e visualizzazione
    parser.add_argument(
        "--export-plots",
        action="store_true",
        help="Exports found solution as a PNG image (default: False)"
    )

    args = parser.parse_args()

    # Calls the solve_single_instance function
    solve_single_instance(
        instance_path=args.instance,
        time_limit=args.time_limit,
        constraint_mode=args.mode,
        num_workers=args.workers,
        allelopathy_threshold=args.allelopathy_threshold,
        export_plots=args.export_plots
    )


if __name__ == '__main__':
    main()
