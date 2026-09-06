import argparse
import sys
import traceback
from ortools.sat.python import cp_model as cp

from utils.parsing import *
from utils.model import *
from utils.output import *


def solve_single_instance(instance_path, time_limit, constraint_mode, num_workers, allelopathy_threshold, export_plots):

    print(f"\nSolving instance : {instance_path}")
    print(f"Solver parameters  : Time Limit={time_limit}s | Mode={constraint_mode} | Workers={num_workers} | Allelopathy Threshold={allelopathy_threshold}")
    
    # Calls the parse_dat_file function from the parsing module
    try:
        K, M, H, a, o, c_min, c_max, d, positive, negative, neutre, file_id = parse_dat_file(instance_path)
    except Exception as e:
        print(f"\n[Parsing Error] {instance_path} cannot be read.")
        return

    # Calls the CP-SAT solver and the build_and_solve function
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

    # In case of optimal or feasible solution, save the found values
    if status in (cp.OPTIMAL, cp.FEASIBLE):
        z_val      = int(solver.objective_value)
        sinergie   = max(z_val, 0)
        conflitti  = max(-z_val, 0)
        status_str = solver.status_name(status)
    else:
        z_val      = sinergie = conflitti = 0
        status_str = solver.status_name(status)

    wall_time = solver.wall_time

    # Print solution
    print("\n" + "="*45)
    print(" RISULTATI OTTIMIZZAZIONE".center(45))
    print("="*45)
    print(f" Stato Soluzione : {status_str}")
    print(f" Tempo (s)       : {wall_time:.3f}")
    print(f" Obiettivo (Z)   : {z_val}")
    print(f" Sinergie        : {sinergie}")
    print(f" Conflitti       : {conflitti}")
    print("="*45)

    # Export plots
    if export_plots and status in (cp.OPTIMAL, cp.FEASIBLE):
        print("\nGenerazione del layout grafico in corso...")
        try:
            save_solution_image(instance_path, solver, presence, start, size, HSI, H, K, DIM_STRIP)
            print("Plot del layout salvato con successo.")
        except Exception as e:
            print(f"[ERRORE di Plotting]: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Tool per l'ottimizzazione del layout colturale su singola istanza."
    )
    
    # Instance file name
    parser.add_argument(
        "instance",
        type=str,
        help="Percorso del file di istanza"
    )
    
    # Constraint mode
    parser.add_argument(
        "--mode",
        type=str,
        default="hard",
        choices=["hard", "soft"],
        help="Modalità vincoli di vicinanza: 'hard' o 'soft' (default: hard)"
    )
    parser.add_argument(
        "--allelopathy-threshold",
        type=int,
        default=-100,
        help="Soglia al di sotto della quale due specie sono incompatibili (default: -100)"
    )
    
    # CP-SAT solver time limit
    parser.add_argument(
        "--time-limit",
        type=int,
        default=60,
        help="Tempo limite di risoluzione in secondi (default: 60)"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Numero di thread paralleli per CP-SAT (default: 4)"
    )
    
    # Output e visualizzazione
    parser.add_argument(
        "--export-plots",
        action="store_true",
        help="Esporta le soluzioni sotto forma di grafico (default: False)"
    )
    
    args = parser.parse_args()

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
