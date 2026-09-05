import argparse
from pathlib import Path
import traceback
from ortools.sat.python import cp_model as cp

from utils.parsing import *
from utils.model import *
from utils.output import *

def solve_single_instance(instance_path, time_limit, constraint_mode, num_workers, allelopathy_threshold, export_plots):
    """
    Risolve una singola istanza del problema e stampa i risultati a schermo,
    adatta per demo lato client o esecuzione stand-alone.
    """
    print(f"\nRisoluzione istanza : {instance_path}")
    print(f"Parametri Solver  : TL={time_limit}s | Mode={constraint_mode} | Workers={num_workers}")
    
    # --- 1. Parsing dell'istanza ---
    try:
        K, M, H, a, o, c_min, c_max, d, positive, negative, neutre, file_id = parse_dat_file(instance_path)
    except Exception as e:
        print(f"\n[ERRORE di Parsing] Impossibile leggere {instance_path}: {e}")
        return

    # --- 2. Risoluzione con CP-SAT ---
    print("\nAvvio del solver CP-SAT...")
    try:
        solver, status, HSI, presence, start, end, size, DIM_STRIP, P = build_and_solve(
            K, M, H, a, o, c_min, c_max, d,
            constraint_mode, allelopathy_threshold,
            num_workers, time_limit
        )
    except Exception as e:
        print(f"\n[ERRORE del Solver]: {e}")
        traceback.print_exc()
        return

    # --- 3. Estrazione Metriche ---
    if status in (cp.OPTIMAL, cp.FEASIBLE):
        z_val      = int(solver.objective_value)
        sinergie   = max(z_val, 0)
        conflitti  = max(-z_val, 0)
        status_str = solver.status_name(status)
    else:
        z_val      = sinergie = conflitti = 0
        status_str = solver.status_name(status)

    wall_time = solver.wall_time

    # --- 4. Report lato Client / Stakeholder ---
    print("\n" + "="*45)
    print(" RISULTATI OTTIMIZZAZIONE".center(45))
    print("="*45)
    print(f" Stato Soluzione : {status_str}")
    print(f" Tempo (s)       : {wall_time:.3f}")
    print(f" Obiettivo (Z)   : {z_val}")
    print(f" Sinergie        : {sinergie}")
    print(f" Conflitti       : {conflitti}")
    print("="*45)

    # --- 5. Export Grafico della Soluzione ---
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
    
    # Argomento obbligatorio
    parser.add_argument(
        "instance", 
        type=str, 
        help="Percorso del file .dat dell'istanza (es. instances/istanza_01.dat)"
    )
    
    # Argomenti opzionali di configurazione del solver
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
    parser.add_argument(
        "--mode", 
        type=str, 
        default="hard", 
        choices=["hard", "soft"], 
        help="Modalità vincoli di vicinanza: 'hard' o 'soft' (default: hard)"
    )
    parser.add_argument(
        "--allelopathy_threshold",
        type=int,
        default="-100",
        help="Soglia al di sotto della quale due specie sono considerate incompatibili"
    )
    parser.add_argument(
    "--export-plots",
    action="store_true",
    help="Esporta le soluzioni sotto forma di grafico."
    )
    
    args = parser.parse_args()

    # Lettura delle impostazioni globali (soglia agronomica e abilitazione plot)
    base_dir = Path(__file__).parent
    config_path = base_dir / 'config.ini'
    
    try:
        allelopathy_threshold, _, export_plots = parse_config_file(config_path)
    except Exception as e:
        print(f"Attenzione: Impossibile leggere {config_path.name} ({e}). Applicati valori di default.")
        allelopathy_threshold = 0
        export_plots = True

    solve_single_instance(
        instance_path=args.instance,
        time_limit=args.time_limit,
        constraint_mode=args.mode,
        num_workers=args.workers,
        allelopathy_threshold=allelopathy_threshold,
        export_plots=export_plots
    )


if __name__ == '__main__':
    main()
