from ortools.sat.python import cp_model as cp
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def map_status(status):
    """
    Maps OR-Tools CP-SAT status codes to standardized integer IDs:
      1: OPTIMAL
      2: FEASIBLE
      3: INFEASIBLE
      4: MODEL_INVALID
      0: UNKNOWN / OTHER
    """
    if status == cp.OPTIMAL:
        return 1
    elif status == cp.FEASIBLE:
        return 2
    elif status == cp.INFEASIBLE:
        return 3
    elif status == cp.MODEL_INVALID:
        return 4
    return 0

# Print solution to terminal
def print_to_terminal(status_str, wall_time, z_val, synergies, conflicts):
    print("\n" + "="*45)
    print(" RESULTS ".center(45))
    print("="*45)
    print(f" Solution State : {status_str}")
    print(f" Time (s)       : {wall_time:.3f}")
    print(f" Objective (Z)   : {z_val}")
    print(f" Synergies        : {synergies}")
    print(f" Conflicts       : {conflicts}")
    print("="*45)


# Plots the solution found and saves it as a PNG image
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
    ax.set_title(f'Total score: {solver.objective_value}')
    handles = [patches.Patch(color=colors[h], label=f'Specie {h+1}') for h in range(H)]
    ax.legend(handles=handles)
    plt.tight_layout()
    plt.savefig(f"plots/output_{instance}.png", dpi=150, bbox_inches='tight')
    plt.close()
