import configparser
import re
import math
from ortools.sat.python import cp_model as cp
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Funzioni di parsing
# Ritorna i parametri dal file di configurazione
def parse_config_file(config_file_name):
    config = configparser.ConfigParser()
    config.read(config_file_name)
    num_workers = int(config['settings']['num_workers'])
    time_limit = int(config['settings']['time_limit'])
    allelopathy_threshold = int(config['settings']['allelopathy_threshold'])
    constraint_mode = config['settings']['constraint_mode']
    return num_workers, time_limit, allelopathy_threshold, constraint_mode

# Ritorna l'elenco delle istanze
def get_instances_list(istancesFileName):
    lines = []
    with open(istancesFileName) as f:
        for l in f:
            lines.append(l.strip())
        return lines

# Ritorna i parametri delle istanze
def parse_dat_file(filepath):
    with open("instances/" + filepath + ".dat", 'r') as file:
        content = file.read()

    def extract_int(name):
        return int(re.search(rf'{name}\s*=\s*(\d+)', content).group(1))

    def extract_int_array(name):
        raw = re.search(rf'{name}\s*=\s*\[(.*?)\]', content, re.DOTALL).group(1)
        return list(map(int, re.findall(r'-?\d+', raw)))

    def extract_2d_array(name):
        raw = re.search(rf'{name}\s*=\s*\[(.*?)\];', content, re.DOTALL).group(1)
        return [list(map(int, re.findall(r'-?\d+', row))) for row in re.findall(r'\[(.*?)\]', raw)]

    K     = extract_int('K')
    M     = extract_int('M')
    H     = extract_int('H')
    a     = extract_2d_array('a')
    o     = extract_int_array('o')
    c_min = extract_int_array('c_min')
    c_max = extract_int_array('c_max')
    d     = extract_int_array('d')

    return K, M, H, a, o, c_min, c_max, d

# Funzioni ausiliarie
# Calcola la distanza minima tra due cluster della stessa specie sulla stessa riga
def get_cluster_distance():
    c_tilde = []
    for h in range(0, H):
        c_tilde.append(min(c_min[k] for k in range(0, H) if k != h))
    return c_tilde

# Calcola la lunghezza di overlap
def overlap_length(start_a, end_a, start_b, end_b, name):
    min_end = model.new_int_var(0, DIM_STRIP, f'min_end_{name}')
    max_start = model.new_int_var(0, DIM_STRIP, f'max_start_{name}')
    diff = model.new_int_var(-DIM_STRIP, DIM_STRIP, f'diff_{name}')
    overlap = model.new_int_var(0, DIM_STRIP, f'overlap_{name}')
    model.add_min_equality(min_end, [end_a, end_b])
    model.add_max_equality(max_start, [start_a, start_b])
    model.add(diff == min_end - max_start)
    model.add_max_equality(overlap, [diff, 0])
    return overlap

# Funzioni di output
# Stampa soluzione sul terminale
def print_solution():
    print(f"\nIstanza: {istance}")
    print(f"Objective: {solver.objective_value}")
    for (h, s, i) in HSI:
        if solver.value(presence[h,s,i]):
            print(f"  z[{h},{s},{i}]: start={solver.value(start[h,s,i])}, "
                f"end={solver.value(end[h,s,i])}, size={solver.value(size[h,s,i])}")

# Stampa soluzione graficamente
def save_solution_image():
    fig, ax = plt.subplots(figsize=(max(12, DIM_STRIP * 0.15), 6))
    fig.canvas.manager.set_window_title('Crop Planting Layout Solver')
    cmap = plt.colormaps.get_cmap('tab10')
    colors = [cmap(h / H) for h in range(H)]
    for (h, s, i) in HSI:
        if solver.value(presence[h,s,i]):
            x = solver.value(start[h,s,i])
            w = solver.value(size[h,s,i])
            ax.add_patch(patches.Rectangle((x, s), w, 1,
                facecolor=colors[h],
                edgecolor='none'))
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
    plt.savefig(f"plots/output_{istance.replace('.dat', '')}.png", dpi=150, bbox_inches='tight')
    plt.close()

num_workers, time_limit, allelopathy_threshold, constraint_mode = parse_config_file('config.ini')
instances = get_instances_list("instances.txt")

for istance in instances:

    K, M, H, a, o, c_min, c_max, d = parse_dat_file(istance)
    DIM_STRIP = math.floor(M / K)

    c_tilde = get_cluster_distance()

    # Numero massimo di cluster di specie h sulla stessa riga
    P = [math.floor((min(d[h]*o[h], DIM_STRIP) + c_tilde[h]) / (c_min[h] + c_tilde[h])) for h in range(H)]
    # Lista di triplette che identificano univocamente il cluster di specie h sulla riga s i-esimo
    HSI = [(h, s, i) for h in range(H) for s in range(K) for i in range(P[h])]

    # Dichiara il modello
    model = cp.CpModel()

    # Variabile decisionale
    z = {}
    start = {}
    end = {}
    size = {}
    presence = {}
    for (h, s, i) in HSI:

        start[h,s,i] = model.new_int_var(0, DIM_STRIP - 1, f"start_h{h}_s{s}_i{i}")
        size[h,s,i] = model.new_int_var(c_min[h], c_max[h], f"size_h{h}_s{s}_i{i}")
        end[h,s,i] = model.new_int_var(1, DIM_STRIP, f"end_h{h}_s{s}_i{i}")
        presence[h,s,i] = model.new_bool_var(f"presence_h{h}_s{s}_i{i}")

        z[h,s,i] = model.new_optional_interval_var(
            start[h,s,i], size[h,s,i], end[h,s,i], presence[h,s,i], f"z_h{h}_s{s}_i{i}"
        )

    # Vincoli
    # La domanda di specie h deve essere soddisfatta
    for h in range(H):
        terms = []
        for s in range(K):
            for i in range(P[h]):
                contrib = model.new_int_var(0, DIM_STRIP, f"contrib_h{h}_s{s}_i{i}")
                model.add_multiplication_equality(contrib, [presence[h,s,i], size[h,s,i]])
                terms.append(contrib)
        model.add(sum(terms) == d[h] * o[h])

    # Dimensioni consentite di ciascun cluster
    for (h, s, i) in HSI:
        allowed = [[v] for v in range(c_min[h], c_max[h] + 1) if v % o[h] == 0]
        model.add_allowed_assignments([size[h,s,i]], allowed)

    # Impedisce l'overlap di cluster sulla stessa strip
    for s in range(K):
        intervals_on_strip = [z[h,s,i] for (h, ss, i) in HSI if ss == s]
        model.add_no_overlap(intervals_on_strip)

    # Verifica che l'hard-constraint mode sia attiva
    if constraint_mode == 'hard':
        # Ricava l'insieme di specie incompatibili
        S = set()
        for h in range(H):
            for k in range(H):
                if a[h][k] <= allelopathy_threshold:
                    S.add((min(h,k), max(h,k)))
        # Aggiunge il vincolo di non adiacenza tra specie incompatibili
        for s in range(K - 1):
            for (h, k) in S:
                intervals_h = [z[h, s, i] for i in range(P[h])]
                intervals_k = [z[k, s+1, j] for j in range(P[k])]
                model.add_no_overlap(intervals_h + intervals_k)
                
                intervals_k2 = [z[k, s, j] for j in range(P[k])]
                intervals_h2 = [z[h, s+1, i] for i in range(P[h])]
                model.add_no_overlap(intervals_k2 + intervals_h2)

    # Symmetry breaking constraint:
    # I cluster devono essere in ordine
    for h in range(H):
        for s in range(K):
            for i in range(1, P[h]):
                #  se il cluster i e' presente, anche i-1 deve esserlo
                model.add(presence[h,s,i-1] == 1).only_enforce_if(presence[h,s,i])
                #  il cluster i inizia dopo la fine del cluster i-1 + c_tilde
                model.add(end[h,s,i-1] + c_tilde[h] <= start[h,s,i]).only_enforce_if(
                    [presence[h,s,i-1], presence[h,s,i]])

    # Funzione obiettivo
    obj_terms = []
    for h in range(H):
        for k in range(H):
            if a[h][k] == 0:
                continue
            for s in range(K - 1):
                for i in range(P[h]):
                    for j in range(P[k]):
                        ol = overlap_length(start[h,s,i], end[h,s,i], start[k,s+1,j], end[k,s+1,j], f"{h}{s}{i}_{k}{s+1}{j}")
                        ol_p1 = model.new_int_var(0, DIM_STRIP, f"ol_p1_{h}{s}{i}_{k}{s+1}{j}")
                        model.add_multiplication_equality(ol_p1, [presence[h,s,i], ol])

                        ol_p2 = model.new_int_var(0, DIM_STRIP, f"ol_p2_{h}{s}{i}_{k}{s+1}{j}")
                        model.add_multiplication_equality(ol_p2, [presence[k,s+1,j], ol_p1])
                        obj_terms.append(a[h][k] * ol_p2)
    model.maximize(sum(obj_terms))

    # Chiama il solver
    solver = cp.CpSolver()
    # Limita il numero di core
    solver.parameters.num_workers = num_workers
    # Imposta un timeout per istanza di 30 secondi
    solver.parameters.max_time_in_seconds = time_limit
    status = solver.solve(model)

    # Se ha trovato una soluzione la stampa su terminale e salva l'immagine nella cartella plots
    if status in (cp.OPTIMAL, cp.FEASIBLE):
        print_solution()
        save_solution_image()
    else:
        print("No solution found:", solver.status_name(status))