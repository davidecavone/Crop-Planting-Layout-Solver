from ortools.sat.python import cp_model as cp
import math

def get_cluster_distance(H, c_min):
    return [min(c_min[k] for k in range(H) if k != h) for h in range(H)]

def overlap_length(model, start_a, end_a, start_b, end_b, name, DIM_STRIP):
    min_end   = model.new_int_var(0, DIM_STRIP, f'min_end_{name}')
    max_start = model.new_int_var(0, DIM_STRIP, f'max_start_{name}')
    diff      = model.new_int_var(-DIM_STRIP, DIM_STRIP, f'diff_{name}')
    overlap   = model.new_int_var(0, DIM_STRIP, f'overlap_{name}')
    model.add_min_equality(min_end, [end_a, end_b])
    model.add_max_equality(max_start, [start_a, start_b])
    model.add(diff == min_end - max_start)
    model.add_max_equality(overlap, [diff, 0])
    return overlap

def map_status(cp_status):
    """Mappa lo status CP-SAT al codice numerico della campagna (-1/0/1)."""
    if cp_status == cp.OPTIMAL:
        return 1
    elif cp_status == cp.FEASIBLE:
        return 0
    else:
        return -1

# Callback prima soluzione
class FirstSolutionCallback(cp.CpSolverSolutionCallback):
    def __init__(self):
        cp.CpSolverSolutionCallback.__init__(self)
        self._first_solution_time = None

    def on_solution_callback(self):
        if self._first_solution_time is None:
            self._first_solution_time = self.wall_time

    def first_solution_time(self):
        return self._first_solution_time

# Creazione e risoluzione del modello
def build_and_solve(K, M, H, a, o, c_min, c_max, d,
                    constraint_mode, allelopathy_threshold,
                    num_workers, time_limit):

    DIM_STRIP = math.floor(M / K)
    c_tilde   = get_cluster_distance(H, c_min)
    P         = [math.floor((min(d[h]*o[h], DIM_STRIP) + c_tilde[h]) / (c_min[h] + c_tilde[h])) for h in range(H)]
    HSI       = [(h, s, i) for h in range(H) for s in range(K) for i in range(P[h])]

    model = cp.CpModel()

    z        = {}
    start    = {}
    end      = {}
    size     = {}
    presence = {}

    for (h, s, i) in HSI:
        start[h,s,i]    = model.new_int_var(0, DIM_STRIP - 1, f"start_h{h}_s{s}_i{i}")
        size[h,s,i]     = model.new_int_var(c_min[h], c_max[h], f"size_h{h}_s{s}_i{i}")
        end[h,s,i]      = model.new_int_var(1, DIM_STRIP, f"end_h{h}_s{s}_i{i}")
        presence[h,s,i] = model.new_bool_var(f"presence_h{h}_s{s}_i{i}")
        z[h,s,i]        = model.new_optional_interval_var(
            start[h,s,i], size[h,s,i], end[h,s,i], presence[h,s,i], f"z_h{h}_s{s}_i{i}"
        )

    # Domanda soddisfatta
    for h in range(H):
        terms = []
        for s in range(K):
            for i in range(P[h]):
                contrib = model.new_int_var(0, DIM_STRIP, f"contrib_h{h}_s{s}_i{i}")
                model.add_multiplication_equality(contrib, [presence[h,s,i], size[h,s,i]])
                terms.append(contrib)
        model.add(sum(terms) == d[h] * o[h])

    # Dimensioni consentite
    for (h, s, i) in HSI:
        allowed = [[v] for v in range(c_min[h], c_max[h] + 1) if v % o[h] == 0]
        model.add_allowed_assignments([size[h,s,i]], allowed)

    # No overlap sulla stessa strip
    for s in range(K):
        intervals_on_strip = [z[h,s,i] for (h, ss, i) in HSI if ss == s]
        model.add_no_overlap(intervals_on_strip)

    # Hard constraints: non adiacenza specie incompatibili
    if constraint_mode == 'hard':
        S = set()
        for h in range(H):
            for k in range(H):
                if a[h][k] <= allelopathy_threshold:
                    S.add((min(h,k), max(h,k)))
        for s in range(K - 1):
            for (h, k) in S:
                model.add_no_overlap([z[h, s, i] for i in range(P[h])] + [z[k, s+1, j] for j in range(P[k])])
                model.add_no_overlap([z[k, s, j] for j in range(P[k])] + [z[h, s+1, i] for i in range(P[h])])

    # Symmetry breaking
    for h in range(H):
        for s in range(K):
            for i in range(1, P[h]):
                model.add(presence[h,s,i-1] == 1).only_enforce_if(presence[h,s,i])
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
                        ol    = overlap_length(model, start[h,s,i], end[h,s,i],
                                               start[k,s+1,j], end[k,s+1,j],
                                               f"{h}{s}{i}_{k}{s+1}{j}", DIM_STRIP)
                        ol_p1 = model.new_int_var(0, DIM_STRIP, f"ol_p1_{h}{s}{i}_{k}{s+1}{j}")
                        model.add_multiplication_equality(ol_p1, [presence[h,s,i], ol])
                        ol_p2 = model.new_int_var(0, DIM_STRIP, f"ol_p2_{h}{s}{i}_{k}{s+1}{j}")
                        model.add_multiplication_equality(ol_p2, [presence[k,s+1,j], ol_p1])
                        obj_terms.append(a[h][k] * ol_p2)
    model.maximize(sum(obj_terms))

    solver   = cp.CpSolver()
    solver.parameters.num_workers          = num_workers
    solver.parameters.max_time_in_seconds  = time_limit
    callback = FirstSolutionCallback()
    status   = solver.solve(model, callback)

    return solver, status, HSI, presence, start, end, size, DIM_STRIP, P