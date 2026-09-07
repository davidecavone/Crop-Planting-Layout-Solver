from ortools.sat.python import cp_model as cp

# Map solver status to a number for campaign analysis
def map_status(cp_status):
    if cp_status == cp.OPTIMAL:
        return 1
    elif cp_status == cp.FEASIBLE:
        return 0
    else:
        return -1

# Callback for the first found solution
class FirstSolutionCallback(cp.CpSolverSolutionCallback):
    def __init__(self):
        cp.CpSolverSolutionCallback.__init__(self)
        self._first_solution_time = None

    def on_solution_callback(self):
        if self._first_solution_time is None:
            self._first_solution_time = self.wall_time

    def first_solution_time(self):
        return self._first_solution_time
