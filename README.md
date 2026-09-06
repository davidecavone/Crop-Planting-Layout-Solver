# Crop Planting Layout Solver

A Constraint Programming solver in Python designed to address the **Crop Planting Layout Problem (CPLP)**. The project models crops spatial allocation using modern CP-SAT formulations.

---

## Problem Overview

The Crop Planting Layout Problem involves assigning agricultural crops to clusters of species across discrete planning horizons while respecting spatial constraints:

- **Spatial Allocation:** Preventing overlaps of clusters.
- **Demand Satisfaction:** Meeting the expected units of species.
- **Objective:** Maximizing beneficial interactions between species.

---

## Tech Stack & Tools

- **Language:** Python 3.10+
- **Solver Engine:** [Google OR-Tools (CP-SAT)](https://developers.google.com/optimization)
- **Data Handling:** Pandas
- **Visualization:** Matplotlib

---

## Repository Structure

```text
.
├── utils/
├── README.md
├── analyze_benchmarks.py
├── main.py
├── requirements.txt
└── run_benchmarks.py
```

---

## Getting Started

### 1. Clone the Repository
```bash
git clone https://github.com/davidecavone/Crop-Planting-Layout-Solver.git
```

### 2. Environment Setup
Create a clean virtual environment and install the dependencies:

```bash
python -m venv venv
# Linux / macOS:
source venv/bin/activate
# Windows:
# .\venv\Scripts\activate

pip install -r requirements.txt
```

---

## Usage

```bash
python main.py <instance> [-h] [--mode {hard,soft}] [--allelopathy-threshold INT]
                          [--time-limit INT] [--workers INT] [--export-plots]
```

---
## References

If you use this work or reference the underlying model and tools, please cite the corresponding sources:

> Tommaso Adamo, L. Colizzi, G. Dimauro, Emanuela Guerriero, and D. Pareo, **"Crop planting layout optimization in sustainable agriculture: A constraint programming approach,"** *Computers and Electronics in Agriculture*, vol. 224, p. 109162, 2024. [https://doi.org/10.1016/j.compag.2024.109162](https://doi.org/10.1016/j.compag.2024.109162)

> Philippe Laborie, Jérôme Rogerie, Paul Shaw, and Petr Vilím, **"IBM ILOG CP Optimizer for scheduling,"** *Constraints*, vol. 23, no. 2, pp. 210–250, 2018. [https://doi.org/10.1007/s10601-018-9281-x](https://doi.org/10.1007/s10601-018-9281-x)

> Philippe Laborie, Jérôme Rogerie, Paul Shaw, and Petr Vilím, **"Interval-Based Language for Modeling Scheduling Problems: An Extension to Constraint Programming,"** in *Algebraic Modeling Systems: Modeling and Solving Real World Optimization Problems*, Applied Optimization, vol. 104, J. Kallrath, Ed., Springer, pp. 111–143, 2012. [https://doi.org/10.1007/978-3-642-23592-4_6](https://doi.org/10.1007/978-3-642-23592-4_6)

> Filippo Focacci, Andrea Lodi, and Michela Milano, **"Mathematical Programming Techniques in Constraint Programming: A Short Overview,"** *Journal of Heuristics*, vol. 8, no. 1, pp. 7–17, 2002.

> K. A. Bybee-Finley and M. R. Ryan, **"Advancing intercropping research and practices in industrialized agricultural landscapes,"** *Agriculture*, vol. 8, no. 6, p. 80, 2018.

> Laurent Perron and Frédéric Didier, **"Google OR-Tools CP-SAT Solver,"** 2024. [https://developers.google.com/optimization/cp/cp_solver](https://developers.google.com/optimization/cp/cp_solver)

> IBM Corporation, **"IBM ILOG CP Optimizer User's Manual,"** Version 22.1.1, 2023.

> Stuart Russell and Peter Norvig, **"Constraint Satisfaction Problems,"** in *Artificial Intelligence: A Modern Approach*, 4th ed., ch. 5, Pearson, 2020.

---
## License

Distributed under the MIT License. See `LICENSE` for more details.
