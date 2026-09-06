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

### Primary Paper

> T. Adamo, L. Colizzi, G. Dimauro, E. Guerriero, and D. Pareo, **"Crop planting layout optimization in sustainable agriculture: A constraint programming approach,"** *Computers and Electronics in Agriculture*, vol. 224, p. 109162, 2024. [https://doi.org/10.1016/j.compag.2024.109162](https://doi.org/10.1016/j.compag.2024.109162)

```bibtex
@article{adamo2024cplp,
  author  = {Adamo, Tommaso and Colizzi, Lucio and Dimauro, Giovanni and Guerriero, Emanuela and Pareo, Deborah},
  title   = {Crop planting layout optimization in sustainable agriculture: A constraint programming approach},
  journal = {Computers and Electronics in Agriculture},
  volume  = {224},
  pages   = {109162},
  year    = {2024},
  doi     = {10.1016/j.compag.2024.109162}
}

@article{laborie2018cpo,
  author  = {Laborie, Philippe and Rogerie, J{\'e}r{\^o}me and Shaw, Paul and Vil{\'\i}m, Petr},
  title   = {{IBM ILOG CP} Optimizer for scheduling},
  journal = {Constraints},
  volume  = {23},
  number  = {2},
  pages   = {210--250},
  year    = {2018},
  doi     = {10.1007/s10601-018-9281-x}
}

@incollection{laborie2012interval,
  author    = {Laborie, Philippe and Rogerie, J{\'e}r{\^o}me and Shaw, Paul and Vil{\'\i}m, Petr},
  title     = {Interval-Based Language for Modeling Scheduling Problems: An Extension to Constraint Programming},
  booktitle = {Algebraic Modeling Systems: Modeling and Solving Real World Optimization Problems},
  editor    = {Kallrath, Josef},
  series    = {Applied Optimization},
  volume    = {104},
  pages     = {111--143},
  publisher = {Springer},
  year      = {2012},
  doi       = {10.1007/978-3-642-23592-4_6}
}

@article{focacci2002mp,
  author  = {Focacci, Filippo and Lodi, Andrea and Milano, Michela},
  title   = {Mathematical Programming Techniques in Constraint Programming: A Short Overview},
  journal = {Journal of Heuristics},
  volume  = {8},
  number  = {1},
  pages   = {7--17},
  year    = {2002}
}

@article{bybeefinley2018intercropping,
  author  = {Bybee-Finley, K. A. and Ryan, M. R.},
  title   = {Advancing intercropping research and practices in industrialized agricultural landscapes},
  journal = {Agriculture},
  volume  = {8},
  number  = {6},
  pages   = {80},
  year    = {2018}
}

@misc{ortools_cpsat,
  author       = {Perron, Laurent and Didier, Fr{\'e}d{\'e}ric},
  title        = {{Google OR-Tools CP-SAT} Solver},
  howpublished = {\url{https://developers.google.com/optimization/cp/cp_solver}},
  year         = {2024}
}

@manual{ibm_cpo_manual,
  title        = {{IBM ILOG CP Optimizer} User's Manual},
  organization = {IBM Corporation},
  year         = {2023},
  note         = {Version 22.1.1}
}

@inbook{russell2020csp,
  author    = {Russell, Stuart and Norvig, Peter},
  title     = {Artificial Intelligence: A Modern Approach},
  edition   = {4},
  chapter   = {5},
  note      = {Constraint Satisfaction Problems},
  year      = {2020},
  publisher = {Pearson}
}
```

---
## License

Distributed under the MIT License. See `LICENSE` for more details.
