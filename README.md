# Crop Planting Layout Solver

A Constraint Programming solver in Python designed to address the **Crop Planting Layout Problem (CPLP)**. The project models agricultural plot scheduling, spatial allocation, and rotation constraints using modern CP-SAT formulations.

---

## 📌 Problem Overview

The Crop Planting Layout Problem involves assigning agricultural crops to land plots (or clusters) across discrete planning horizons while respecting operational, biological, and spatial constraints:

- **Spatial Allocation & Overlaps:** Preventing resource conflicts on identical plots/clusters over time (`noOverlap`).
- **Agronomic Constraints:** Respecting crop growth cycles, fallow intervals, and compatibility sequences.
- **Demand Satisfaction:** Meeting periodic harvest yields and volume targets.
- **Objective:** Minimizing operational costs, idle land periods, or schedule delays (makespan/penalty objectives).

---

## 🛠 Tech Stack & Tools

- **Language:** Python 3.10+
- **Solver Engine:** [Google OR-Tools (CP-SAT)](https://developers.google.com/optimization) *(or IBM ILOG CP Optimizer / Docplex)*
- **Data Handling:** NumPy, Pandas
- **Visualization:** Matplotlib / Seaborn (Gantt charts & schedule layouts)

---

## 📂 Repository Structure

```text
├── data/
│   ├── instances/          # Benchmark instances / raw input datasets (JSON/CSV)
│   └── solutions/          # Output assignments and scheduling logs
├── src/
│   ├── model.py            # Core CP formulation (variables, intervals, constraints)
│   ├── parser.py           # Instance reader and parameter mapper
│   ├── solver.py           # Engine invocation and parameter tuning
│   └── visualize.py        # Schedule/Gantt chart generation
├── notebooks/              # Exploratory data analysis and model validation
├── tests/                  # Unit tests for constraints and edge cases
├── requirements.txt        # Project dependencies
└── README.md
