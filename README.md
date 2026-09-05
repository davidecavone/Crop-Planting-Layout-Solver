# Crop Planting Layout Solver

A Constraint Programming solver in Python designed to address the **Crop Planting Layout Problem (CPLP)**. The project models crops spatial allocation using modern CP-SAT formulations.

---

## 📌 Problem Overview

The Crop Planting Layout Problem involves assigning agricultural crops to clusters of species across discrete planning horizons while respecting spatial constraints:

- **Spatial Allocation:** Preventing overlaps of clusters.
- **Demand Satisfaction:** Meeting the expected units of species.
- **Objective:** Maximizing beneficial interactions between species.

---

## 🛠 Tech Stack & Tools

- **Language:** Python 3.10+
- **Solver Engine:** [Google OR-Tools (CP-SAT)](https://developers.google.com/optimization)
- **Data Handling:** Pandas
- **Visualization:** Matplotlib

---

## 📂 Repository Structure

```text
├── instances/            # Benchmark problem instance files (.dat files)
├── plots/                # Output visual layouts
├── results/              # Output CSV logs
├── utils/                # Utility modules (instances data parsing, output functions)
├── analyze_results.py    # Computational campaign analyzer
├── config.ini            # Global solver configurations and execution parameters
├── cplp_ortools.py       # Core OR-Tools CP-SAT formulation logic
├── instances.txt         # List of instance names for experimental runs
├── requirements.txt      # Project dependencies
├── solver.py             # Computational campaign script
└── README.md
```

---

## 🚀 Getting Started

### 1. Clone the Repository
```bash
git clone [https://github.com/davidecavone/Crop-Planting-Layout-Solver.git](https://github.com/davidecavone/Crop-Planting-Layout-Solver.git)
cd Crop-Planting-Layout-Solver
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

## 💻 Usage

### Running the Solver on a Dataset

Run the main runner script with a specified dataset instance:

```bash
python -m src.solver --input data/instances/instance_01.json --timeout 60
```

### Command-Line Arguments

| Flag | Description | Default |
| :--- | :--- | :--- |
| `--input`, `-i` | Path to the problem instance file (`.json` or `.csv`) | *Required* |
| `--timeout`, `-t` | Maximum search wall-time in seconds | `60` |
| `--threads` | Number of worker threads for parallel search | `8` |
| `--plot` | Generate visual Gantt chart / spatial allocation plot | `False` |

---

## 📊 Visualizing Results

The solver can output a spatial plot with the found solution.

---

## 📈 Benchmarks & Performance

| Instance | Plots / Clusters | Crops / Tasks | Solver Time | Status | Objective Value |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `small_01` | 4 | 12 | 0.42s | Optimal | 1,240 |
| `med_03` | 10 | 45 | 5.81s | Optimal | 4,890 |
| `large_01`| 25 | 150 | 60.00s | Feasible | 14,200 |

---

## 📜 License

Distributed under the MIT License. See `LICENSE` for more details.
