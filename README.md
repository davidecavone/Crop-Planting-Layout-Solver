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

## License

Distributed under the MIT License. See `LICENSE` for more details.
