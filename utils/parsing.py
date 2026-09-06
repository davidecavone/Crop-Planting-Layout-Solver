import re
from pathlib import Path

# ---------------------------------------------------------------------------
# Instances batch parser (used for benchmarks)
# ---------------------------------------------------------------------------
def parse_instances_list(instances_file_name):
    with open(instances_file_name) as f:
        return [l.strip() for l in f if l.strip()]

# ---------------------------------------------------------------------------
# Instance parser (.dat)
# ---------------------------------------------------------------------------
def parse_dat_file(filepath):
    path_obj = Path(filepath)

    # Aggiunge .dat se omesso
    if path_obj.suffix != '.dat':
        path_obj = path_obj.with_suffix('.dat')

    # Se non esiste al path dato, fallback nella cartella instances/
    if not path_obj.exists():
        base_dir = Path(__file__).parent.parent
        fallback_path = base_dir / "instances" / path_obj.name
        if fallback_path.exists():
            path_obj = fallback_path
        else:
            raise FileNotFoundError(f"File non trovato né in '{filepath}' né in '{fallback_path}'")

    with open(path_obj, 'r') as f:
        content = f.read()

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

    positive = sum(1 for h in range(H) for k in range(h+1, H) if a[h][k] > 0)
    negative = sum(1 for h in range(H) for k in range(h+1, H) if a[h][k] < 0)
    neutre   = sum(1 for h in range(H) for k in range(h+1, H) if a[h][k] == 0)

    try:
        file_id = int(path_obj.stem.split('_')[-1])
    except (ValueError, IndexError):
        file_id = 0

    return K, M, H, a, o, c_min, c_max, d, positive, negative, neutre, file_id
