import configparser
import re
from pathlib import Path

# Parsing del file di configurazione
def parse_config_file(config_file_name):
    config = configparser.ConfigParser()
    config.read(config_file_name)
    allelopathy_threshold = int(config['settings']['allelopathy_threshold'])
    export_results        = config['settings'].getboolean('export_results')
    export_plots          = config['settings'].getboolean('export_plots')
    return allelopathy_threshold, export_results, export_plots

# Parsing della lista delle istanze
def parse_instances_list(instances_file_name):
    with open(instances_file_name) as f:
        return [l.strip() for l in f if l.strip()]

# Parsing dei parametri delle istanze
def parse_dat_file(filepath):
    base_dir = Path(__file__).parent.parent
    with open(base_dir / "instances" / (filepath + ".dat"), 'r') as f:
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

    # file_id: ultima parte del nome istanza (es. I_2_6_33_1_1 -> 1)
    file_id = int(filepath.split('_')[-1])

    return K, M, H, a, o, c_min, c_max, d, positive, negative, neutre, file_id
