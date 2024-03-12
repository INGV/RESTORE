import os

filename = 'input_file.txt'

def load_config(file):
    with open(os.path.join(os.getcwd(), filename), 'r') as file:
        lines = file.readlines()

    # Extract configuration parameters
    name_catalog = str(lines[0].split( )[1])    
    delimiter = str(lines[1].split( )[1])
    size = int(lines[2].split( )[1])
    step = int(lines[3].split( )[1])
    st_dev_multiplier = int(lines[4].split( )[1])
    Sigma = float(lines[5].split( )[1])
    sbin = float(lines[6].split( )[1])
    fault_length_multiplier = int(lines[7].split( )[1])
    t_end_quiet = str(lines[8].split( )[1]) 
    if str(lines[9].split( )[1]) == 'None':
        b = None
    else: 
        b = int(lines[9].split( )[1])
    alpha = float(lines[10].split( )[1])
    if str(lines[11].split( )[1]) == 'None':
        mc = None
    else: 
        mc = float(lines[11].split( )[1]) 
    depth_distribution = str(lines[12].split( )[1])  
    p0 = [float(val) for val in lines[13].split( )[1].split(',')]

    # Create configuration dictionary
    config = {
        'name_catalog': name_catalog,
        'delimiter': delimiter,
        'size': size,
        'step': step,
        'st_dev_multiplier': st_dev_multiplier,
        'Sigma': Sigma,
        'sbin': sbin,
        'fault_length_multiplier': fault_length_multiplier,
        't_end_quiet': t_end_quiet,
        'b': b,
        'alpha': alpha,
        'mc': mc,
        'depth_distribution': depth_distribution,
        'p0': p0
    }
    
    return config
