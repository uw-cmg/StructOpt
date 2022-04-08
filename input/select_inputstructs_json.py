import os
import sys
import json
import random

def get_json_data(json_path):
    with open(json_path,'rb') as f:
        params = json.load(f)
        f.close()
    return params

def write_json_data(dict, json_path_out):
    with open(json_path_out,'w') as r:
        json.dump(dict, r, indent=4, sort_keys=False)
    r.close()

def get_input_structures(input_strucs_path, n=10):
    structs = list()
    for dirpath, dirnames, filenames in os.walk(input_strucs_path):
        for file in filenames:
            if 'xyz' in file:
                structs.append(file)
    if n <= len(structs):
        selected_structs = random.sample(structs, n)
    else:
        raise ValueError('the selected number of structures n is larger than the total number of the repository')
    selected_structs = [ input_strucs_path+'/'+i for i in selected_structs]
    return selected_structs

json_path = './TiO2_600.comb.json'
#json_path_out = './TiO2_600_out.comb.json'
input_strucs_path = '/home/mengjun/TiO2_structopt/input_struct/ana600_comb_relaxed'
num_of_individuals = 10
params = get_json_data(json_path)

for json_path_out in ['./1/TiO2_600.comb.json','./2/TiO2_600.comb.json', './3/TiO2_600.comb.json', './4/TiO2_600.comb.json', './5/TiO2_600.comb.json']:
    selected_structs = get_input_structures(input_strucs_path, n=num_of_individuals)
    params['generators']['read_xyz']['number_of_individuals'] = num_of_individuals
    params['generators']['read_xyz']['kwargs'] = selected_structs
    write_json_data(params, json_path_out)
