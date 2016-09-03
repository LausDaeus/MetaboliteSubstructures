import src.utils.verbose as verbose
import os


def read_smiles(file_names):
    """
    Takes in a txt file, reads each line and stores the result in a list.
    Returns another list when using_lists
    :param file_name: path to smiles file
    :param using_lists: flag to create a list that keeps each list the smiles string is in.
    :return: a tuple of (smiles strings, list_memberships)
    """

    list_membership = []
    list_counter = 0
    smiles_set = []


    for file in file_names:
        with open(file,'r') as input:
            for line in input:
                if not line.strip():
                    list_counter += 1
                else:
                    smiles_set.append(line)
                    list_membership.append(list_counter)

        list_counter += 1

    verbose.info("data_input",list_counter, len(smiles_set))

    return (smiles_set, list_membership)

def read_spectra(file_name):
    with open(file_name,'r') as input:
        energy_type = "unknown energy"
        spectra = {}
        for line in input:
            if not any(char.isdigit() for char in line):
                energy_type = line.rstrip()
                spectra[energy_type] = []
            elif not line in ['\n','\r\n']:
                mass, intensity = line.split(" ")
                spectra[energy_type].append([float(mass),float(intensity)])
            else:
                continue

        return spectra

def read_file(file_name):
    with open(file_name, mode='r') as reader:
        return [line.rstrip() for line in reader]

def write_representative_structure(structure_list, file_name):
    with open(file_name, mode='wb') as writer:
        for smiles in structure_list:
            writer.write(smiles+"\n")

def write_to_file(text, file_name):
    with open(file_name, mode='wb') as writer:
        writer.write(text+"\n")

def check_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)