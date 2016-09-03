from src.cs_algorithms.characteristic_substructure import CSAlgorithm
import src.utils.verbose as verbose
import src.utils.fileIO as io
import os
from src.parsing.smiles_molecule_parser import Parser

def find_most_representative_molecule(args, output_dir):

    if len(args.smiles_files) % 2 != 0:
        print "Number of input smiles files NEED to be an even number"
        print "Given:",args.smiles_files
        print "Number:",len(args.smiles_files)
        return

    file_path, extension = os.path.splitext(args.molecule_check_list)
    output_directory = os.path.dirname(output_dir)+"/"+os.path.basename(file_path)

    verbose.info("total_time_start")
    set_length = len(args.smiles_files)/2

    cs_A = setup_CS("CS A", args.smiles_files[:set_length], args)
    cs_B = setup_CS("CS B", args.smiles_files[set_length:], args)
    candidate_smiles = io.read_smiles([args.molecule_check_list])[0]

    best_score = -1
    best_candidate = None
    best_smiles = None
    best_score_counter = 0


    for smiles in candidate_smiles:
        smiles = smiles.rstrip()
        candidate = Parser().parse_smiles(smiles)
        candidate.set_list_membership(0)
        score_A = assess_representativeness(candidate, cs_A)
        score_B = assess_representativeness(candidate, cs_B)

        score = calculate_score(score_A, score_B)
        print smiles, score
        if not best_candidate or score > best_score:
            best_candidate = candidate
            best_score = score
            best_smiles = smiles

    if best_candidate:
        smiles_A = get_CS(cs_A, best_candidate)
        smiles_B = get_CS(cs_B, best_candidate)
        results = "Best candidate: "+best_smiles
        results += "\nScore: "+str(best_score)
        results += "\n"+cs_A.name+":"+smiles_A
        results += "\n"+cs_B.name+":"+smiles_B

        print "\n"+results
        io.check_dir(output_directory)
        io.write_to_file(results, output_directory + '/representative_molecule.txt')

        if args.images_of_molecules:
            import src.drawing.draw_molecule as draw
            draw.draw_smiles(best_smiles, output_directory+"/representative_mol.png")
            draw.draw_smiles(smiles_A, output_directory+"/cs_A.png")
            draw.draw_smiles(smiles_B, output_directory+"/cs_B.png")
    else:
        print "No candidate found."

    verbose.info("total_time_finish")

def calculate_score(score_A, score_B):
    mean = (score_A+score_B)/2.0
    variance = ((mean-score_A)**2+(mean-score_B)**2)/2.0
    st_dev = variance**0.5
    return mean - st_dev

def get_CS(cs, best_candidate):
    cs.molecules.append(best_candidate)
    paths = cs.find_graphs_paths()
    cs.find_characteristic_substructure(paths)
    smiles = cs.characteristic_substructure.to_smiles()
    cs.clear_all_structures()
    return smiles

def assess_representativeness(candidate, cs):

    cs.molecules.append(candidate)

    verbose.info("start_paths")
    paths= cs.find_graphs_paths()
    verbose.info("finish_paths")

    rep_structs = cs.find_all_representative_structures(paths, name=cs.name)
    score = score_structures(rep_structs)
    cs.clear_all_structures()
    cs.molecules.pop()
    return score

def score_structures(rep_structs_A):

    return sum([len(x.vertices())*rep_structs_A[x] for x in rep_structs_A])

def setup_CS(name, smiles_files, args):
    cs = CSAlgorithm(**vars(args))
    cs.name = name
    smiles_set, list_memberships = io.read_smiles(smiles_files)
    cs.list_number = len(set(list_memberships))
    if cs.list_number > 1:
        cs.using_lists = True

    cs.add_molecules(smiles_set, list_memberships)

    return cs