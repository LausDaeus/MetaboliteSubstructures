import argparse
import sys
import os

#This adds the necessary directories to find the other modules.
current = os.getcwd()
outer = os.path.dirname(os.getcwd())
sys.path.append(current)
sys.path.append(outer)



import src.utils.verbose as verbose
import src.utils.fileIO as io
import src.parsing.cfm_annotate_parser as cfm_parse
from src.cs_algorithms.characteristic_substructure import CSAlgorithm
import src.cs_algorithms.representative_molecule as rep_mol

#Declare output folder.
output_directory = os.path.abspath(os.path.join(os.path.join(__file__, os.pardir), os.pardir))+"/output_data"

def parse_arguments():
    parser = argparse.ArgumentParser(description="This application finds characteristic substructures, "
                                                 "matches them to fragmentation spectra,"
                                                 "and finding a representative molecule from multiple lists.")
    sub = parser.add_subparsers()
    cs_parser = sub.add_parser("cs", help="Find characteristic substructure and optionally use fragmentation comparison.")
    cs_parser.add_argument("smiles_files", help="1 or more files containing SMILES strings.", nargs="+")
    cs_parser.set_defaults(which='cs')
    cs_parser.add_argument("-r", "--representative_struct", action="store_true",
                        help="Output all of the representative path structures as well.")
    cs_parser.add_argument("-t", "--threshold", type=float, metavar="THRESHOLD",
                        help="the relative frequency of structures in the molecules. Choose percentage from  0..1. "
                             "\nDefault is 0.8.")
    cs_parser.add_argument("-lt", "--list_threshold", type=float, metavar="THRESHOLD",
                        help="relative frequency of structures among a set of lists. Choose percentage from  0..1."
                         "\nDefault is 0.8.")
    cs_parser.add_argument("-iso", "--isomorphism_factor", type=float, metavar="THRESHOLD",
                    help="sets how many path structures need to be present for adding multiple structures." \
                         " Choose percentage from  0..1. \nDefault is 0.8.")
    cs_parser.add_argument("-f", "--fragmentation_pattern", type=str, metavar="FILE_NAME",
                        help="Add a file for a fragmentation pattern of the format specified by the README. Will be "
                             "tested against the deduced characteristic substructure and"
                             "a peak graph + fragments will be returned.")
    cs_parser.add_argument("-nr", "--not_random", action="store_true",
                        help="Ensures tie-breaking is consistent by using a secondary ordering of arbitrary IDs,"
                             " dependent on molecule generation order.")
    cs_parser.add_argument("-ls","--length_start",type=int, metavar="LENGTH", help="Specify starting path length of CS"
                                                                                     " algorithm.")
    cs_parser.add_argument("-le","--length_end",type=int, metavar="LENGTH", help="Specify finishing path length of CS"
                                                                                     " algorithm.")
    cs_parser.add_argument("-s","--step",type=int, metavar="LENGTH", help="Specify the stepping length for"
                                                                             " the algorithm, once it has found a CS.")
    cs_parser.add_argument("-cfm", "--cfm_min_score", type=float, metavar="SCORE",
                        help="Specify a minimum score to highlight matched fragments. \nDefault is 10.")
    cs_parser.add_argument("-img", "--image_of_CS", action="store_true",
                        help="Select if an image of the characteristics substructure should be outputted.")
    cs_parser.add_argument("-o", "--output_name", type=str, metavar="OUTPUT_FOLDER_NAME",
                        help="Add a folder name to use for the output folder. "
                             "Will otherwise default to first SMILES file name.")
    cs_parser.add_argument("-v", "--verbose", action="store_true",
                        help="prints additional info & statistics while processing.")

    molecule_parser = sub.add_parser("rm", help="Find a the best-matching molecule from a list, when building 2 CS.")
    molecule_parser.set_defaults(which='rm')
    molecule_parser.add_argument("smiles_files", help="2 or more files containing SMILES strings. "
                                                      "Number of files has to be even.", nargs="+")
    molecule_parser.add_argument("molecule_check_list", type=str, metavar="LIST_FILE_NAME",
                        help="File containing list of SMILES strings to try adding to the 2 CS as individual molecules.")
    molecule_parser.add_argument("-t", "--threshold", type=float, metavar="THRESHOLD",
                        help="the relative frequency of structures in the molecules. Choose percentage from  0..1. "
                             "\nDefault is 0.8.")
    molecule_parser.add_argument("-lt", "--list_threshold", type=float, metavar="THRESHOLD",
                        help="relative frequency of structures among a set of lists. Choose percentage from  0..1."
                         "\nDefault is 0.8.")
    molecule_parser.add_argument("-nr", "--not_random", action="store_true",
                        help="Ensures tie-breaking is consistent by using a secondary ordering of arbitrary IDs,"
                             " dependent on molecule generation order.")
    molecule_parser.add_argument("-ls","--length_start",type=int, metavar="LENGTH", help="Specify starting path length of CS"
                                                                                     " algorithm.")
    molecule_parser.add_argument("-le","--length_end",type=int, metavar="LENGTH", help="Specify finishing path length of CS"
                                                                                     " algorithm.")
    molecule_parser.add_argument("-s","--step",type=int, metavar="LENGTH", help="Specify the stepping length for"
                                                                             " the algorithm, once it has found a CS.")
    molecule_parser.add_argument("-o", "--output_name", type=str, metavar="OUTPUT_FOLDER_NAME",
                        help="Add a folder name to use for the output folder. "
                             "Will otherwise default to first SMILES file name.")
    molecule_parser.add_argument("-v", "--verbose", action="store_true",
                        help="prints additional info & statistics while processing.")
    molecule_parser.add_argument("-img", "--images_of_molecules", action="store_true",
                        help="Select if an image of the chosen molecule, CS1 and CS2 should be outputted.")

    return parser.parse_args()

def main():
    """
    Creates the algorithm object with the correct parameters and calls methods based on the command line arguments.

    Uses the argparse module to parse the command line arguments that are given
    The file name for the SMILES file is mandatory, the flags and threshold specification are optional
    :return: None
    """


    global output_directory
    #Argument Parsing
    args = parse_arguments()

    #Configuring options
    if args.verbose:
        verbose.initialise()

    if args.output_name:
        output_directory +="/"+args.output_name
    else:
        file_path, extension = os.path.splitext(args.smiles_files[0])
        output_directory += "/"+os.path.basename(file_path)

    if args.which == "rm":
        rep_mol.find_most_representative_molecule(args, output_directory)
        return

    cs = CSAlgorithm(**vars(args))

    smiles_set, list_memberships = io.read_smiles(args.smiles_files)
    cs.list_number = len(set(list_memberships))

    if cs.list_number > 1:
        cs.using_lists = True

    verbose.info("total_time_start")
    cs.add_molecules(smiles_set, list_memberships)

    verbose.info("start_paths")
    paths = cs.find_graphs_paths()
    verbose.info("finish_paths")

    if args.representative_struct:
        all_structures = cs.find_all_representative_structures(paths)
        print 'Representative Structures:'
        if all_structures:
            io.check_dir(output_directory)
            structure_frequency_list = [molecule.to_smiles()+" "
                                        +str(all_structures[molecule]) for molecule in all_structures]
            for structure in structure_frequency_list:
                print structure
            io.write_representative_structure(structure_frequency_list,output_directory+'/representative_structures.txt')
        else:
            print "None found."


    c_structure = cs.find_characteristic_substructure(paths)

    #Display the structures in the command line
    print 'Characteristic Substructure:'
    if c_structure:
        smiles = c_structure.to_smiles()
        print smiles

        io.check_dir(output_directory)
        io.write_to_file(smiles, output_directory + '/cs.txt')

        if args.image_of_CS:
            import src.drawing.draw_molecule as draw
            draw.draw_smiles(smiles, output_directory+"/cs.png")

        if args.fragmentation_pattern:
            import src.fragmentation.fragmentation_matcher as frag
            import src.drawing.plot_spectra as plot
            result = frag.call_cfm_peak_annotation(smiles, args.fragmentation_pattern, "CS")
            io.write_to_file(result, output_directory+'/cfm_output.txt')
            print result
            cfm_dict = cfm_parse.parse_cfm(io.read_file(output_directory+'/cfm_output.txt'))
            spec_dict = io.read_spectra(args.fragmentation_pattern)

            if args.cfm_min_score:
                plot.plot_cfm_data(cfm_dict, spec_dict, output_directory, args.cfm_min_score)
            else:
                plot.plot_cfm_data(cfm_dict, spec_dict, output_directory)
    else:
        print "None found."
    verbose.info("total_time_finish")

if __name__ == '__main__':
    main()