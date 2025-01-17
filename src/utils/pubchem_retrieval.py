import pubchempy
import os
import re


def get_Mary_data():
    """
    The code used to retrieve a molecules SMILES string from the PubChem database using its CID
    """

    lipids = [440885, 46173186, 46173313, 50909837, 46173444, 25246183, 25246208, 25200834, 46173394, 25246224, 3083382,
              25202130, 25244473, 51351771, 51351791, 46173153, 50909852, 46173409, 127960, 13831140, 440886, 160295]

    thymines = [1135, 452067, 451954, 6439704, 135557, 6450954, 452068, 452063, 452715, 607039, 196362, 6477693, 370631,
                406592, 6477695, 6477692, 492029, 465896, 445213, 495405, 374907, 457298, 485384, 452946]

    benzothiadiazines = [62940, 43148, 22425, 2910, 2348, 2343, 2122, 12933, 5560, 4870, 4121, 44154271, 44152755, 44151672,
                         6441852, 871720, 216293, 198367, 194167, 188359, 174783, 173791, 107748, 72070, 71652]

    adenines = [3083432, 3083316, 3082029, 3081390, 3080770, 3080762, 3036950, 703739, 466837, 465383, 440867, 25246029,
                15938965, 12358355, 7059571, 7058055, 41211, 34768, 32014, 10238, 9687, 6083, 6076, 1913, 224, 50909893,
                44134557, 25244014, 25201135, 23615303, 23615194, 23421209, 22848660, 16078938, 9578273, 9589376, 6992262,
                6452236, 6449870, 6426627, 5748329, 5491933, 5399013, 4617095]

    glucosinolates = [23682211, 25244590, 46173878, 9548633, 9548621, 9576240, 9576416, 656498, 9548605, 656539, 656538,
                      9576241, 46173877, 44237373, 656547, 46173876, 44237368, 5281133, 656555, 656557, 46173880, 44237257,
                      656543, 656523, 44237260, 656527, 6325242, 9548892, 25244892, 441524, 46173879, 6442557, 5281135,
                      656568, 46173882, 9576738, 46173875, 656545, 656525, 7098673, 17756749, 656541, 656531, 46173881,
                      44237206, 5281136, 44237203, 9548619, 5281138, 46173884, 25245521, 6443008, 5281134, 656537, 25244538,
                      25244201, 17756744, 5281139, 46173883, 25246161, 25245774, 25244220, 25243874, 9548618, 656562, 656548]

    guanines = [764, 374910, 160219, 129161, 133387, 161069, 145817, 129136, 406591, 130450, 478537, 25082899, 471293,
                485625, 485629, 485626, 485624, 195385]

    peroxides = [5497123, 5464098, 641668, 637882, 45266618, 16760624, 16219283, 5311493, 445049, 1035, 45027791, 45027789,
                 44202131, 44145773, 25245484, 25244877, 25244708, 23690934, 22169438, 16394563, 6476300, 6543478, 6450800,
                 6454765]

    cytosines = [455597, 5276954, 597, 452713, 441224, 374908, 492030, 492031, 500131, 6473860, 471292, 477169, 477168,
                 467421, 455604, 455603, 455602, 455598, 16727509, 477170, 455605, 455601]

    uraciles = [1174, 6194, 18323, 9412, 5386, 68216, 13712, 5899, 5360852, 5282192, 688297, 6971263, 6029, 69672, 13268,
                208432, 453162, 1177, 452948, 125110, 3067786, 55281, 456513, 54929]

    chlorothiazides = [127085, 116034, 107748, 71656, 62940, 2348, 2720, 5560, 3639, 50987261, 44151672, 44147212, 24847808,
                       23717274, 11354874, 3083286, 3083063, 242921, 216293, 193444, 188359, 174783, 172393, 159328]

    pregnadienes = [63043, 63042, 63041, 62961, 45006164, 44266812, 24867475, 20054915, 23671691, 11957468, 11954369, 16490,
                    9793, 9782, 9642, 5876, 9571040, 6714002, 6713977, 6452749, 5388957, 5388959, 656804, 633091, 443958,
                    443936]

    acetylaminofluorenes = [5897, 22469, 5896, 168033, 135827, 130776, 130694, 119334, 108117, 22722, 19347, 17270]

    compounds = {'lipids': lipids, 'thymines': thymines, 'benzothiadiazines': benzothiadiazines, 'adenines': adenines,
                 'glucosinolates': glucosinolates, 'guanines': guanines, 'peroxides': peroxides, 'cytosines': cytosines,
                 'uraciles': uraciles, 'chlorothiazides': chlorothiazides, 'pregnadienes': pregnadienes,
                 'acetylaminofluorenes': acetylaminofluorenes}

    for cids in compounds:
        writer = open(cids + '.txt', mode='wb')
        for cid in compounds[cids]:
            c = pubchempy.Compound.from_cid(cid)
            smiles = c.canonical_smiles
            writer.write(str(smiles) + '\n')
        writer.close()

def get_pubchem_by_name(dir):
    files = [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f))]
    for in_file in files:
        #print in_file
        file_name, extension = os.path.splitext(str(in_file))
        if "smiles" not in file_name:
            out_file = file_name+"_smiles"+extension
            print "processing:",in_file,"writing to:", out_file
            with open(dir+in_file,'r') as input, open(dir+out_file,'w') as output:
                mol_info = [re.split(',|:',line)[-2:] for line in input]
                #print mol_info
                smiles_string = ""
                for molecule in mol_info:
                    query_results = pubchempy.get_compounds(molecule[0],"name")
                    if len(query_results) != 1:
                        print "Query for:",molecule[0], "yielded",len(query_results), "results"

                        if len(query_results) > 1:
                            print "Using 1st result: ", query_results[0].canonical_smiles
                    if len(query_results) > 0:
                        smiles_string += query_results[0].canonical_smiles+'\n'

                output.write(smiles_string)


get_pubchem_by_name("evaluation_data/fragmentation_and_lists/molecules/")
