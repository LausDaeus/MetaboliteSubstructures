import subprocess
import os

cfm_directory = os.path.abspath(os.path.join(os.path.join(__file__, os.pardir),os.pardir))+"\\fragmentation"
def_param_file = cfm_directory+"\\param_output0.log"
def_config_file = cfm_directory+"\\param_config.txt"

def call_cfm_peak_annotation(smiles, spectrum_file, molecule_id, ppm_mass_tol=10, abs_mass_tol=0.01,
                           param_file=def_param_file,
                           config_file=def_config_file):

    output =  subprocess.Popen([cfm_directory+"\\cfm-annotate.exe", smiles,
                            spectrum_file, molecule_id,
                            str(ppm_mass_tol),
                            str(abs_mass_tol),
                            param_file,
                            config_file], stdout=subprocess.PIPE).communicate()[0]

    return output

def main():
    pass
if __name__ == '__main__':
    main()