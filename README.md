#Characteristic Substructure Application

If you're interested in using this/working on it, message me and I shall explain the details.

##Description:

The application can take the characteristic substructure (CS) of a list of 1-n SMILES strings as well as take in 1 fragmentation 
spectrum to compare to the CS. Additionally it may also take in 2 sets of list of SMILES strings and 1 single list to determine
the most correlating molecule in the single list.

##Execution
    All features of the application can be run through src/main.py.

##Input
    The following input categories are supported:
    python main.py cs <files and flags>
    python main.py rm <files and flags>

    where:
    cs = Characteristic substructure
    rm = Representative molecule

    ##cs inputs:
    To receive help:
    python main.py cs -h

    To create

    There is a range of commands, all of which are detailed in the respective paper. Here are some examples of what can be run:
    python main.py cs -h                                                             ###For help on characteristic substructure
    python main.py rm -h                                                             ###For help on representative molecules
    python main.py cs -v ..\test_data\acetylaminofluorenes.txt                       ###Example of CS algorithm
    python main.py cs -img ..\test_data\acetylaminofluorenes.txt ..\test_data\benzothiadiazines.txt -o
    ###A multi list example, saved to output_data\demonstration



##Data

##Output

    All outputs are saved to output_data. The output will always be placed in a folder whose name is determined by
    either the first SMILES file or any name supplied with the -o flag. E.g for:
    "python main.py cs ..\test_data\acetylaminofluorenes.txt", folder name = acetylaminofluorenes
    "python main.py cs -o example ..\test_data\acetylaminofluorenes.txt" folder name = example

##Dependencies

    NetworkX is required to run the CS algorithm
    MatplotLib, Numpy are required if the user wishes to draw the molecule
    PubChemPy is required to run the pubchem_retrieval module
    rdkit is required to plot molecules.

##Environment
    
    Has been tested in Windows 8.1-
    CFM-ID software may not work with Linux distributions
