# Characteristic Substructure of Metabolites Application

This project is an implementation of a [method](http://drops.dagstuhl.de/opus/volltexte/2012/3715/)([PDF](http://drops.dagstuhl.de/opus/volltexte/2012/3715/pdf/4.pdf)) that should help identify [metabolites](https://en.wikipedia.org/wiki/Metabolite). To do this it requires:
* A [fragmentation pattern](https://en.wikipedia.org/wiki/Fragmentation_(mass_spectrometry)) of a molecule.
* Potential candidate molecules obtained by performing a lookup with the spectral data on a major database/search engine of chemicals such as [HMDB](http://www.hmdb.ca/metabolites) or [PubChem](https://pubchem.ncbi.nlm.nih.gov/).

It will then produce a _Characteristic Substructure_ that is a representative molecule of all the input molecules. This can then be optionally fed into a tool such as [CFM-ID](http://cfmid.wishartlab.com/) that breaks the molecule back down into the fragmentation pattern via machine learning & heuristics: Ideally, if the original pattern matches this one, then the algorithm has produced a good representation.

As this is a lot of information here is a visualisation of how the application is intended to be used.

![Alt text](/readme_img/research_context_v2.png?raw=true "Optional Title")

More details to how the algorithms work and the general flow of the application are detailed [here]([a relative link](docs/characteristic_substructure.pdf))

## Getting Started

This section explains the requirements of the application and how to get it running.

### Prerequisites

The application has the following dependencies:
* [Python 2](https://www.python.org/downloads/) - it is the main language this is implemented in.
* [Python Enum](https://pypi.python.org/pypi/enum34) - needed for compatibility purposes.
* [NetworkX](https://networkx.github.io/) - a graph library that is used to create the Characteristic Substructure.

#### Optionally
* [MatplotLib](https://matplotlib.org/users/installing.html), [numpy](https://www.scipy.org/scipylib/download.html) and [rdkit](http://www.rdkit.org/) are required to draw molecules.
* [PubChemPy](https://pubchempy.readthedocs.io/en/latest/guide/install.html) - required if you want to do lookups on the PubChem database.
* [CFM-ID](http://cfmid.wishartlab.com/) - used to create a fragmentation pattern from a molecule. An older version is already included by default for testing purposes (is a Windows binary, so will probably not work with Linux distributions).

### Installing

To be able to install the application the aforementioned dependencies are required.

From there, all that is needed is to download the `src` folder and run `src/main.py`.

## Running

The application is given inputs in the command line and has 2 modes:

```
main.py [-h] {cs,rm}

{cs,rm}
  cs        Find characteristic substructure and optionally use
            fragmentation comparison.
  rm        Find a the best-matching molecule from a list, when building 2
            CS.
```

The `cs` creates a characteristic substructure from a molecule list and will draw it if the given libraries are installed.
The `rm` is an experimental process that allows the closest fitting molecule between several lists of candidate molecues via CS.
### Examples

There are some example files of metabolites & spectral patterns that can be used as inputs in `evaluation_data` and results are output into `output_data`.

## Contributing

Please contact [me](https://github.com/NiklasZ) if you want to contribute to this project. There are possible enhancements in regard of the characteristic substructure and heuristic choices.


## Authors

* [NiklasZ](https://github.com/NiklasZ)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* The University of Glasgow's [Computing Department](https://www.gla.ac.uk/schools/computing/)
* [CFM-ID](http://cfmid.wishartlab.com/)

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
