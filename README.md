# Characteristic Substructure of Metabolites Application

This project is an implementation of a [method](http://drops.dagstuhl.de/opus/volltexte/2012/3715/)([PDF](http://drops.dagstuhl.de/opus/volltexte/2012/3715/pdf/4.pdf)) that should help identify [metabolites](https://en.wikipedia.org/wiki/Metabolite). To do this it requires:
* A [fragmentation pattern](https://en.wikipedia.org/wiki/Fragmentation_(mass_spectrometry)) of a molecule.
* Potential candidate molecules obtained by performing a lookup with the spectral data on a major database/search engine of chemicals such as [HMDB](http://www.hmdb.ca/metabolites) or [Pubchem](https://pubchem.ncbi.nlm.nih.gov/)

It will then produce a _Characteristic Substructure_ that is a representative molecule of all the input molecules. This can then be optionally drawn as well as fed into a [tool](http://cfmid.wishartlab.com/) that breaks the molecule back down into the fragmentation pattern. Ideally, if the original pattern matches this one, then the algorithm has produced a good representation.

As this is a lot of information here is a visualisation of how the application is intended to be used.

![Alt text](/readme_img/research_context_v2.png?raw=true "Optional Title")

More details to individual aspects follow in further sections.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags).

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc


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
