[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.32016.svg)](http://dx.doi.org/10.5281/zenodo.32016)

# AlloPred

This is the code and training/testing data used at the AlloPred web server, which predicts allosteric pockets on proteins from a PDB format file:

http://www.sbg.bio.ic.ac.uk/allopred/home

The method is described here: Greener, JG and Sternberg, MJE, AlloPred: prediction of allosteric pockets on proteins using normal mode perturbation analysis, *BMC Bioinformatics*, 2015, 16(335). [Link to paper](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0771-1).

Downloading this code lets you run AlloPred locally. It should work on all Unix-based systems. With modification it could be made to work on Windows.

Please contact Joe Greener (joe.greener13 at imperial.ac.uk) for support, queries and suggestions. Alternatively, open an issue on GitHub.


## Requirements

* Python 2.7 with the [NumPy](http://www.numpy.org/) and [ProDy](http://prody.csb.pitt.edu/) packages installed.
* fpocket v2.0, which can be downloaded from [here](http://fpocket.sourceforge.net/). Follow the installation instructions to compile the executables.
* SVM-light, which can be downloaded from [here](http://svmlight.joachims.org/). Follow the installation instructions to compile the executables.


## Usage

Follow these steps to set up AlloPred - the shell commands are for bash:

1. Download the files and extract them as usual, or clone the repository.

2. The environmental variables `$ALLOPRED_DIR` and `$SVM_LIGHT_DIR` need to be set as the filepaths to the AlloPred directory and the SVM-light directory respectively:
  ```
  export ALLOPRED_DIR=/path/to/allopred/
  export SVM_LIGHT_DIR=/path/to/svm_light/
  ```
  Consider adding these lines to your profile so you don't have to run them every session.

Follow these step to run AlloPred:

1. Obtain a PDB format file (`in_file.pdb`), e.g. from the [Protein Data Bank](http://www.rcsb.org/pdb/home/home.do).

2. Create a one-line file (`act_res.txt`) containing the active site residues of the protein. The format is `10:A,11:B` for residue 10 on chain A and residue 11 on chain B. These can be found using resources such as the [Catalytic Site Atlas](http://www.ebi.ac.uk/thornton-srv/databases/CSA/). An example PDB file and active residue file can be found in the example directory of AlloPred.

3. Run fpocket v2.0 on the PDB file:
  ```
  fpocket -f in_file.pdb
  ```
  This assumes `fpocket` is on the path. This produces the directory `in_file_out`. AlloPred is optimised on the default fpocket parameters but you can change these in accordance with the fpocket documentation if you wish.

4. The following command, from the directory containing `in_file.pdb` and `in_file_out`, runs the AlloPred pipeline:
  ```
  python $ALLOPRED_DIR/run_allopred.py in_file act_res.txt
  ```
  The arguments are the input file prefix and the path to the active site residue file. Running the `run_allopred.py` script with fewer than 2 arguments returns these instructions for the command.

5. The output files are:
  * `in_file.out`: the AlloPred output file containing the input parameters and the values for each pocket in order of AlloPred ranking.
  * `in_file.svm`: the SVM input file in the SVM-light format.


## Other files

* `dataset` contains information on the training and testing sets.
* `example` contains the inputs and outputs of an example run using the PDB entry with ID 1FX2.
* `svm_model.txt` is the optimised SVM built on the whole training set.


## Reproducibility

For reference, here are steps to reproduce the results in Figure 3 of the paper (link above):

1. Obtain PDB files from the test set list in `dataset/test_set.tsv`.
2. Keep only the chains for each file given in `dataset/test_set.tsv`.
3. For each protein carry out steps 2-4 above using the active site residues given in `dataset/test_set.tsv`.
4. Compare the residues for the top predicted pockets in `in_file.out` to the allosteric residues given in `dataset/test_set.tsv`.


## Note

On some systems the normal mode calculation step fails (`More than 6 zero eigenvalues are calculated`). This appears to be due to an error in `scipy.linalg.eigh` which causes an error in the `calcModes` function of ProDy.

Thank you for using AlloPred!
