# CIM-2017

This is the R code associated to the article: 

Rodríguez-Fernández, Víctor, Antonio Gonzalez-Pardo, and David Camacho. “Modelling Behaviour in UAV Operations Using Higher Order Double Chain Markov Models.” IEEE Computational Intelligence Magazine. In press (2017)

You can find a link to the article here:
[TBD]

## Requirements

1. Download *R 3.4.1* and *RStudio* (a version of RStudio supporting *R 3.4.1* is needed)
2. Open the file *CIM-2017.RProj* with RStudio
3. Type `install.packages("packrat")` to install package *packrat*, for managing dependencies
4. Type `packrat::restore()` to load all the project dependencies from remote repositories. This process will take several minutes.
5. Type `install.packages("./MmgraphR_0.1.tar.gz",repos = NULL, type = "source")` to install the package ``"MmgraphR"
from a local source tarball.

## Usage

The project structure has been created using the R package *ProjectTemplate*. For more details about ProjectTemplate, see http://projecttemplate.net

* The input data can be found in the `data` folder
* Every analysis script can be found in the folder `src`.
* The output data from the scripts is stored in the `output` folder

The scripts found under the `src` folder are named with numerical prefixes (01-, 02-, ...) to indicate the proper order of execution in order to reproduce the results of the paper. Below are briefly described the contents of each script:

* `1-Markovian_Comparison-InputData.R`: Load the data from the `data` folder to be used in the rest of the scripts.
* `2-Markovian_Comparison-DCMM`: train and evaluate DCMMs using the march package.
* `3-Markovian_Comparison-Rank_Aggregation`: performs a Rank Aggregation process over a grid of evaluated DCMMs.
* `4-Markovian_Comparison-CIM-All_in_one_Rank_Aggregation`: Study of the predictability/interpretability importance in the rank aggregation.

## Contact

For any questions about the use of the code please contact me by email or add a new issue in this repository.