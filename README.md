[![Travis Build Status](https://travis-ci.org/tuxette/RNAseqNet.svg?branch=master)](https://travis-ci.org/tuxette/RNAseqNet)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/tuxette/rnaseqnet?branch=master&svg=true)](https://ci.appveyor.com/project/tuxette/rnaseqnet)
[![Coverage Status](https://img.shields.io/codecov/c/github/tuxette/RNAseqNet/master.svg)](https://codecov.io/github/tuxette/RNAseqNet?branch=master)

# RNAseqNet

Infer log-linear Poisson Graphical Model with an auxiliary dataset. Hot-deck
multiple imputation method is used to improve the reliability of the inference
ith an auxiliary dataset. Standard log-linear Poisson GM can also be used for 
the inference and the StARS criterion is implemented to drive the selection of
the regularization parameter.

See [[Citation file]](./inst/CITATION) for citation details.
