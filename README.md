# NNLOCAL

A parton-level Monte Carlo program for color-singlet production in hadron collisions at NNLO accuracy. For a description of the computational method and code, see [arXiv:2412.21028](https://arxiv.org/abs/2412.21028).

## Authors

* Vittorio Del Duca [Vittorio.Del.Duca@cern.ch]
* Claude Duhr [cduhr@uni-bonn.de]
* Levente Fekésházy [levente.fekeshazy@desy.de]
* Flavio Guadagni [guadagni.flavio@gmail.com]
* Pooja Mukherjee [pooja.mukherjee@desy.de]
* Gábor Somogyi [somogyi.gabor@wigner.hun-ren.hu]
* Sam Van Thurenhout [sam.van.thurenhout@wigner.hun-ren.hu]
* Francesco Tramontano [francesco.tramontano@unina.it]

## Installation

### Dependencies

An installation of [LHAPDF](https://www.lhapdf.org) is required.

### Installation

1. Clone the repository
```
git clone https://github.com/nnlocal/nnlocal.git
```
2. Navigate to the project directory
```
cd nnlocal
```
3. Compile
```
make
```

### Running

1. Prepare the run by editing the input file `./bin/testrun-H/input.DAT`. The most important parameters that are set in this file are described below. 

#### Serial mode

2. Run the exacutable 
```
./bin/nnlocal
```

#### Parallel mode

2. Set the `parallel` flag to `1` in `./bin/testrun-H/input.DAT`

3. Set up the parameters for parallel running by editing the script `./bin/testrun-H/runpar.sh`. The most important parameters in the script are descibed below.

4. Run the script
```
./bin/testrun-H/runpar.sh
```

## Preparing the input

### Basic parameters in `input.DAT`

* `nproc`: the process ID number. Currently only the $pp \to H$ process is implemented, for which `nproc = 710`.
* `order`: the order in the strong coupling relative to the Born process. Hence, `order = 0,1,2` correspond to the LO, NLO and NNLO computations.
* `part`: specifies which part of the full computation to perform. Possible values are the following.
  1. `born`: include all contributions up to the given `order` that have Born-like (i.e., $2 \to 1$) kinematics, e.g., the double virtual contribution at NNLO; 
  2. `virt`: include all contributions up to the given `order` that have Born + one parton kinematics, e.g., the real-virtual contribution at NNLO; 
  3. `real`: include all contributions up to the given `order` that have Born + two parton kinematics, e.g., the double real contribution at NNLO; 
  4. `tota`: include all contributions up to the given `order`.
* `sqrts`: the total center of mass energy of the hadron-hadron collision in GeV.
* `hmass`: the mass of the Higgs boson $m_H$ in GeV.
* `scale`: the renormalization scale $\mu_R$ in GeV.
* `facscale`: the factorization scale $\mu_F$ in GeV.
* `itmx1`: the total number of iterations used for grid refinement in serial running mode.
* `ncall1`: the total number of evaluations per iteration during grid refinement.
* `itmx2`: the total number of iterations for collecting results after the grid has been set up.
* `ncall2`: the total number of evaluations per iteration during result collection.
* `(ncall virt)/(ncall born)`: if `part = tota` is selected, then all partonic contributions are evaluated during the run. In this case, `ncall1` and `ncall2` give the number of evaluations per iteration during the grid refinement and collection stages for the contributions with Born-like kinematics. However it is usually necessary to run higher-multiplicity partonic processes with more points and this parameter provides a way to increase the total number of points by multiplying `ncall1` and `ncall2` with the value set here for contributions with Born + one parton kinematics. If the value of `part` is something other than `tota`, this parameter is inactive.
* `(ncall real)/(ncall born)`: same as above, for Born + two parton kinematics.  
* `parallel`: specify whether to run in serial mode (`0`) or parallel mode (`1`), see below.

### Basic parameters in `runpar.sh`

* `ncores`: the number of CPU cores that the user wishes to use simultaneously.
* `nprocessesgrid`: the total number of instances to be used during each iteration step of grid refinement.
* `nprocessesaccu`: the total number of instances to be used during the collection of results after the grid has been set up.
* `maxgrid`: the number of iteration steps used for grid refinement.
