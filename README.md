# PidTrd

## General information

The PidTrd framework allows to identify hadrons based on the
Trd-dEdx. It assigns a probability to be a certain particle species to
every track depending on its Trd-dEdx and its momentum. 

## Pre-requirements

### Root

ROOT6 is needed for installation:

https://root.cern/install/build_from_source/

Follow instructions
    
### AnalysisTree (optional)

AnalysisTree is needed for usage the interface based on AnalysisTree.

https://github.com/HeavyIonAnalysis/AnalysisTree

Follow instructions

## Installation

Clone PidTrd

    git clone git@github.com:SusanneGlaessel/PidTrd
    
Source ROOT

    source /path-to-root/install/bin/thisroot.sh
    
Export AnalysisTree libraries

    export AnalysisTree_DIR=/path-to-analysistree/install/lib/cmake/AnalysisTree
    
Install PidTrd
    
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/path-to-install-trdpid /path-to-source-trdpid
    make -j install
    
## First run

Step 1 and 2 produce the mc-histograms and probabilities that are
required to identify the tracks. These two steps only need to be
performed once.

Step 3 identfifies the track from an analysistree and writes them into
a new analysistree.

### 1) Creating the MC-input histograms

    Run path-to-trdpid/build/interface/create_mcinput filelist.sh
    outputpath mcfilename jobId

The input for the creation of the MC-histrograms must be an
analysistree including mc-information in the format jobId.analysistree.root.

The path and name of the analysistree is read in via
filelist.txt. Only one analysistree is read in one run. A separate
filelist is required for every run.

For valid MC-histograms a set of run needs to be performed and the
resulting files for every run need to get merged into file with:

    hadd haddId.outfilename.root *. outfilename.root

### 2) Creating the Getter with probabilities

    Run path-to-trdpid/build/interface/create_getter path mcfilename haddId

### 3) Identifying tracks

    Run path-to-trdpid/build/interface/main filelist.sh outputfile
    truncation_mode probability_mode min_hits

The path and name of the analysistree that should be analysed is read in via
filelist.txt.

truncation_mode: Modes for calculation of energy loss for up to 4 layers:
=0: <dEdx> average over all hits
=1-4: Select hits with lowest dEdx:
=1: 1 hit, =2: 2 hits, =3: 3 hits, =4: 4 hits

probability_mode: 
=0: total probability - probability based on particle multiplicites
=1: likelihood - probability based on dEdx-distribution of particle
species

min_hits: Min. number of required hits per track
