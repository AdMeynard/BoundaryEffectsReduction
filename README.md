## Overview
This repository contains all the files necessary to perform the algorithm **BoundEffRed** (Boundary Effects Reduction). This is an an efficient forecasting approach for the real-time reduction of boundary effects in time-frequency representations.

![](Results/RTexample.gif)

***WARNING:*** This algorithm has been designed to work optimally with MATLAB R2020a. We do not guarantee that it will work with earlier versions.

## Contents

The folder `signals` contains some biomedical and synthetic signals used to implement the above-mentioned algorithm. The following folders are specifically related to the implementation of BoundEffRed:

* `Algorithm` contains the MATLAB functions enabling the implementation of BoundEffRed, using differents extensions schemes such as SigExt, EDMD, or GPR (see paper for details).
* `TimeFrequencyScaleRep` contains the functions generating Time-Frequency and Time-Scale representations such as STFT, Synchrosqueezing Transform, Reassignement, or ConceFT.
* `Scripts` contains the scripts to perform the experiments detailed in the paper. Some corresponding results are directly provided in subfolder `Results`.

The folder `Notes` contains the paper and the associated Supplementary Material.

# Authors

Authors: Adrien Meynard and Hau-Tieng Wu  
Contact email: adrien.meynard@duke.edu
