## Overview
This repository contains all the files necessary to perform the algorithm **BoundEffRed** (Boundary Effects Reduction). This is an efficient forecasting approach for the real-time reduction of boundary effects in time-frequency representations.

![](Animations/WithBoundEffRed.gif) ![](Animations/WithoutBoundEffRed.gif)

***WARNING:*** This algorithm has been designed to work optimally with MATLAB R2020a. We do not guarantee that it will work with earlier versions.

## Contents

The folder `Signals` contains some biomedical and synthetic signals used to implement the above-mentioned algorithm. The following folders are specifically related to the implementation of BoundEffRed:

* `Algorithm` contains the MATLAB functions enabling the implementation of BoundEffRed, using differents extensions schemes such as SigExt, EDMD, or GPR (see paper for details).
* `TimeFrequencyScaleRep` contains the functions generating Time-Frequency and Time-Scale representations such as STFT, Synchrosqueezing Transform, Reassignment, or ConceFT.
* `Scripts` contains the scripts to perform the experiments detailed in the paper. Some corresponding results are directly provided in subfolder `Results`.

The folder `Paper` contains the article and the associated Supplementary Materials.

# Authors

Authors: Adrien Meynard and Hau-Tieng Wu  
Contact email: adrien.meynard@duke.edu
