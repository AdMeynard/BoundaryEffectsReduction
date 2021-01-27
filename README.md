## Overview
This repository contains all the files necessary to perform the algorithm **BoundEffRed** (Boundary Effects Reduction). This is an efficient forecasting approach for the real-time reduction of boundary effects in time-frequency representations.

<table style="width:100%">
  <tr>
    <th width="50%"><img src="Animations/WithoutBoundEffRed.gif" width=450 height=326></th>
    <th width="50%"><img src="Animations/WithBoundEffRed.gif" width=450 height=326></th>
  </tr>
  <tr>
    <th width="50%">Ordinary syncrosqueezing transform of a cardiac signal</th>
    <th width="50%">Boundary-free syncrosqueezing transform of the same signal, obtained via BoundEffRed</th>
  </tr>
 </table>

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
