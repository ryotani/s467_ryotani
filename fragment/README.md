# Fragment calibration macros

A set of macros to calibrate the reconstruction of PIDs.


## Concept

Using full statistics of data for each target runs, the calibration processes for the position correlation, the Brho reconstruction, and the mass reconstruction are proceeded.


## mktree_tofw_frs.C

The first step is to merge the runs for each secondary target setups into one file. By selecting the position of target based, a merged file can be created.


## fragment_calib_twim.C

** For this process, empty target run should be analysed. **
As the X position at MW3 should be correlate with the position and angle before the GLAD magnet, in addition to the reduced momentum of the particle, the contributions from the position X, estimated from the average of MW1X and MW2X, and the X angle, theta reconstructed with TwinMUSIC, are fitted. The result, MW3X_mod, will be output in the output directory.


## fragment_reco_twim.C

** For this process, empty target run should be analysed. MW3X_mod should be updated according to the previous step. **

Assuming the brho value doesn't change, the linear correlation with MW3X_mod and the brho will be fitted. Whith this, a preliminary PID plot can be obtained. The beta correction yet to be optimised.


## tofw_beta_offset_nofrsgate.C

** This should be run for each target setting. **

Assuming the most of the particle doesn't change the mass and chage, the A/Q value will be corrected paddle-by-paddle of TofW. By applying several condiions as velocity and charges in MUSIC detctors, the probable elastic channel will be gated. Then the correlation with angle of the Twim will be taken to obtain the coefficient for the flight length and the time offset will be fitted. This code is still under development.