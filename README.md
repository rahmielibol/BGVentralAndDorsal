# Basal Ganglia Model
from Ventral (limbic) circuit to Dorsal (planning, action selection) circuit.

A Computational Model from Single Cell to Circuit Level

How to cite: Elibol, R., Şengör, N.S. Ventral Striatum Modifies the Activity of Basal Ganglia Circuit: A Computational Model (2022). https://doi.org/10.21203/rs.3.rs-1594417/v1

A computational model is established for Basal Ganglia circuits in Python, using BRIAN2 library. The results are verified the validity of the model by showing the consistency of simulation results with the empirical data. So, to run the code, one should download BRIAN2 (https://briansimulator.org/)

First, single neurons are considered and their model behavior, then synaptic currents are considered. The results are given with single neuron behavior, synaptic currents, raster plots and local field potentials. There are three scenarios, scenario 0 is for testing whether all is working OK, and scenario 1 and scenario 2 stimuli are given and the overall behavior of the nucleus accumbens and BG circiuts are observed.

As there are random variables, the results of two runs would not be exactly same, but very similar and behavior would be consistent.

Further detail is written as comment lines in the code.



