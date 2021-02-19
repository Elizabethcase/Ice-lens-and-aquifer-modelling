This folder contains the code to run a simulation for the evolution of firn temperature, porosity, and saturation.

The main code is `RunCode.m’ which does the time-stepping, applies the boundary conditions, and calls the functions that are required to do portions of the computation. The current configuration solves the nondimensional equations described in Appendix B for 10 units of time with a surface flux of -0.5 and an accumulation rate of 0.09. The various parameters are described and comments on lines 19-26. 

In the plots that appear when the code is run there are many lines that represent the nondimensional quantities:

Yellow : temperature
Black : enthalpy
Green : porosity
Blue : total water
Red : saturation
Cyan : water flux
Magenta : water pressure

Certain variables within `RunCode.m’ that may be of interest:

T : total simulation time
plot_amount : the number of timesteps between plots
save_amount : the number of timesteps between solution saves
AccumulationRate : the nondimensional meters of ice accumulated per year (or (1-phi_0) times the meters of snow accumulated)
Qbar : constant nondimensional surface flux
type : the style of compaction. The options are a pore-closure model or the Herron-Langway 	empirical model. 

