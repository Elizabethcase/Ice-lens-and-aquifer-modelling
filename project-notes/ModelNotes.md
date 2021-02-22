### Model code

^c629cb
 The code is described in the README.txt in the code folder in projects/Meyer and Hewitt Model
#### ReadMe

    This folder contains the code to run a simulation for the evolution of firn temperature, porosity, and saturation.
    
    The main code is `RunCode.m` which does the time-stepping, applies the boundary conditions, and calls the functions that are required to do portions of the computation. The current configuration solves the nondimensional equations described in Appendix B for 10 units of time with a surface flux of -0.5 and an accumulation rate of 0.09. The various parameters are described and comments on lines 19-26.
    
    In the plots that appear when the code is run there are many lines that represent the nondimensional quantities:
    
    Yellow : temperature Black : enthalpy Green : porosity Blue : total water Red : saturation Cyan : water flux Magenta : water pressure
    
    Certain variables within `RunCode.m` that may be of interest:
    
    T : total simulation time plot\_amount : the number of timesteps between plots save\_amount : the number of timesteps between solution saves AccumulationRate : the nondimensional meters of ice accumulated per year (or (1-phi\_0) times the meters of snow accumulated) Qbar : constant nondimensional surface flux type : the style of compaction. The options are a pore-closure model or the Herron-Langway empirical model.
    
#### Descriptions of each file

1.  `AdvectiveFlux.m`

- The advection of a field due to vertical velocity
	- if velocity is positive, receives field from previous vertical step
	- If velocity is negative, receives field from 
        
        Notes: possibly these are on offset grids, which would make this make a bit more sense
        
        **Questions**
        
        1.  Not sure what fadvp vs. fadvm are other than offset? And u field.. what does that =?
    2.  ` CompactionFunction.m `
        
        -   Function takes in W (total water e.g. sum of ice and liquid water), phi (snow porosity), and theta (snow temperature) and outputs compaction
        -   Options for compaction:
            -   `none` (zeros)
            -   `poreclosure`
                -   $\frac{\partial u}{\partial z}=2A\frac{\phi}{1-\phi}\sigma^n$ (Glen's flowlaw)
                
    3.  `conversiontoemthalpy.m`
        
    4.  `conversiontotemperature.m`
        
    5.  `EnthalpyDiffusiveFlux.m`
        
    6.  `FullySaturatedWaterPressure.m`
        
    7.  `SaturatedDiffusiveFlux.m`
        
    8.  `TemperatureDiffusiveFlux.m`

### Surface Energy Balance and Surface Mass Balance

To understand model, need to understand two things:

1.  **boundary conditions** —> what are the appropriate surface energy magnitudes & surface mass balances for JIRP over the course of a year?
2.  **non-dimennsionalization** —> how is each term non-dimensionalized, and how is this reflected in the figures
    [[@abe-ouchiInsolationdriven100000year2013]]
    -   **Surface Energy Balance**
        
        $$(1) \qquad -\bar{K}\frac{\partial T}{\partial z} = - (1-\alpha)S_w - L_w + \epsilon \sigma T^4 - \chi (T_a-T) -\rho_w c_i a (T_a -T) - \rho_w c_w R (T_a - T) + \rho_w \mathcal{L} M $$
        
        -   _Table of Surface Energy Balance terms & possible citations_

| Term                                    | Description                                                                                                    | Citations                                                                                               |
| --------------------------------------- | -------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| $-\bar{K}\frac{\partial T}{\partial z}$ | conduction, $(\bar{K} = thermal conductivity, (1-\phi)K)$                                                      | Quakenbush T and Wendler G Ultraviolet CA) and Short-Wave Radiation on the Juneau Icefield, Alaska. , 6 |
| $-L_w$	                              | incoming longwave radiation                                                                                    |                                                                                                         |
| $(1-\alpha)S_w$                         | incoming shortwave radiation, $(\alpha$ = albedo, $S_w$  = shortwave radiation)                                |                                                                                                         |
| $\epsilon \sigma T^4$                   | outgoing longwave radiation, $(\epsilon$ = emissivity, $\sigma$ = Stefan-Boltzman constant, $T$ = temperature) |                                                                                                         |
| $\chi (T_a-T)$                          | turbulant heat transfer, ($\chi$ = coefficient)                                                                |                                                                                                         |
| $-\rho_w c_i a (T_a -T)$                | sensible heat flux from solid precipitation ($c_i$ = heat capacity for ice, $a$ = accumulation rate)           |                                                                                                         |
| $\rho_w c_w R (T_a - T)$ |	sensible heat flux from liquid precipitation, ($c_w$  = heat capacity for water, $R$ = rainfall rate)	||
| $\rho_w \mathcal{L} M$|	latent heat flux associated with melting, ($\mathcal{L}$ = latent heat from melting ice, $M$ = melt rate	||

Equation (1) can be [[Linearization | linearized]] around the melting temperature to obtain 
        
$$(2) \qquad Q(t) = (1-\alpha)S_w + L_w - \epsilon\sigma T_m^4 + \chi(T_a - T_m)+\rho_w c_i a (T_a - T_m) + \rho_w c_w R(T_a-T_m)$$
        
We need to figure out what each of these terms should be for JIRP over the course of a year
        
-   **Surface Mass Balance**
        
-   **Non-dimensionalization**
        
Non-dimensionalize by $l = \frac{Q_0 t_0}{\rho \mathcal{L}}$ and $t=t_0$
        
Write Temperature as : $T=T_m + \Delta T \theta$ and scale by $\Delta T = Q_0/h$ where h = $\chi + 4\epsilon \sigma T^3_m$

### Compaction model
    
Struggling with the math. To summarize so far:
    
Firn compacts due to overburden pressure, temperature and grain size. The presence of liquid water only changes densification inasmuch as: 1) where there is water, T = 0, 2) water releases heat when it refreezes, absorbs heat when it melts. Meltwater is an important way the climate plays a role in ice sheet dynamics. Meltwater storage and porosity (e.g. densification) are closely linked — meltwater can only be stored where there is pore space, otherwise generates ponding/runoff.
    
Pieces of the model
    
1.  **Percolation** through porous ice is controlled by
        -   mass conservation — how much ice, water, and air are there, and how does melting contribute to increasing or decreasing amount of ice/water.
        -   Darcy's law — flow of water through pore space)
        -   Permeability — depends on the Kozeny relationship (%VF of pore space and grain size determine how easy(?) it is for water to flow.)
        -   Water pressure
            -   When flow is partially saturated, capillary forces drive flow and water pressure equals capillary pressure
            -   When flow is fully saturated, water pressure is constrained by mass conservation
    
Point of confusion - In section 2.2, authors use HL form of c to write out equation in terms of porosity (% volume fraction of air). Why is this allowed? E.g. they claim
    
$$\frac{\partial\phi}{\partial t} +u_i \nabla \phi = -c\phi$$
    
is the same as
    
$$\frac{\partial(1-\phi)\rho_i}{\partial t} +u_i \nabla (1-\phi)\rho_i = -c(\rho_i - (1-\phi)\rho_i)$$
    
UPDATE: I realized the RHS was actually $c*(\rho_i-\rho)$ in HL so everything works out now. All the $\rho_i$ drop out, and so do the effects of the 1 in the $1-\phi$ on the LHS, but on the RHS, we end up with an extra $c$ term where $c$ is defined as

$$ c = \begin{cases} 11 a \exp{\frac{-1222}{T}} &\text{if } \phi > 0.4 \ 575 \sqrt{a} \exp{\frac{-2574}{T}} &\text{if } \phi \leq 0.4 \end{cases} $$
    
 This is (sic: different?) from [[Herron and Langway - 1980 - Firn Densification An Empirical Model.pdf | Herron and Langway]] ([[@herronFirnDensificationEmpirical1980 | notes]]) 6a and 6b, where the values in the exponent have been adjusted to include R .
 
 [[Densification]]

### Non-dimensionalization
    
standard units -

L = length; T = time; $\theta$ = temperature; M = mass

| Key Parameter | Units        | To non-dimensionalized  | From non-dimensionalized | 
| ------------- | ------------ | ----------------------- | ------------------------ | 
| $\dot b$      | $L / T$      | $\rho \mathcal{L}/Q\_0$ | $Q\_0/\rho \mathcal{L}$  | 
| $Q$           | $ML^2T^{−3}$ |                         |                          |  
| time          | T            | $1/t\_0$                | $t\_0$                   |


I started in JIRP Firn Aquifer/MH-Model/convertoparams.m to try to understand how these ND'd numbers change with JIRP conditions, and how to reproduce a stacked aquifer.

Received some SEB data/a paper suggestion from Seth (see [[McNeil et al. - 2020 - Explaining mass balance and retreat dichotomies at.pdf]] and associated data/ [[@mcneilExplainingMassBalance2020 | notes]]).

I downloaded that data and found the average temperature at C18. I need to find

-  elevation of the divide
-   lapse rate (see paper)
-   to calculate average temperature at the divide

dT here is more like 15. Assuming h (turbulent heat transfer) stays the same (which it might not, I have no idea how this is calculated, this gives a Q0 (bar?) of 222. Not sure if this will change things much.

-    how do you calculate turbulent heat transfer
-   how does the 2 (annual mean forcing) relate to Q0 and where does it show up in the model
    -   is Qbar 400 W m-2 in the model?
    -   what is Qbar at Juneau?
    -   what is Q0 (e.g. maximum change in Qbar) at Juneau?
-   What other info do I need to answer these questions?

$\mathcal{U} = \rho g k\_0 t\_0 / lu$

$S = \mathcal{L}/c\_p \Delta T$

$Pe = \rho c\_p l^2/ K t\_0$

$B = \rho g d\_p l / \gamma$

where

$\rho =\rho\_i = 917\: kg\: m^{-3}$

$g = 9.81\: m\: s^{-2}$