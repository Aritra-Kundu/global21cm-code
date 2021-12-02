# global21cm-code
A simple Python code to generate the global 21 cm signal. Includes the standard cosmology and also the interaction between dark matter and baryons. The later is incorporated in order to explain the amplitude of the EDGES signal (Bowman et al. 2018). Currently, the there are some limitations of the code which we are working on to solve. It uses a polynomial fit to model the redshift evolution of the X-ray emissivity, which we obtained from another separate code. And also the coupling process due to Lyman-alpha photons is absent from the code.

There are two separate codes.

global_21cm.py - Simulates the redshift evolution of the gas kinetic temperature, hydrogen ionization fraction, dark matter temperature and relative velocity between dark matter and baryons.

global_21m_constrain_DM_b.py - Provides the contraint on the mass of dark matter particles and the interaction cross-section.


Future plans:-

1. To include the Lyman-alpha coupling and have the complete standard picture (currenly under development).
2. Model the X-ray heating part within the code itself.
3. Include the effects of primordial magnetic field (currently under development).
4. Create a more user-friendly version of the code.
