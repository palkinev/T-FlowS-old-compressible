"""
A freely-propagating, premixed hydrogen flat flame with multicomponent
transport properties.
"""

import cantera as ct
import numpy as np

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'CH4:0.375, O2:1, N2:3.76'  # premixed gas composition
width = 0.03  # m
loglevel = 1  # amount of diagnostic output (0 to 8)
tol_ss    = [1.0e-13, 1.0e-14]        # [rtol atol] for steady-state problem
tol_ts    = [1.0e-13, 1.0e-14]        # [rtol atol] for time stepping


# IdealGasMix object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture
gas = ct.Solution('gri30.xml')
gas.TPX = Tin, p, reactants
initial_grid = [0.0, 0.001, 0.01, 0.02, 0.03, 0.04, 0.049, 0.05] # m

# Set up flame object
f = ct.FreeFlame(gas, grid=initial_grid, width=width)
# Set solver tolerance
f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)

# inlet parameters
#f.inlet.mdot = mdot_f
f.inlet.X=reactants
f.inlet.T=Tin
f.energy_enabled = False
f.set_refine_criteria(ratio=3, slope=0.1, curve=0.2)
f.show_solution()

# Solve with mixture-averaged transport model
f.transport_model = 'Mix'
f.solve(loglevel=loglevel, auto=True)

# Solve with the energy equation enabled
f.save('h2_adiabatic.xml', 'mix', 'solution with mixture-averaged transport')
f.show_solution()
print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.u[0]))

# Solve with multi-component transport properties
f.transport_model = 'Multi'
#f.soret_enabled = True #  enable calculation of the Soret diffusion term
f.solve(loglevel) # don't use 'auto' on subsequent solves
f.show_solution()
print('multicomponent flamespeed = {0:7f} m/s'.format(f.u[0]))
f.save('h2_adiabatic.xml','multi', 'solution with multicomponent transport')

# write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('h2_adiabatic.csv', quiet=False)
