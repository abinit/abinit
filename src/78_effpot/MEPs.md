# Multibinit Enhancement Proposals

## changes
- add label property to potentials (for debugging)
   For each potential, it has a default label, which can be mofified by set_label function.
   After it is finalized, it becomes "Destroyed potential".

- add get_delta_energy method to potential. (For monte carlo)
    get_delta_energy( S, tau, lwf, ispin, dS, ilatt, dtau, ilwf, dlwf)
    all params are optional.
    This is limited to one move per step (simple metropolis-hastings)
- 
