#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Z2 invariant calculation for BiTeI at 0 and 5 GPa

import os
import z2pack

if not os.path.exists('./results_tz2_3'):
    os.makedirs('./results_tz2_3')

# creating the System object
# We will create one for each pressure.
BiTeI_0gpa = z2pack.fp.System(
    input_files=['input_tz2_3/tz2_3_0gpa.abi', 'input_tz2_3/wannier90.win' ],
    kpt_fct=z2pack.fp.kpoint.abinit,
    kpt_path='tz2_3_0gpa.abi',
    command='abinit tz2_3_0gpa.abi >& log',
    executable='/bin/bash'
)

BiTeI_5gpa = z2pack.fp.System(
    input_files=['input_tz2_3/tz2_3_5gpa.abi', 'input_tz2_3/wannier90.win' ],
    kpt_fct=z2pack.fp.kpoint.abinit,
    kpt_path='tz2_3_5gpa.abi',
    command='abinit tz2_3_5gpa.abi >& log',
    executable='/bin/bash'
)

# calculating the WCC
# Invariant at fixed k_z planes
result_0 = z2pack.surface.run(system=BiTeI_0gpa, surface=lambda s, t: [s / 2, t,0.0], num_lines=20, save_file = './results_tz2_3/BiTeI_0.msgpack', load=True)
result_1 = z2pack.surface.run(system=BiTeI_0gpa, surface=lambda s, t: [s / 2, t,0.5], num_lines=20, save_file = './results_tz2_3/BiTeI_1.msgpack', load=True)
result_2 = z2pack.surface.run(system=BiTeI_5gpa, surface=lambda s, t: [s / 2, t,0.0], num_lines=20, save_file = './results_tz2_3/BiTeI_2.msgpack', load=True)
result_3 = z2pack.surface.run(system=BiTeI_5gpa, surface=lambda s, t: [s / 2, t,0.5], num_lines=20, save_file = './results_tz2_3/BiTeI_3.msgpack', load=True)

print('for 0 GPa:')
print('z2 topological invariant at kz = 0.0: {0}'.format(z2pack.invariant.z2(result_0)))
print('z2 topological invariant at kz = 0.5: {0}'.format(z2pack.invariant.z2(result_1)))
print('for 5 GPa:')
print('z2 topological invariant at kz = 0.0: {0}'.format(z2pack.invariant.z2(result_2)))
print('z2 topological invariant at kz = 0.5: {0}'.format(z2pack.invariant.z2(result_3)))
