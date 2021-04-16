# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

ininame    = ['dx','dy','dz']
inivalue   = ['0.5','0.5','0.5']

ininame.append('latitude')
inivalue.append('-70.')

ininame.append('nx')
inivalue.append('63')
ininame.append('ny')
inivalue.append('63')
ininame.append('nz')
inivalue.append('128')

ininame.append('ocean')
inivalue.append('.T.')

ininame.append('surface_pressure')
inivalue.append('1013.25')

ininame.append('reference_state')
inivalue.append('single_value')

ininame.append('initializing_actions')
inivalue.append('set_constant_profiles')

ininame.append('alpha_surface')
inivalue.append('0')
ininame.append('slope_parallel_gradients')
inivalue.append('.T.')
ininame.append('slope_offset')
inivalue.append('.F.')
ininame.append('ambient_density_for_buoyancy')
inivalue.append('.T.')

ininame.append('dp_external')
inivalue.append('.T.')
ininame.append('dpdxy')
inivalue.append('0,0')
ininame.append('dp_level_b')
inivalue.append('0')

ininame.append('ug_surface')
inivalue.append('0.0')
ininame.append('vg_surface')
inivalue.append('0.0')

ininame.append('pt_surface')
inivalue.append('273.13')
ininame.append('pt_vertical_gradient')
inivalue.append('0.0')
ininame.append('pt_vertical_gradient_level')
inivalue.append('0.0')

ininame.append('sa_surface')
inivalue.append('34.5')
ininame.append('sa_vertical_gradient')
inivalue.append('0.0')
ininame.append('sa_vertical_gradient_level')
inivalue.append('0.0')

ininame.append('use_top_fluxes')
inivalue.append('.T.')
ininame.append('use_surface_fluxes')
inivalue.append('.F.')

ininame.append('constant_flux_layer')
inivalue.append('top')
ininame.append('most_method')
inivalue.append('mcphee')
ininame.append('drag_coeff')
inivalue.append('0.003')
ininame.append('most_xy_av')
inivalue.append('.F.')
ininame.append('zeta_min')
inivalue.append('5.0')
ininame.append('k_offset')
inivalue.append('3')

ininame.append('rayleigh_damping_factor')
inivalue.append('0.001')
ininame.append('rayleigh_damping_height')
inivalue.append('0')

ininame.append('top_heatflux')
inivalue.append('0.0')
ininame.append('surface_heatflux')
inivalue.append('9999999.9')

ininame.append('top_momentumflux_u')
inivalue.append('9999999.9')
ininame.append('top_momentumflux_v')
inivalue.append('9999999.9')

ininame.append('top_salinityflux')
inivalue.append('0.0')
ininame.append('bottom_salinityflux')
inivalue.append('9999999.9')

ininame.append('bc_uv_b')
inivalue.append('neumann')
ininame.append('bc_uv_t')
inivalue.append('dirichlet')

ininame.append('bc_pt_b')
inivalue.append('dirichlet')
ininame.append('bc_pt_t')
inivalue.append('neumann')

ininame.append('bc_sa_b')
inivalue.append('dirichlet')
ininame.append('bc_sa_t')
inivalue.append('neumann')

ininame.append('bc_s_b')
inivalue.append('dirichlet')
ininame.append('bc_s_t')
inivalue.append('neumann')

ininame.append('bc_p_b')
inivalue.append('neumann')
ininame.append('bc_p_t')
inivalue.append('neumann')

ininame.append('end_time')
inivalue.append('0.')

ininame.append('averaging_interval')
inivalue.append('0.')
ininame.append('averaging_interval_pr')
inivalue.append('3600.')

ininame.append('create_disturbances')
inivalue.append('.T.')

ininame.append('disturbance_energy_limit')
inivalue.append('1e-3')
ininame.append('disturbance_level_t')
inivalue.append('-10')
ininame.append('disturbance_level_b')
inivalue.append('-10')
ininame.append('dt_disturb')
inivalue.append('3600')

ininame.append('dt_run_control')
inivalue.append('0.')

ininame.append('dt_data_output')
inivalue.append('3600.')
ininame.append('dt_data_output_av')
inivalue.append('3600.')
ininame.append('dt_dopr')
inivalue.append('3600.')

ininame.append('section_xy')
inivalue.append('-9999')
ininame.append('section_yz')
inivalue.append('-9999')
ininame.append('section_xz')
inivalue.append('-9999')

ininame.append('termination_time_needed')
inivalue.append('3600')

ininame.append('data_output')
inivalue.append("'e', 'pt', 'sa', 'u', 'v', 'w', 'rho_ocean','melt*_xy','shf*_xy','sasws*_xy','u*_xy','ol*_xy','pt1*_xy','pt_io*_xy','sa1*_xy','sa_io*_xy'")
#inivalue.append("'e', 'pt', 'sa', 'u', 'v', 'w', 'rho_ocean'")

ininame.append('data_output_pr')
inivalue.append("'e', 'e*', '#pt', '#sa', 'p', 'hyp',"\
                "'#u', '#v', 'w', 'prho', 'rho_ocean',"\
                "'wu', 'w\"u\"', 'w*u*', 'wv', 'w\"v\"', 'w*v*',"\
                "'wpt', 'w\"pt\"', 'w*pt*', 'wsa', 'w\"sa\"', 'w*sa*',"\
                "'w*e*', 'u*2', 'v*2', 'w*2', 'pt*2','w*3', 'Sw',"\
                "'w*2pt*', 'w*pt*2', 'w*u*u*:dz', 'w*p*:dz', "\
                "'alpha_T', 'beta_S','km', 'kh', 'l'") 

