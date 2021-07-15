PREFIX_T=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_dpdy_021/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
PREFIX_SLOPE=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_alpha_surface_006/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
SUFFIX_OLD_T_AV=_cmp_dT0.2_dT0.4_43hr_tav13_42zmax_z_profile.png
SUFFIX_OLD_T=_cmp_dT0.2_dT0.4_40hr_42zmax_z_profile.png
SUFFIX_NEW_T_AV=_cmp_dT_43h_tav13h_z_profile.png
SUFFIX_NEW_T=_cmp_dT_40h_z_profile.png
SUFFIX_OLD_SLOPE=_cmp_slope0.10_slope0.50_40hr_42zmax_z_profile.png
SUFFIX_OLD_SLOPE_AV=_cmp_slope0.10_slope0.50_43hr_tav13.0_42zmax_z_profile.png
SUFFIX_NEW_SLOPE=_cmp_dslope_40h_z_profile.png
SUFFIX_NEW_SLOPE_AV=_cmp_dslope_43h_tav13h_z_profile.png

#cp $PREFIX_T/vel_var_ratio$SUFFIX_OLD_T ../Figures/vel_var_ratio$SUFFIX_NEW_T
#cp $PREFIX_T/k_all$SUFFIX_OLD_T_AV ../Figures/k$SUFFIX_NEW_T_AV
#cp $PREFIX_SLOPE/vel_var_ratio$SUFFIX_OLD_SLOPE ../Figures/vel_var_ratio$SUFFIX_NEW_SLOPE
#cp $PREFIX_SLOPE/k_all$SUFFIX_OLD_SLOPE_AV ../Figures/k$SUFFIX_NEW_SLOPE_AV
cp $PREFIX_T/km_eff_cmp_dT0.2_dT0.4_43hr_tav13_24zmax_z_profile.png ../Figures/km$SUFFIX_NEW_T_AV
cp $PREFIX_SLOPE/km_eff_cmp_slope0.10_slope0.50_43hr_tav13.0_24zmax_z_profile.png ../Figures/km$SUFFIX_NEW_SLOPE_AV
