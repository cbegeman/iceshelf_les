PREFIX_T=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_dpdy_021/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
PREFIX_SLOPE=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_alpha_surface_006/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
SUFFIX_OLD_T=_cmp_dT0.2_dT0.4_43hr_tav13_42zmax_z_profile.png
SUFFIX_OLD_T_SCALE=_cmp_dT0.2_dT0.4_43hr_tav13_scale_z_profile.png
SUFFIX_NEW_T=_cmp_dT_43h_tav13h_z_profile.png
SUFFIX_NEW_T_SCALE=_cmp_dT_43h_tav13h_dTscale_z_profile.png
SUFFIX_OLD_SLOPE=_cmp_slope0.10_slope0.50_43hr_tav13.0_42zmax_z_profile.png
SUFFIX_OLD_SLOPE_SCALE=_cmp_slope0.10_slope0.50_43hr_tav13.0_scale_z_profile.png
SUFFIX_NEW_SLOPE=_cmp_dslope_43h_tav13h_z_profile.png
SUFFIX_NEW_SLOPE_SCALE=_cmp_dslope_43h_tav13h_slopescale_z_profile.png

cp $PREFIX_T/wpt$SUFFIX_OLD_T ../Figures/heatflux$SUFFIX_NEW_T
cp $PREFIX_T/wpt$SUFFIX_OLD_T_SCALE ../Figures/heatflux$SUFFIX_NEW_T_SCALE
cp $PREFIX_T/wsa$SUFFIX_OLD_T ../Figures/wsa$SUFFIX_NEW_T
cp $PREFIX_T/momflux_z$SUFFIX_OLD_T ../Figures/momflux$SUFFIX_NEW_T
cp $PREFIX_SLOPE/wpt$SUFFIX_OLD_SLOPE ../Figures/heatflux$SUFFIX_NEW_SLOPE
cp $PREFIX_SLOPE/wpt$SUFFIX_OLD_SLOPE_SCALE ../Figures/heatflux$SUFFIX_NEW_SLOPE_SCALE
cp $PREFIX_SLOPE/wsa$SUFFIX_OLD_SLOPE ../Figures/wsa$SUFFIX_NEW_SLOPE
cp $PREFIX_SLOPE/momflux_z$SUFFIX_OLD_SLOPE ../Figures/momflux$SUFFIX_NEW_SLOPE
