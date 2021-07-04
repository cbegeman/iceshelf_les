PREFIX_T=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_dpdy_021/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
PREFIX_SLOPE=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_alpha_surface_006/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
SUFFIX_OLD_T=_cmp_dT0.2_dT0.4_43hr_tav13_42zmax_z_profile.png
SUFFIX_NEW_T=_cmp_dT_43h_tav13h_z_profile.png
SUFFIX_OLD_SLOPE=_cmp_slope0.10_slope0.50_43hr_tav13.0_42zmax_z_profile.png
SUFFIX_NEW_SLOPE=_cmp_dslope_43h_tav13h_z_profile.png

cp $PREFIX_T/ptfar_diff$SUFFIX_OLD_T ../Figures/pt$SUFFIX_NEW_T
cp $PREFIX_T/safar_diff$SUFFIX_OLD_T ../Figures/sa$SUFFIX_NEW_T
cp $PREFIX_T/velocity$SUFFIX_OLD_T ../Figures/velocity$SUFFIX_NEW_T
cp $PREFIX_SLOPE/ptfar_diff$SUFFIX_OLD_SLOPE ../Figures/pt$SUFFIX_NEW_SLOPE
cp $PREFIX_SLOPE/safar_diff$SUFFIX_OLD_SLOPE ../Figures/sa$SUFFIX_NEW_SLOPE
cp $PREFIX_SLOPE/velocity$SUFFIX_OLD_SLOPE ../Figures/velocity$SUFFIX_NEW_SLOPE
