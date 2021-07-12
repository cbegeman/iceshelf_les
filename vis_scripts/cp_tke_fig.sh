PREFIX_T=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_dpdy_021/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
PREFIX_SLOPE=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_alpha_surface_006/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
SUFFIX_OLD_T=_cmp_dT0.2_dT0.4_43hr_tav13_42zmax_z_profile.png
SUFFIX_NEW_T=_cmp_dT_43h_tav13h_z_profile.png
SUFFIX_OLD_SLOPE=_cmp_slope0.10_slope0.50_43hr_tav13.0_42zmax_z_profile.png
SUFFIX_NEW_SLOPE=_cmp_dslope_43h_tav13h_z_profile.png

cp $PREFIX_T/e_res$SUFFIX_OLD_T ../Figures/eres$SUFFIX_NEW_T
cp $PREFIX_T/Fshear$SUFFIX_OLD_T ../Figures/Fshear$SUFFIX_NEW_T
cp $PREFIX_T/Fbuoy_uw$SUFFIX_OLD_T ../Figures/Fbuoy$SUFFIX_NEW_T
cp $PREFIX_T/Ftrans$SUFFIX_OLD_T ../Figures/Ftrans$SUFFIX_NEW_T

cp $PREFIX_SLOPE/e_res$SUFFIX_OLD_SLOPE ../Figures/eres$SUFFIX_NEW_SLOPE
cp $PREFIX_SLOPE/Fshear$SUFFIX_OLD_SLOPE ../Figures/Fshear$SUFFIX_NEW_SLOPE
cp $PREFIX_SLOPE/Fbuoy_uw$SUFFIX_OLD_SLOPE ../Figures/Fbuoy$SUFFIX_NEW_SLOPE
cp $PREFIX_SLOPE/Ftrans$SUFFIX_OLD_SLOPE ../Figures/Ftrans$SUFFIX_NEW_SLOPE
