PREFIX_T=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_dpdy_021/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
PREFIX_SLOPE=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_alpha_surface_006/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
SUFFIX_OLD_T=_cmp_dT0.2_dT0.4_tlim2-50_t.png
SUFFIX_NEW_T=_cmp_dT_t.png
SUFFIX_OLD_SLOPE=_cmp_slope0.10_slope0.50_tlim2-50_t.png
SUFFIX_NEW_SLOPE=_cmp_dslope_t.png

cp $PREFIX_T/melt$SUFFIX_OLD_T ../Figures/melt$SUFFIX_NEW_T
cp $PREFIX_T/e_res$SUFFIX_OLD_T ../Figures/tke$SUFFIX_NEW_T
cp $PREFIX_T/us$SUFFIX_OLD_T ../Figures/us$SUFFIX_NEW_T
cp $PREFIX_SLOPE/melt$SUFFIX_OLD_SLOPE ../Figures/melt$SUFFIX_NEW_SLOPE
cp $PREFIX_SLOPE/e_res$SUFFIX_OLD_SLOPE ../Figures/tke$SUFFIX_NEW_SLOPE
cp $PREFIX_SLOPE/us$SUFFIX_OLD_SLOPE ../Figures/us$SUFFIX_NEW_SLOPE
