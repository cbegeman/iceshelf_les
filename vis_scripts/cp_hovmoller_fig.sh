PREFIX_T=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_dpdy_021/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
PREFIX_SLOPE=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_dT_020/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
SUFFIX_OLD_T=_dT0.4_tlim2-50_42zmax_hovmoller.png
SUFFIX_NEW_T=_hovmoller.png
SUFFIX_OLD_SLOPE=$SUFFIX_OLD_T
SUFFIX_NEW_SLOPE=$SUFFIX_NEW_T

cp $PREFIX_T/e_res$SUFFIX_OLD_T ../Figures/e_dT01$SUFFIX_NEW_T
cp $PREFIX_SLOPE/e_res$SUFFIX_OLD_SLOPE ../Figures/e_dt06$SUFFIX_NEW_SLOPE
