#PREFIX_T=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_dpdy_021/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
PREFIX_T=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_dT_020/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
PREFIX_SLOPE=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_alpha_surface_006/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
SUFFIX_OLD_T=_dT0.4_tlim2-50_42zmax_contour_diss_hovmoller.png
SUFFIX_NEW_T=_hovmoller.png
SUFFIX_OLD_SLOPE=_slope0.50_tlim2-50_42zmax_contour_diss_hovmoller.png
SUFFIX_NEW_SLOPE=_hovmoller.png

cp $PREFIX_T/e_res$SUFFIX_OLD_T ../Figures/dT4/e$SUFFIX_NEW_T
cp $PREFIX_SLOPE/e_res$SUFFIX_OLD_SLOPE ../Figures/dslope2/e$SUFFIX_NEW_SLOPE
