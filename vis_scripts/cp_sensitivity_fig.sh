PREFIX_T=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_dpdy_021/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
PREFIX_SLOPE=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_alpha_surface_006/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
SUFFIX_OLD_T=_cmp_dT0.2_dT0.4_tlim43-43_tav13.png
SUFFIX_NEW_T=_cmp_dT_43h_tav13h.png
SUFFIX_OLD_SLOPE=_cmp_slope0.10_slope0.50_tlim43-43_tav13.png
SUFFIX_NEW_SLOPE=_cmp_dslope_43h_tav13h.png

#cp $PREFIX_T/melt_thermal_driving_infty$SUFFIX_OLD_T ../Figures/melt_dT$SUFFIX_NEW_T
cp $PREFIX_T/melt_thermal_driving_infty_cmp_dT0.2_dT0.4_tlim43-43_tav13_tcycles.png ../Figures/melt_dT_cmp_dT_43h_tav13h_cycles.png
#cp $PREFIX_T/gamma_T_2m_thermal_driving_infty$SUFFIX_OLD_T ../Figures/gammaT_dT$SUFFIX_NEW_T
#cp $PREFIX_SLOPE/melt_sin_alpha$SUFFIX_OLD_SLOPE ../Figures/melt_slope$SUFFIX_NEW_SLOPE
#cp $PREFIX_SLOPE/gamma_T_2m_us$SUFFIX_OLD_SLOPE ../Figures/gammaT_us$SUFFIX_NEW_SLOPE
