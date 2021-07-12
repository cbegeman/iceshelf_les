PREFIX_T=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_dpdy_021/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
PREFIX_SLOPE=/turquoise/usr/projects/climate/cbegeman/palm/jobs/test_ocean_melt_batch_alpha_surface_006/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml
SUFFIX_Y=_xz_y64_zmax42_t40.png
SUFFIX_Z=_xy_z1_zmax1_t40.png

cp $PREFIX_T/pt$SUFFIX_Y ../Figures/pt_slope1$SUFFIX_Y
cp $PREFIX_T/sa$SUFFIX_Y ../Figures/sa_slope1$SUFFIX_Y
cp $PREFIX_T/v$SUFFIX_Z ../Figures/v_slope1$SUFFIX_Z
cp $PREFIX_SLOPE/pt$SUFFIX_Y ../Figures/pt_slope01$SUFFIX_Y
cp $PREFIX_SLOPE/sa$SUFFIX_Y ../Figures/sa_slope01$SUFFIX_Y
cp $PREFIX_SLOPE/v$SUFFIX_Z ../Figures/v_slope01$SUFFIX_Z
