/*
 * Part P5: Declare the function to evaluate the scalar wave RHSs
 */
void rhs_eval(rfm_struct *restrict rfmstruct,const paramstruct *restrict params,
                 const REAL *restrict in_gfs, REAL *restrict rhs_gfs) {
#include "./set_Cparameters-SIMD.h"

  #pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++) {
    #include "rfm_files/rfm_struct__SIMD_outer_read2.h"
    for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++) {
      #include "rfm_files/rfm_struct__SIMD_outer_read1.h"
      for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0 += SIMD_width) {
        #include "rfm_files/rfm_struct__SIMD_inner_read0.h"
        {
           /*
            * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
            */
           /*
            *  Original SymPy expressions:
            *  "[const REAL_SIMD_ARRAY bb_dD0 = invdx0*(-bb_i0m1_i1_i2/2 + bb_i0p1_i1_i2/2),
            *    const REAL_SIMD_ARRAY bb_dD1 = invdx1*(-bb_i0_i1m1_i2/2 + bb_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY bb_dD2 = invdx2*(-bb_i0_i1_i2m1/2 + bb_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY bb_dKOD0 = invdx0*(-3*bb/8 + bb_i0m1_i1_i2/4 - bb_i0m2_i1_i2/16 + bb_i0p1_i1_i2/4 - bb_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY bb_dKOD1 = invdx1*(-3*bb/8 + bb_i0_i1m1_i2/4 - bb_i0_i1m2_i2/16 + bb_i0_i1p1_i2/4 - bb_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY bb_dKOD2 = invdx2*(-3*bb/8 + bb_i0_i1_i2m1/4 - bb_i0_i1_i2m2/16 + bb_i0_i1_i2p1/4 - bb_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY bp_dD0 = invdx0*(-bp_i0m1_i1_i2/2 + bp_i0p1_i1_i2/2),
            *    const REAL_SIMD_ARRAY bp_dD1 = invdx1*(-bp_i0_i1m1_i2/2 + bp_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY bp_dD2 = invdx2*(-bp_i0_i1_i2m1/2 + bp_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY bp_dKOD0 = invdx0*(-3*bp/8 + bp_i0m1_i1_i2/4 - bp_i0m2_i1_i2/16 + bp_i0p1_i1_i2/4 - bp_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY bp_dKOD1 = invdx1*(-3*bp/8 + bp_i0_i1m1_i2/4 - bp_i0_i1m2_i2/16 + bp_i0_i1p1_i2/4 - bp_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY bp_dKOD2 = invdx2*(-3*bp/8 + bp_i0_i1_i2m1/4 - bp_i0_i1_i2m2/16 + bp_i0_i1_i2p1/4 - bp_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY fgauge_dD0 = invdx0*(-fgauge_i0m1_i1_i2/2 + fgauge_i0p1_i1_i2/2),
            *    const REAL_SIMD_ARRAY fgauge_dD1 = invdx1*(-fgauge_i0_i1m1_i2/2 + fgauge_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY fgauge_dD2 = invdx2*(-fgauge_i0_i1_i2m1/2 + fgauge_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY fgauge_dKOD0 = invdx0*(-3*fgauge/8 + fgauge_i0m1_i1_i2/4 - fgauge_i0m2_i1_i2/16 + fgauge_i0p1_i1_i2/4 - fgauge_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY fgauge_dKOD1 = invdx1*(-3*fgauge/8 + fgauge_i0_i1m1_i2/4 - fgauge_i0_i1m2_i2/16 + fgauge_i0_i1p1_i2/4 - fgauge_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY fgauge_dKOD2 = invdx2*(-3*fgauge/8 + fgauge_i0_i1_i2m1/4 - fgauge_i0_i1_i2m2/16 + fgauge_i0_i1_i2p1/4 - fgauge_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY fgaugep_dD0 = invdx0*(-fgaugep_i0m1_i1_i2/2 + fgaugep_i0p1_i1_i2/2),
            *    const REAL_SIMD_ARRAY fgaugep_dD1 = invdx1*(-fgaugep_i0_i1m1_i2/2 + fgaugep_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY fgaugep_dD2 = invdx2*(-fgaugep_i0_i1_i2m1/2 + fgaugep_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY fgaugep_dKOD0 = invdx0*(-3*fgaugep/8 + fgaugep_i0m1_i1_i2/4 - fgaugep_i0m2_i1_i2/16 + fgaugep_i0p1_i1_i2/4 - fgaugep_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY fgaugep_dKOD1 = invdx1*(-3*fgaugep/8 + fgaugep_i0_i1m1_i2/4 - fgaugep_i0_i1m2_i2/16 + fgaugep_i0_i1p1_i2/4 - fgaugep_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY fgaugep_dKOD2 = invdx2*(-3*fgaugep/8 + fgaugep_i0_i1_i2m1/4 - fgaugep_i0_i1_i2m2/16 + fgaugep_i0_i1_i2p1/4 - fgaugep_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY gg_dD0 = invdx0*(-gg_i0m1_i1_i2/2 + gg_i0p1_i1_i2/2),
            *    const REAL_SIMD_ARRAY gg_dD1 = invdx1*(-gg_i0_i1m1_i2/2 + gg_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY gg_dD2 = invdx2*(-gg_i0_i1_i2m1/2 + gg_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY gg_dKOD0 = invdx0*(-3*gg/8 + gg_i0m1_i1_i2/4 - gg_i0m2_i1_i2/16 + gg_i0p1_i1_i2/4 - gg_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY gg_dKOD1 = invdx1*(-3*gg/8 + gg_i0_i1m1_i2/4 - gg_i0_i1m2_i2/16 + gg_i0_i1p1_i2/4 - gg_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY gg_dKOD2 = invdx2*(-3*gg/8 + gg_i0_i1_i2m1/4 - gg_i0_i1_i2m2/16 + gg_i0_i1_i2p1/4 - gg_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY gp_dD0 = invdx0*(-gp_i0m1_i1_i2/2 + gp_i0p1_i1_i2/2),
            *    const REAL_SIMD_ARRAY gp_dD1 = invdx1*(-gp_i0_i1m1_i2/2 + gp_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY gp_dD2 = invdx2*(-gp_i0_i1_i2m1/2 + gp_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY gp_dKOD0 = invdx0*(-3*gp/8 + gp_i0m1_i1_i2/4 - gp_i0m2_i1_i2/16 + gp_i0p1_i1_i2/4 - gp_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY gp_dKOD1 = invdx1*(-3*gp/8 + gp_i0_i1m1_i2/4 - gp_i0_i1m2_i2/16 + gp_i0_i1p1_i2/4 - gp_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY gp_dKOD2 = invdx2*(-3*gp/8 + gp_i0_i1_i2m1/4 - gp_i0_i1_i2m2/16 + gp_i0_i1_i2p1/4 - gp_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY psibD_dD00 = invdx0*(-psibD0_i0m1_i1_i2/2 + psibD0_i0p1_i1_i2/2),
            *    const REAL_SIMD_ARRAY psibD_dD01 = invdx1*(-psibD0_i0_i1m1_i2/2 + psibD0_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY psibD_dD02 = invdx2*(-psibD0_i0_i1_i2m1/2 + psibD0_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY psibD_dD11 = invdx1*(-psibD1_i0_i1m1_i2/2 + psibD1_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY psibD_dD22 = invdx2*(-psibD2_i0_i1_i2m1/2 + psibD2_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY psibD_dKOD00 = invdx0*(-3*psibD0/8 + psibD0_i0m1_i1_i2/4 - psibD0_i0m2_i1_i2/16 + psibD0_i0p1_i1_i2/4 - psibD0_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY psibD_dKOD01 = invdx1*(-3*psibD0/8 + psibD0_i0_i1m1_i2/4 - psibD0_i0_i1m2_i2/16 + psibD0_i0_i1p1_i2/4 - psibD0_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY psibD_dKOD02 = invdx2*(-3*psibD0/8 + psibD0_i0_i1_i2m1/4 - psibD0_i0_i1_i2m2/16 + psibD0_i0_i1_i2p1/4 - psibD0_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY psibD_dKOD10 = invdx0*(-3*psibD1/8 + psibD1_i0m1_i1_i2/4 - psibD1_i0m2_i1_i2/16 + psibD1_i0p1_i1_i2/4 - psibD1_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY psibD_dKOD11 = invdx1*(-3*psibD1/8 + psibD1_i0_i1m1_i2/4 - psibD1_i0_i1m2_i2/16 + psibD1_i0_i1p1_i2/4 - psibD1_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY psibD_dKOD12 = invdx2*(-3*psibD1/8 + psibD1_i0_i1_i2m1/4 - psibD1_i0_i1_i2m2/16 + psibD1_i0_i1_i2p1/4 - psibD1_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY psibD_dKOD20 = invdx0*(-3*psibD2/8 + psibD2_i0m1_i1_i2/4 - psibD2_i0m2_i1_i2/16 + psibD2_i0p1_i1_i2/4 - psibD2_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY psibD_dKOD21 = invdx1*(-3*psibD2/8 + psibD2_i0_i1m1_i2/4 - psibD2_i0_i1m2_i2/16 + psibD2_i0_i1p1_i2/4 - psibD2_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY psibD_dKOD22 = invdx2*(-3*psibD2/8 + psibD2_i0_i1_i2m1/4 - psibD2_i0_i1_i2m2/16 + psibD2_i0_i1_i2p1/4 - psibD2_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY psifD_dD00 = invdx0*(-psifD0_i0m1_i1_i2/2 + psifD0_i0p1_i1_i2/2),
            *    const REAL_SIMD_ARRAY psifD_dD01 = invdx1*(-psifD0_i0_i1m1_i2/2 + psifD0_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY psifD_dD02 = invdx2*(-psifD0_i0_i1_i2m1/2 + psifD0_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY psifD_dD11 = invdx1*(-psifD1_i0_i1m1_i2/2 + psifD1_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY psifD_dD22 = invdx2*(-psifD2_i0_i1_i2m1/2 + psifD2_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY psifD_dKOD00 = invdx0*(-3*psifD0/8 + psifD0_i0m1_i1_i2/4 - psifD0_i0m2_i1_i2/16 + psifD0_i0p1_i1_i2/4 - psifD0_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY psifD_dKOD01 = invdx1*(-3*psifD0/8 + psifD0_i0_i1m1_i2/4 - psifD0_i0_i1m2_i2/16 + psifD0_i0_i1p1_i2/4 - psifD0_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY psifD_dKOD02 = invdx2*(-3*psifD0/8 + psifD0_i0_i1_i2m1/4 - psifD0_i0_i1_i2m2/16 + psifD0_i0_i1_i2p1/4 - psifD0_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY psifD_dKOD10 = invdx0*(-3*psifD1/8 + psifD1_i0m1_i1_i2/4 - psifD1_i0m2_i1_i2/16 + psifD1_i0p1_i1_i2/4 - psifD1_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY psifD_dKOD11 = invdx1*(-3*psifD1/8 + psifD1_i0_i1m1_i2/4 - psifD1_i0_i1m2_i2/16 + psifD1_i0_i1p1_i2/4 - psifD1_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY psifD_dKOD12 = invdx2*(-3*psifD1/8 + psifD1_i0_i1_i2m1/4 - psifD1_i0_i1_i2m2/16 + psifD1_i0_i1_i2p1/4 - psifD1_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY psifD_dKOD20 = invdx0*(-3*psifD2/8 + psifD2_i0m1_i1_i2/4 - psifD2_i0m2_i1_i2/16 + psifD2_i0p1_i1_i2/4 - psifD2_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY psifD_dKOD21 = invdx1*(-3*psifD2/8 + psifD2_i0_i1m1_i2/4 - psifD2_i0_i1m2_i2/16 + psifD2_i0_i1p1_i2/4 - psifD2_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY psifD_dKOD22 = invdx2*(-3*psifD2/8 + psifD2_i0_i1_i2m1/4 - psifD2_i0_i1_i2m2/16 + psifD2_i0_i1_i2p1/4 - psifD2_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY psigD_dD00 = invdx0*(-psigD0_i0m1_i1_i2/2 + psigD0_i0p1_i1_i2/2),
            *    const REAL_SIMD_ARRAY psigD_dD01 = invdx1*(-psigD0_i0_i1m1_i2/2 + psigD0_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY psigD_dD02 = invdx2*(-psigD0_i0_i1_i2m1/2 + psigD0_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY psigD_dD11 = invdx1*(-psigD1_i0_i1m1_i2/2 + psigD1_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY psigD_dD22 = invdx2*(-psigD2_i0_i1_i2m1/2 + psigD2_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY psigD_dKOD00 = invdx0*(-3*psigD0/8 + psigD0_i0m1_i1_i2/4 - psigD0_i0m2_i1_i2/16 + psigD0_i0p1_i1_i2/4 - psigD0_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY psigD_dKOD01 = invdx1*(-3*psigD0/8 + psigD0_i0_i1m1_i2/4 - psigD0_i0_i1m2_i2/16 + psigD0_i0_i1p1_i2/4 - psigD0_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY psigD_dKOD02 = invdx2*(-3*psigD0/8 + psigD0_i0_i1_i2m1/4 - psigD0_i0_i1_i2m2/16 + psigD0_i0_i1_i2p1/4 - psigD0_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY psigD_dKOD10 = invdx0*(-3*psigD1/8 + psigD1_i0m1_i1_i2/4 - psigD1_i0m2_i1_i2/16 + psigD1_i0p1_i1_i2/4 - psigD1_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY psigD_dKOD11 = invdx1*(-3*psigD1/8 + psigD1_i0_i1m1_i2/4 - psigD1_i0_i1m2_i2/16 + psigD1_i0_i1p1_i2/4 - psigD1_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY psigD_dKOD12 = invdx2*(-3*psigD1/8 + psigD1_i0_i1_i2m1/4 - psigD1_i0_i1_i2m2/16 + psigD1_i0_i1_i2p1/4 - psigD1_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY psigD_dKOD20 = invdx0*(-3*psigD2/8 + psigD2_i0m1_i1_i2/4 - psigD2_i0m2_i1_i2/16 + psigD2_i0p1_i1_i2/4 - psigD2_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY psigD_dKOD21 = invdx1*(-3*psigD2/8 + psigD2_i0_i1m1_i2/4 - psigD2_i0_i1m2_i2/16 + psigD2_i0_i1p1_i2/4 - psigD2_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY psigD_dKOD22 = invdx2*(-3*psigD2/8 + psigD2_i0_i1_i2m1/4 - psigD2_i0_i1_i2m2/16 + psigD2_i0_i1_i2p1/4 - psigD2_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY psiuD_dD00 = invdx0*(-psiuD0_i0m1_i1_i2/2 + psiuD0_i0p1_i1_i2/2),
            *    const REAL_SIMD_ARRAY psiuD_dD01 = invdx1*(-psiuD0_i0_i1m1_i2/2 + psiuD0_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY psiuD_dD02 = invdx2*(-psiuD0_i0_i1_i2m1/2 + psiuD0_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY psiuD_dD11 = invdx1*(-psiuD1_i0_i1m1_i2/2 + psiuD1_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY psiuD_dD22 = invdx2*(-psiuD2_i0_i1_i2m1/2 + psiuD2_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY psiuD_dKOD00 = invdx0*(-3*psiuD0/8 + psiuD0_i0m1_i1_i2/4 - psiuD0_i0m2_i1_i2/16 + psiuD0_i0p1_i1_i2/4 - psiuD0_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY psiuD_dKOD01 = invdx1*(-3*psiuD0/8 + psiuD0_i0_i1m1_i2/4 - psiuD0_i0_i1m2_i2/16 + psiuD0_i0_i1p1_i2/4 - psiuD0_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY psiuD_dKOD02 = invdx2*(-3*psiuD0/8 + psiuD0_i0_i1_i2m1/4 - psiuD0_i0_i1_i2m2/16 + psiuD0_i0_i1_i2p1/4 - psiuD0_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY psiuD_dKOD10 = invdx0*(-3*psiuD1/8 + psiuD1_i0m1_i1_i2/4 - psiuD1_i0m2_i1_i2/16 + psiuD1_i0p1_i1_i2/4 - psiuD1_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY psiuD_dKOD11 = invdx1*(-3*psiuD1/8 + psiuD1_i0_i1m1_i2/4 - psiuD1_i0_i1m2_i2/16 + psiuD1_i0_i1p1_i2/4 - psiuD1_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY psiuD_dKOD12 = invdx2*(-3*psiuD1/8 + psiuD1_i0_i1_i2m1/4 - psiuD1_i0_i1_i2m2/16 + psiuD1_i0_i1_i2p1/4 - psiuD1_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY psiuD_dKOD20 = invdx0*(-3*psiuD2/8 + psiuD2_i0m1_i1_i2/4 - psiuD2_i0m2_i1_i2/16 + psiuD2_i0p1_i1_i2/4 - psiuD2_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY psiuD_dKOD21 = invdx1*(-3*psiuD2/8 + psiuD2_i0_i1m1_i2/4 - psiuD2_i0_i1m2_i2/16 + psiuD2_i0_i1p1_i2/4 - psiuD2_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY psiuD_dKOD22 = invdx2*(-3*psiuD2/8 + psiuD2_i0_i1_i2m1/4 - psiuD2_i0_i1_i2m2/16 + psiuD2_i0_i1_i2p1/4 - psiuD2_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY up_dD0 = invdx0*(-up_i0m1_i1_i2/2 + up_i0p1_i1_i2/2),
            *    const REAL_SIMD_ARRAY up_dD1 = invdx1*(-up_i0_i1m1_i2/2 + up_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY up_dD2 = invdx2*(-up_i0_i1_i2m1/2 + up_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY up_dKOD0 = invdx0*(-3*up/8 + up_i0m1_i1_i2/4 - up_i0m2_i1_i2/16 + up_i0p1_i1_i2/4 - up_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY up_dKOD1 = invdx1*(-3*up/8 + up_i0_i1m1_i2/4 - up_i0_i1m2_i2/16 + up_i0_i1p1_i2/4 - up_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY up_dKOD2 = invdx2*(-3*up/8 + up_i0_i1_i2m1/4 - up_i0_i1_i2m2/16 + up_i0_i1_i2p1/4 - up_i0_i1_i2p2/16),
            *    const REAL_SIMD_ARRAY uu_dD0 = invdx0*(-uu_i0m1_i1_i2/2 + uu_i0p1_i1_i2/2),
            *    const REAL_SIMD_ARRAY uu_dD1 = invdx1*(-uu_i0_i1m1_i2/2 + uu_i0_i1p1_i2/2),
            *    const REAL_SIMD_ARRAY uu_dD2 = invdx2*(-uu_i0_i1_i2m1/2 + uu_i0_i1_i2p1/2),
            *    const REAL_SIMD_ARRAY uu_dKOD0 = invdx0*(-3*uu/8 + uu_i0m1_i1_i2/4 - uu_i0m2_i1_i2/16 + uu_i0p1_i1_i2/4 - uu_i0p2_i1_i2/16),
            *    const REAL_SIMD_ARRAY uu_dKOD1 = invdx1*(-3*uu/8 + uu_i0_i1m1_i2/4 - uu_i0_i1m2_i2/16 + uu_i0_i1p1_i2/4 - uu_i0_i1p2_i2/16),
            *    const REAL_SIMD_ARRAY uu_dKOD2 = invdx2*(-3*uu/8 + uu_i0_i1_i2m1/4 - uu_i0_i1_i2m2/16 + uu_i0_i1_i2p1/4 - uu_i0_i1_i2p2/16)]"
            */
           const REAL_SIMD_ARRAY gg_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(GGGF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY gg_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(GGGF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY gg_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(GGGF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY gg_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(GGGF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY gg_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(GGGF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY gg_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(GGGF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY gg = ReadSIMD(&in_gfs[IDX4S(GGGF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY gg_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(GGGF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY gg_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(GGGF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY gg_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(GGGF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY gg_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(GGGF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY gg_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(GGGF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY gg_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(GGGF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY gp_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(GPGF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY gp_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(GPGF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY gp_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(GPGF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY gp_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(GPGF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY gp_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(GPGF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY gp_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(GPGF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY gp = ReadSIMD(&in_gfs[IDX4S(GPGF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY gp_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(GPGF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY gp_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(GPGF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY gp_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(GPGF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY gp_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(GPGF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY gp_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(GPGF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY gp_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(GPGF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY bb_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(BBGF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY bb_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(BBGF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY bb_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(BBGF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY bb_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(BBGF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY bb_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(BBGF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY bb_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(BBGF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY bb = ReadSIMD(&in_gfs[IDX4S(BBGF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY bb_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(BBGF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY bb_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(BBGF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY bb_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(BBGF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY bb_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(BBGF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY bb_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(BBGF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY bb_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(BBGF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY bp_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(BPGF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY bp_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(BPGF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY bp_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(BPGF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY bp_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(BPGF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY bp_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(BPGF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY bp_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(BPGF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY bp = ReadSIMD(&in_gfs[IDX4S(BPGF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY bp_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(BPGF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY bp_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(BPGF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY bp_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(BPGF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY bp_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(BPGF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY bp_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(BPGF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY bp_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(BPGF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY uu_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY uu_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY uu_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY uu_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY uu_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY uu_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY uu = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY uu_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY uu_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY uu_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY uu_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY uu_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY uu_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY up_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(UPGF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY up_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(UPGF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY up_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(UPGF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY up_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(UPGF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY up_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UPGF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY up_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UPGF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY up = ReadSIMD(&in_gfs[IDX4S(UPGF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY up_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UPGF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY up_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UPGF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY up_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(UPGF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY up_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(UPGF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY up_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(UPGF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY up_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(UPGF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY fgauge_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEGF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY fgauge_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(FGAUGEGF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY fgauge_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEGF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY fgauge_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEGF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY fgauge_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEGF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY fgauge_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEGF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY fgauge = ReadSIMD(&in_gfs[IDX4S(FGAUGEGF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY fgauge_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEGF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY fgauge_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEGF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY fgauge_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEGF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY fgauge_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEGF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY fgauge_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(FGAUGEGF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY fgauge_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEGF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY fgaugep_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEPGF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY fgaugep_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(FGAUGEPGF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY fgaugep_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEPGF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY fgaugep_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEPGF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY fgaugep_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEPGF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY fgaugep_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEPGF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY fgaugep = ReadSIMD(&in_gfs[IDX4S(FGAUGEPGF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY fgaugep_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEPGF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY fgaugep_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEPGF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY fgaugep_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEPGF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY fgaugep_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEPGF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY fgaugep_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(FGAUGEPGF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY fgaugep_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(FGAUGEPGF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY chi = ReadSIMD(&in_gfs[IDX4S(CHIGF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY chip = ReadSIMD(&in_gfs[IDX4S(CHIPGF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY chipp = ReadSIMD(&in_gfs[IDX4S(CHIPPGF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY psigD0_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(PSIGD0GF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY psigD0_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(PSIGD0GF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY psigD0_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD0GF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY psigD0_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD0GF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY psigD0_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD0GF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY psigD0_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD0GF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY psigD0 = ReadSIMD(&in_gfs[IDX4S(PSIGD0GF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY psigD0_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD0GF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY psigD0_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD0GF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY psigD0_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD0GF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY psigD0_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD0GF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY psigD0_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(PSIGD0GF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY psigD0_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(PSIGD0GF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY psigD1_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(PSIGD1GF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY psigD1_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(PSIGD1GF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY psigD1_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD1GF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY psigD1_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD1GF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY psigD1_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD1GF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY psigD1_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD1GF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY psigD1 = ReadSIMD(&in_gfs[IDX4S(PSIGD1GF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY psigD1_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD1GF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY psigD1_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD1GF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY psigD1_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD1GF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY psigD1_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD1GF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY psigD1_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(PSIGD1GF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY psigD1_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(PSIGD1GF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY psigD2_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(PSIGD2GF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY psigD2_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(PSIGD2GF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY psigD2_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD2GF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY psigD2_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD2GF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY psigD2_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD2GF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY psigD2_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD2GF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY psigD2 = ReadSIMD(&in_gfs[IDX4S(PSIGD2GF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY psigD2_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD2GF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY psigD2_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD2GF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY psigD2_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD2GF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY psigD2_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIGD2GF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY psigD2_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(PSIGD2GF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY psigD2_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(PSIGD2GF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY psibD0_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(PSIBD0GF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY psibD0_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(PSIBD0GF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY psibD0_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD0GF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY psibD0_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD0GF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY psibD0_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD0GF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY psibD0_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD0GF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY psibD0 = ReadSIMD(&in_gfs[IDX4S(PSIBD0GF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY psibD0_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD0GF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY psibD0_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD0GF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY psibD0_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD0GF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY psibD0_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD0GF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY psibD0_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(PSIBD0GF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY psibD0_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(PSIBD0GF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY psibD1_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(PSIBD1GF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY psibD1_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(PSIBD1GF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY psibD1_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD1GF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY psibD1_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD1GF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY psibD1_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD1GF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY psibD1_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD1GF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY psibD1 = ReadSIMD(&in_gfs[IDX4S(PSIBD1GF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY psibD1_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD1GF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY psibD1_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD1GF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY psibD1_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD1GF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY psibD1_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD1GF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY psibD1_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(PSIBD1GF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY psibD1_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(PSIBD1GF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY psibD2_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(PSIBD2GF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY psibD2_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(PSIBD2GF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY psibD2_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD2GF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY psibD2_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD2GF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY psibD2_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD2GF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY psibD2_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD2GF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY psibD2 = ReadSIMD(&in_gfs[IDX4S(PSIBD2GF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY psibD2_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD2GF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY psibD2_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD2GF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY psibD2_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD2GF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY psibD2_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIBD2GF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY psibD2_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(PSIBD2GF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY psibD2_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(PSIBD2GF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY psiuD0_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(PSIUD0GF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY psiuD0_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(PSIUD0GF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY psiuD0_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD0GF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY psiuD0_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD0GF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY psiuD0_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD0GF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY psiuD0_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD0GF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY psiuD0 = ReadSIMD(&in_gfs[IDX4S(PSIUD0GF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY psiuD0_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD0GF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY psiuD0_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD0GF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY psiuD0_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD0GF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY psiuD0_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD0GF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY psiuD0_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(PSIUD0GF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY psiuD0_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(PSIUD0GF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY psiuD1_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(PSIUD1GF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY psiuD1_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(PSIUD1GF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY psiuD1_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD1GF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY psiuD1_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD1GF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY psiuD1_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD1GF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY psiuD1_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD1GF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY psiuD1 = ReadSIMD(&in_gfs[IDX4S(PSIUD1GF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY psiuD1_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD1GF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY psiuD1_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD1GF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY psiuD1_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD1GF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY psiuD1_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD1GF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY psiuD1_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(PSIUD1GF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY psiuD1_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(PSIUD1GF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY psiuD2_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(PSIUD2GF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY psiuD2_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(PSIUD2GF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY psiuD2_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD2GF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY psiuD2_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD2GF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY psiuD2_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD2GF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY psiuD2_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD2GF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY psiuD2 = ReadSIMD(&in_gfs[IDX4S(PSIUD2GF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY psiuD2_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD2GF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY psiuD2_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD2GF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY psiuD2_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD2GF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY psiuD2_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIUD2GF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY psiuD2_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(PSIUD2GF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY psiuD2_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(PSIUD2GF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY psifD0_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(PSIFD0GF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY psifD0_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(PSIFD0GF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY psifD0_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD0GF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY psifD0_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD0GF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY psifD0_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD0GF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY psifD0_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD0GF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY psifD0 = ReadSIMD(&in_gfs[IDX4S(PSIFD0GF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY psifD0_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD0GF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY psifD0_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD0GF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY psifD0_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD0GF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY psifD0_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD0GF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY psifD0_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(PSIFD0GF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY psifD0_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(PSIFD0GF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY psifD1_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(PSIFD1GF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY psifD1_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(PSIFD1GF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY psifD1_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD1GF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY psifD1_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD1GF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY psifD1_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD1GF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY psifD1_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD1GF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY psifD1 = ReadSIMD(&in_gfs[IDX4S(PSIFD1GF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY psifD1_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD1GF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY psifD1_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD1GF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY psifD1_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD1GF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY psifD1_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD1GF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY psifD1_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(PSIFD1GF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY psifD1_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(PSIFD1GF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY psifD2_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(PSIFD2GF, i0,i1,i2-2)]);
           const REAL_SIMD_ARRAY psifD2_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(PSIFD2GF, i0,i1,i2-1)]);
           const REAL_SIMD_ARRAY psifD2_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD2GF, i0,i1-2,i2)]);
           const REAL_SIMD_ARRAY psifD2_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD2GF, i0,i1-1,i2)]);
           const REAL_SIMD_ARRAY psifD2_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD2GF, i0-2,i1,i2)]);
           const REAL_SIMD_ARRAY psifD2_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD2GF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY psifD2 = ReadSIMD(&in_gfs[IDX4S(PSIFD2GF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY psifD2_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD2GF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY psifD2_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD2GF, i0+2,i1,i2)]);
           const REAL_SIMD_ARRAY psifD2_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD2GF, i0,i1+1,i2)]);
           const REAL_SIMD_ARRAY psifD2_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(PSIFD2GF, i0,i1+2,i2)]);
           const REAL_SIMD_ARRAY psifD2_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(PSIFD2GF, i0,i1,i2+1)]);
           const REAL_SIMD_ARRAY psifD2_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(PSIFD2GF, i0,i1,i2+2)]);
           const REAL_SIMD_ARRAY rcoord = ReadSIMD(&in_gfs[IDX4S(RCOORDGF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY Radius = ReadSIMD(&in_gfs[IDX4S(RADIUSGF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY Rprime = ReadSIMD(&in_gfs[IDX4S(RPRIMEGF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY sinth = ReadSIMD(&in_gfs[IDX4S(SINTHGF, i0,i1,i2)]);
           const REAL_SIMD_ARRAY costh = ReadSIMD(&in_gfs[IDX4S(COSTHGF, i0,i1,i2)]);
           const double tmpFDPart1_NegativeOne_ = -1.0;
           const REAL_SIMD_ARRAY FDPart1_NegativeOne_ = ConstSIMD(tmpFDPart1_NegativeOne_);

           const double tmpFDPart1_Rational_1_16 = 1.0/16.0;
           const REAL_SIMD_ARRAY FDPart1_Rational_1_16 = ConstSIMD(tmpFDPart1_Rational_1_16);

           const double tmpFDPart1_Rational_1_2 = 1.0/2.0;
           const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(tmpFDPart1_Rational_1_2);

           const double tmpFDPart1_Rational_1_4 = 1.0/4.0;
           const REAL_SIMD_ARRAY FDPart1_Rational_1_4 = ConstSIMD(tmpFDPart1_Rational_1_4);

           const double tmpFDPart1_Rational_3_8 = 3.0/8.0;
           const REAL_SIMD_ARRAY FDPart1_Rational_3_8 = ConstSIMD(tmpFDPart1_Rational_3_8);

           const REAL_SIMD_ARRAY FDPart1_0 = MulSIMD(FDPart1_Rational_1_2, invdx0);
           const REAL_SIMD_ARRAY FDPart1_1 = MulSIMD(FDPart1_Rational_1_2, invdx1);
           const REAL_SIMD_ARRAY FDPart1_2 = MulSIMD(FDPart1_Rational_1_2, invdx2);
           const REAL_SIMD_ARRAY FDPart1_3 = MulSIMD(FDPart1_Rational_3_8, bb);
           const REAL_SIMD_ARRAY FDPart1_4 = MulSIMD(FDPart1_Rational_3_8, bp);
           const REAL_SIMD_ARRAY FDPart1_5 = MulSIMD(FDPart1_Rational_3_8, fgauge);
           const REAL_SIMD_ARRAY FDPart1_6 = MulSIMD(FDPart1_Rational_3_8, fgaugep);
           const REAL_SIMD_ARRAY FDPart1_7 = MulSIMD(FDPart1_Rational_3_8, gg);
           const REAL_SIMD_ARRAY FDPart1_8 = MulSIMD(FDPart1_Rational_3_8, gp);
           const REAL_SIMD_ARRAY FDPart1_9 = MulSIMD(FDPart1_Rational_3_8, psibD0);
           const REAL_SIMD_ARRAY FDPart1_10 = MulSIMD(FDPart1_Rational_3_8, psibD1);
           const REAL_SIMD_ARRAY FDPart1_11 = MulSIMD(FDPart1_Rational_3_8, psibD2);
           const REAL_SIMD_ARRAY FDPart1_12 = MulSIMD(FDPart1_Rational_3_8, psifD0);
           const REAL_SIMD_ARRAY FDPart1_13 = MulSIMD(FDPart1_Rational_3_8, psifD1);
           const REAL_SIMD_ARRAY FDPart1_14 = MulSIMD(FDPart1_Rational_3_8, psifD2);
           const REAL_SIMD_ARRAY FDPart1_15 = MulSIMD(FDPart1_Rational_3_8, psigD0);
           const REAL_SIMD_ARRAY FDPart1_16 = MulSIMD(FDPart1_Rational_3_8, psigD1);
           const REAL_SIMD_ARRAY FDPart1_17 = MulSIMD(FDPart1_Rational_3_8, psigD2);
           const REAL_SIMD_ARRAY FDPart1_18 = MulSIMD(FDPart1_Rational_3_8, psiuD0);
           const REAL_SIMD_ARRAY FDPart1_19 = MulSIMD(FDPart1_Rational_3_8, psiuD1);
           const REAL_SIMD_ARRAY FDPart1_20 = MulSIMD(FDPart1_Rational_3_8, psiuD2);
           const REAL_SIMD_ARRAY FDPart1_21 = MulSIMD(FDPart1_Rational_3_8, up);
           const REAL_SIMD_ARRAY FDPart1_22 = MulSIMD(FDPart1_Rational_3_8, uu);
           const REAL_SIMD_ARRAY bb_dD0 = MulSIMD(FDPart1_0, SubSIMD(bb_i0p1_i1_i2, bb_i0m1_i1_i2));
           const REAL_SIMD_ARRAY bb_dD1 = MulSIMD(FDPart1_1, SubSIMD(bb_i0_i1p1_i2, bb_i0_i1m1_i2));
           const REAL_SIMD_ARRAY bb_dD2 = MulSIMD(FDPart1_2, SubSIMD(bb_i0_i1_i2p1, bb_i0_i1_i2m1));
           const REAL_SIMD_ARRAY bb_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(bb_i0m1_i1_i2, bb_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(bb_i0m2_i1_i2, bb_i0p2_i1_i2), FDPart1_3)));
           const REAL_SIMD_ARRAY bb_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(bb_i0_i1m1_i2, bb_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(bb_i0_i1m2_i2, bb_i0_i1p2_i2), FDPart1_3)));
           const REAL_SIMD_ARRAY bb_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(bb_i0_i1_i2m1, bb_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(bb_i0_i1_i2m2, bb_i0_i1_i2p2), FDPart1_3)));
           const REAL_SIMD_ARRAY bp_dD0 = MulSIMD(FDPart1_0, SubSIMD(bp_i0p1_i1_i2, bp_i0m1_i1_i2));
           const REAL_SIMD_ARRAY bp_dD1 = MulSIMD(FDPart1_1, SubSIMD(bp_i0_i1p1_i2, bp_i0_i1m1_i2));
           const REAL_SIMD_ARRAY bp_dD2 = MulSIMD(FDPart1_2, SubSIMD(bp_i0_i1_i2p1, bp_i0_i1_i2m1));
           const REAL_SIMD_ARRAY bp_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(bp_i0m1_i1_i2, bp_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(bp_i0m2_i1_i2, bp_i0p2_i1_i2), FDPart1_4)));
           const REAL_SIMD_ARRAY bp_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(bp_i0_i1m1_i2, bp_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(bp_i0_i1m2_i2, bp_i0_i1p2_i2), FDPart1_4)));
           const REAL_SIMD_ARRAY bp_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(bp_i0_i1_i2m1, bp_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(bp_i0_i1_i2m2, bp_i0_i1_i2p2), FDPart1_4)));
           const REAL_SIMD_ARRAY fgauge_dD0 = MulSIMD(FDPart1_0, SubSIMD(fgauge_i0p1_i1_i2, fgauge_i0m1_i1_i2));
           const REAL_SIMD_ARRAY fgauge_dD1 = MulSIMD(FDPart1_1, SubSIMD(fgauge_i0_i1p1_i2, fgauge_i0_i1m1_i2));
           const REAL_SIMD_ARRAY fgauge_dD2 = MulSIMD(FDPart1_2, SubSIMD(fgauge_i0_i1_i2p1, fgauge_i0_i1_i2m1));
           const REAL_SIMD_ARRAY fgauge_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(fgauge_i0m1_i1_i2, fgauge_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(fgauge_i0m2_i1_i2, fgauge_i0p2_i1_i2), FDPart1_5)));
           const REAL_SIMD_ARRAY fgauge_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(fgauge_i0_i1m1_i2, fgauge_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(fgauge_i0_i1m2_i2, fgauge_i0_i1p2_i2), FDPart1_5)));
           const REAL_SIMD_ARRAY fgauge_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(fgauge_i0_i1_i2m1, fgauge_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(fgauge_i0_i1_i2m2, fgauge_i0_i1_i2p2), FDPart1_5)));
           const REAL_SIMD_ARRAY fgaugep_dD0 = MulSIMD(FDPart1_0, SubSIMD(fgaugep_i0p1_i1_i2, fgaugep_i0m1_i1_i2));
           const REAL_SIMD_ARRAY fgaugep_dD1 = MulSIMD(FDPart1_1, SubSIMD(fgaugep_i0_i1p1_i2, fgaugep_i0_i1m1_i2));
           const REAL_SIMD_ARRAY fgaugep_dD2 = MulSIMD(FDPart1_2, SubSIMD(fgaugep_i0_i1_i2p1, fgaugep_i0_i1_i2m1));
           const REAL_SIMD_ARRAY fgaugep_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(fgaugep_i0m1_i1_i2, fgaugep_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(fgaugep_i0m2_i1_i2, fgaugep_i0p2_i1_i2), FDPart1_6)));
           const REAL_SIMD_ARRAY fgaugep_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(fgaugep_i0_i1m1_i2, fgaugep_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(fgaugep_i0_i1m2_i2, fgaugep_i0_i1p2_i2), FDPart1_6)));
           const REAL_SIMD_ARRAY fgaugep_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(fgaugep_i0_i1_i2m1, fgaugep_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(fgaugep_i0_i1_i2m2, fgaugep_i0_i1_i2p2), FDPart1_6)));
           const REAL_SIMD_ARRAY gg_dD0 = MulSIMD(FDPart1_0, SubSIMD(gg_i0p1_i1_i2, gg_i0m1_i1_i2));
           const REAL_SIMD_ARRAY gg_dD1 = MulSIMD(FDPart1_1, SubSIMD(gg_i0_i1p1_i2, gg_i0_i1m1_i2));
           const REAL_SIMD_ARRAY gg_dD2 = MulSIMD(FDPart1_2, SubSIMD(gg_i0_i1_i2p1, gg_i0_i1_i2m1));
           const REAL_SIMD_ARRAY gg_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(gg_i0m1_i1_i2, gg_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(gg_i0m2_i1_i2, gg_i0p2_i1_i2), FDPart1_7)));
           const REAL_SIMD_ARRAY gg_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(gg_i0_i1m1_i2, gg_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(gg_i0_i1m2_i2, gg_i0_i1p2_i2), FDPart1_7)));
           const REAL_SIMD_ARRAY gg_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(gg_i0_i1_i2m1, gg_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(gg_i0_i1_i2m2, gg_i0_i1_i2p2), FDPart1_7)));
           const REAL_SIMD_ARRAY gp_dD0 = MulSIMD(FDPart1_0, SubSIMD(gp_i0p1_i1_i2, gp_i0m1_i1_i2));
           const REAL_SIMD_ARRAY gp_dD1 = MulSIMD(FDPart1_1, SubSIMD(gp_i0_i1p1_i2, gp_i0_i1m1_i2));
           const REAL_SIMD_ARRAY gp_dD2 = MulSIMD(FDPart1_2, SubSIMD(gp_i0_i1_i2p1, gp_i0_i1_i2m1));
           const REAL_SIMD_ARRAY gp_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(gp_i0m1_i1_i2, gp_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(gp_i0m2_i1_i2, gp_i0p2_i1_i2), FDPart1_8)));
           const REAL_SIMD_ARRAY gp_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(gp_i0_i1m1_i2, gp_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(gp_i0_i1m2_i2, gp_i0_i1p2_i2), FDPart1_8)));
           const REAL_SIMD_ARRAY gp_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(gp_i0_i1_i2m1, gp_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(gp_i0_i1_i2m2, gp_i0_i1_i2p2), FDPart1_8)));
           const REAL_SIMD_ARRAY psibD_dD00 = MulSIMD(FDPart1_0, SubSIMD(psibD0_i0p1_i1_i2, psibD0_i0m1_i1_i2));
           const REAL_SIMD_ARRAY psibD_dD01 = MulSIMD(FDPart1_1, SubSIMD(psibD0_i0_i1p1_i2, psibD0_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psibD_dD02 = MulSIMD(FDPart1_2, SubSIMD(psibD0_i0_i1_i2p1, psibD0_i0_i1_i2m1));
           const REAL_SIMD_ARRAY psibD_dD11 = MulSIMD(FDPart1_1, SubSIMD(psibD1_i0_i1p1_i2, psibD1_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psibD_dD22 = MulSIMD(FDPart1_2, SubSIMD(psibD2_i0_i1_i2p1, psibD2_i0_i1_i2m1));
           const REAL_SIMD_ARRAY psibD_dKOD00 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psibD0_i0m1_i1_i2, psibD0_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psibD0_i0m2_i1_i2, psibD0_i0p2_i1_i2), FDPart1_9)));
           const REAL_SIMD_ARRAY psibD_dKOD01 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psibD0_i0_i1m1_i2, psibD0_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psibD0_i0_i1m2_i2, psibD0_i0_i1p2_i2), FDPart1_9)));
           const REAL_SIMD_ARRAY psibD_dKOD02 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psibD0_i0_i1_i2m1, psibD0_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psibD0_i0_i1_i2m2, psibD0_i0_i1_i2p2), FDPart1_9)));
           const REAL_SIMD_ARRAY psibD_dKOD10 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psibD1_i0m1_i1_i2, psibD1_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psibD1_i0m2_i1_i2, psibD1_i0p2_i1_i2), FDPart1_10)));
           const REAL_SIMD_ARRAY psibD_dKOD11 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psibD1_i0_i1m1_i2, psibD1_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psibD1_i0_i1m2_i2, psibD1_i0_i1p2_i2), FDPart1_10)));
           const REAL_SIMD_ARRAY psibD_dKOD12 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psibD1_i0_i1_i2m1, psibD1_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psibD1_i0_i1_i2m2, psibD1_i0_i1_i2p2), FDPart1_10)));
           const REAL_SIMD_ARRAY psibD_dKOD20 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psibD2_i0m1_i1_i2, psibD2_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psibD2_i0m2_i1_i2, psibD2_i0p2_i1_i2), FDPart1_11)));
           const REAL_SIMD_ARRAY psibD_dKOD21 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psibD2_i0_i1m1_i2, psibD2_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psibD2_i0_i1m2_i2, psibD2_i0_i1p2_i2), FDPart1_11)));
           const REAL_SIMD_ARRAY psibD_dKOD22 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psibD2_i0_i1_i2m1, psibD2_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psibD2_i0_i1_i2m2, psibD2_i0_i1_i2p2), FDPart1_11)));
           const REAL_SIMD_ARRAY psifD_dD00 = MulSIMD(FDPart1_0, SubSIMD(psifD0_i0p1_i1_i2, psifD0_i0m1_i1_i2));
           const REAL_SIMD_ARRAY psifD_dD01 = MulSIMD(FDPart1_1, SubSIMD(psifD0_i0_i1p1_i2, psifD0_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psifD_dD02 = MulSIMD(FDPart1_2, SubSIMD(psifD0_i0_i1_i2p1, psifD0_i0_i1_i2m1));
           const REAL_SIMD_ARRAY psifD_dD11 = MulSIMD(FDPart1_1, SubSIMD(psifD1_i0_i1p1_i2, psifD1_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psifD_dD22 = MulSIMD(FDPart1_2, SubSIMD(psifD2_i0_i1_i2p1, psifD2_i0_i1_i2m1));
           const REAL_SIMD_ARRAY psifD_dKOD00 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psifD0_i0m1_i1_i2, psifD0_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psifD0_i0m2_i1_i2, psifD0_i0p2_i1_i2), FDPart1_12)));
           const REAL_SIMD_ARRAY psifD_dKOD01 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psifD0_i0_i1m1_i2, psifD0_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psifD0_i0_i1m2_i2, psifD0_i0_i1p2_i2), FDPart1_12)));
           const REAL_SIMD_ARRAY psifD_dKOD02 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psifD0_i0_i1_i2m1, psifD0_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psifD0_i0_i1_i2m2, psifD0_i0_i1_i2p2), FDPart1_12)));
           const REAL_SIMD_ARRAY psifD_dKOD10 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psifD1_i0m1_i1_i2, psifD1_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psifD1_i0m2_i1_i2, psifD1_i0p2_i1_i2), FDPart1_13)));
           const REAL_SIMD_ARRAY psifD_dKOD11 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psifD1_i0_i1m1_i2, psifD1_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psifD1_i0_i1m2_i2, psifD1_i0_i1p2_i2), FDPart1_13)));
           const REAL_SIMD_ARRAY psifD_dKOD12 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psifD1_i0_i1_i2m1, psifD1_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psifD1_i0_i1_i2m2, psifD1_i0_i1_i2p2), FDPart1_13)));
           const REAL_SIMD_ARRAY psifD_dKOD20 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psifD2_i0m1_i1_i2, psifD2_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psifD2_i0m2_i1_i2, psifD2_i0p2_i1_i2), FDPart1_14)));
           const REAL_SIMD_ARRAY psifD_dKOD21 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psifD2_i0_i1m1_i2, psifD2_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psifD2_i0_i1m2_i2, psifD2_i0_i1p2_i2), FDPart1_14)));
           const REAL_SIMD_ARRAY psifD_dKOD22 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psifD2_i0_i1_i2m1, psifD2_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psifD2_i0_i1_i2m2, psifD2_i0_i1_i2p2), FDPart1_14)));
           const REAL_SIMD_ARRAY psigD_dD00 = MulSIMD(FDPart1_0, SubSIMD(psigD0_i0p1_i1_i2, psigD0_i0m1_i1_i2));
           const REAL_SIMD_ARRAY psigD_dD01 = MulSIMD(FDPart1_1, SubSIMD(psigD0_i0_i1p1_i2, psigD0_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psigD_dD02 = MulSIMD(FDPart1_2, SubSIMD(psigD0_i0_i1_i2p1, psigD0_i0_i1_i2m1));
           const REAL_SIMD_ARRAY psigD_dD11 = MulSIMD(FDPart1_1, SubSIMD(psigD1_i0_i1p1_i2, psigD1_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psigD_dD22 = MulSIMD(FDPart1_2, SubSIMD(psigD2_i0_i1_i2p1, psigD2_i0_i1_i2m1));
           const REAL_SIMD_ARRAY psigD_dKOD00 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psigD0_i0m1_i1_i2, psigD0_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psigD0_i0m2_i1_i2, psigD0_i0p2_i1_i2), FDPart1_15)));
           const REAL_SIMD_ARRAY psigD_dKOD01 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psigD0_i0_i1m1_i2, psigD0_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psigD0_i0_i1m2_i2, psigD0_i0_i1p2_i2), FDPart1_15)));
           const REAL_SIMD_ARRAY psigD_dKOD02 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psigD0_i0_i1_i2m1, psigD0_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psigD0_i0_i1_i2m2, psigD0_i0_i1_i2p2), FDPart1_15)));
           const REAL_SIMD_ARRAY psigD_dKOD10 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psigD1_i0m1_i1_i2, psigD1_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psigD1_i0m2_i1_i2, psigD1_i0p2_i1_i2), FDPart1_16)));
           const REAL_SIMD_ARRAY psigD_dKOD11 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psigD1_i0_i1m1_i2, psigD1_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psigD1_i0_i1m2_i2, psigD1_i0_i1p2_i2), FDPart1_16)));
           const REAL_SIMD_ARRAY psigD_dKOD12 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psigD1_i0_i1_i2m1, psigD1_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psigD1_i0_i1_i2m2, psigD1_i0_i1_i2p2), FDPart1_16)));
           const REAL_SIMD_ARRAY psigD_dKOD20 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psigD2_i0m1_i1_i2, psigD2_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psigD2_i0m2_i1_i2, psigD2_i0p2_i1_i2), FDPart1_17)));
           const REAL_SIMD_ARRAY psigD_dKOD21 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psigD2_i0_i1m1_i2, psigD2_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psigD2_i0_i1m2_i2, psigD2_i0_i1p2_i2), FDPart1_17)));
           const REAL_SIMD_ARRAY psigD_dKOD22 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psigD2_i0_i1_i2m1, psigD2_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psigD2_i0_i1_i2m2, psigD2_i0_i1_i2p2), FDPart1_17)));
           const REAL_SIMD_ARRAY psiuD_dD00 = MulSIMD(FDPart1_0, SubSIMD(psiuD0_i0p1_i1_i2, psiuD0_i0m1_i1_i2));
           const REAL_SIMD_ARRAY psiuD_dD01 = MulSIMD(FDPart1_1, SubSIMD(psiuD0_i0_i1p1_i2, psiuD0_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psiuD_dD02 = MulSIMD(FDPart1_2, SubSIMD(psiuD0_i0_i1_i2p1, psiuD0_i0_i1_i2m1));
           const REAL_SIMD_ARRAY psiuD_dD11 = MulSIMD(FDPart1_1, SubSIMD(psiuD1_i0_i1p1_i2, psiuD1_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psiuD_dD22 = MulSIMD(FDPart1_2, SubSIMD(psiuD2_i0_i1_i2p1, psiuD2_i0_i1_i2m1));
           const REAL_SIMD_ARRAY psiuD_dKOD00 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psiuD0_i0m1_i1_i2, psiuD0_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psiuD0_i0m2_i1_i2, psiuD0_i0p2_i1_i2), FDPart1_18)));
           const REAL_SIMD_ARRAY psiuD_dKOD01 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psiuD0_i0_i1m1_i2, psiuD0_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psiuD0_i0_i1m2_i2, psiuD0_i0_i1p2_i2), FDPart1_18)));
           const REAL_SIMD_ARRAY psiuD_dKOD02 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psiuD0_i0_i1_i2m1, psiuD0_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psiuD0_i0_i1_i2m2, psiuD0_i0_i1_i2p2), FDPart1_18)));
           const REAL_SIMD_ARRAY psiuD_dKOD10 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psiuD1_i0m1_i1_i2, psiuD1_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psiuD1_i0m2_i1_i2, psiuD1_i0p2_i1_i2), FDPart1_19)));
           const REAL_SIMD_ARRAY psiuD_dKOD11 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psiuD1_i0_i1m1_i2, psiuD1_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psiuD1_i0_i1m2_i2, psiuD1_i0_i1p2_i2), FDPart1_19)));
           const REAL_SIMD_ARRAY psiuD_dKOD12 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psiuD1_i0_i1_i2m1, psiuD1_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psiuD1_i0_i1_i2m2, psiuD1_i0_i1_i2p2), FDPart1_19)));
           const REAL_SIMD_ARRAY psiuD_dKOD20 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psiuD2_i0m1_i1_i2, psiuD2_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psiuD2_i0m2_i1_i2, psiuD2_i0p2_i1_i2), FDPart1_20)));
           const REAL_SIMD_ARRAY psiuD_dKOD21 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psiuD2_i0_i1m1_i2, psiuD2_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psiuD2_i0_i1m2_i2, psiuD2_i0_i1p2_i2), FDPart1_20)));
           const REAL_SIMD_ARRAY psiuD_dKOD22 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(psiuD2_i0_i1_i2m1, psiuD2_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(psiuD2_i0_i1_i2m2, psiuD2_i0_i1_i2p2), FDPart1_20)));
           const REAL_SIMD_ARRAY up_dD0 = MulSIMD(FDPart1_0, SubSIMD(up_i0p1_i1_i2, up_i0m1_i1_i2));
           const REAL_SIMD_ARRAY up_dD1 = MulSIMD(FDPart1_1, SubSIMD(up_i0_i1p1_i2, up_i0_i1m1_i2));
           const REAL_SIMD_ARRAY up_dD2 = MulSIMD(FDPart1_2, SubSIMD(up_i0_i1_i2p1, up_i0_i1_i2m1));
           const REAL_SIMD_ARRAY up_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(up_i0m1_i1_i2, up_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(up_i0m2_i1_i2, up_i0p2_i1_i2), FDPart1_21)));
           const REAL_SIMD_ARRAY up_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(up_i0_i1m1_i2, up_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(up_i0_i1m2_i2, up_i0_i1p2_i2), FDPart1_21)));
           const REAL_SIMD_ARRAY up_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(up_i0_i1_i2m1, up_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(up_i0_i1_i2m2, up_i0_i1_i2p2), FDPart1_21)));
           const REAL_SIMD_ARRAY uu_dD0 = MulSIMD(FDPart1_0, SubSIMD(uu_i0p1_i1_i2, uu_i0m1_i1_i2));
           const REAL_SIMD_ARRAY uu_dD1 = MulSIMD(FDPart1_1, SubSIMD(uu_i0_i1p1_i2, uu_i0_i1m1_i2));
           const REAL_SIMD_ARRAY uu_dD2 = MulSIMD(FDPart1_2, SubSIMD(uu_i0_i1_i2p1, uu_i0_i1_i2m1));
           const REAL_SIMD_ARRAY uu_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(uu_i0m1_i1_i2, uu_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(uu_i0m2_i1_i2, uu_i0p2_i1_i2), FDPart1_22)));
           const REAL_SIMD_ARRAY uu_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(uu_i0_i1m1_i2, uu_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(uu_i0_i1m2_i2, uu_i0_i1p2_i2), FDPart1_22)));
           const REAL_SIMD_ARRAY uu_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(uu_i0_i1_i2m1, uu_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(uu_i0_i1_i2m2, uu_i0_i1_i2p2), FDPart1_22)));
           /*
            * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
            */
           /*
            *  Original SymPy expressions:
            *  "[const REAL_SIMD_ARRAY __RHS_exp_0 = 5764607523034235*gg_dKOD0/288230376151711744 + psigD0/2 + 5764607523034235*gg_dKOD1/(288230376151711744*f0_of_xx0) + 5764607523034235*gg_dKOD2/(288230376151711744*f0_of_xx0*f1_of_xx1) - chip*gg/(2*chi) + gp/(2*chi),
            *    const REAL_SIMD_ARRAY __RHS_exp_1 = 5764607523034235*gp_dKOD0/288230376151711744 + (Rprime*chip*psigD0 - Rprime*chipp*gg + Rprime*chip**2*gg/chi - Rprime*chip*gp/chi - 2*chi*psigD0*rcoord/(1 - rcoord**2) - 2*chip*gg*rcoord/(1 - rcoord**2) - chip*gg_dD0 - chip*psigD0/2 + 2*gp*rcoord/(1 - rcoord**2) + gp_dD0 + chip**2*gg/(2*chi) - chip*gp/(2*chi) + Rprime*chi*costh*psigD1/(Radius**2*sinth) + Rprime*chi*psigD_dD11/Radius**2 + Rprime*chi*psigD_dD22/(Radius**2*sinth))/(2*Rprime - 1) + 5764607523034235*gp_dKOD1/(288230376151711744*f0_of_xx0) + 5764607523034235*gp_dKOD2/(288230376151711744*f0_of_xx0*f1_of_xx1),
            *    const REAL_SIMD_ARRAY __RHS_exp_2 = -bb*chip/(2*chi) + 5764607523034235*bb_dKOD0/288230376151711744 + 5764607523034235*bb_dKOD1/(288230376151711744*f0_of_xx0) + 5764607523034235*bb_dKOD2/(288230376151711744*f0_of_xx0*f1_of_xx1) + bp/(2*chi) + psibD0/2,
            *    const REAL_SIMD_ARRAY __RHS_exp_3 = 5764607523034235*bp_dKOD0/288230376151711744 + 5764607523034235*bp_dKOD1/(288230376151711744*f0_of_xx0) + 5764607523034235*bp_dKOD2/(288230376151711744*f0_of_xx0*f1_of_xx1) + (-Rprime*bb*chipp + Rprime*bb*chip**2/chi - Rprime*bp*chip/chi + Rprime*chip*psibD0 - Rprime*psifD0/2 - Rprime*psigD0**2/4 + Rprime*chip*fgauge/(2*chi) + Rprime*chip*gg*psigD0/(2*chi) - Rprime*fgaugep/(2*chi) - Rprime*gp*psigD0/(2*chi) - Rprime*chip**2*gg**2/(4*chi**2) + Rprime*chip*gg*gp/(2*chi**2) - Rprime*gp**2/(4*chi**2) - 2*bb*chip*rcoord/(1 - rcoord**2) + bb*chip**2/(2*chi) - bb_dD0*chip + 2*bp*rcoord/(1 - rcoord**2) - bp*chip/(2*chi) + bp_dD0 - 2*chi*psibD0*rcoord/(1 - rcoord**2) - chip*psibD0/2 + Rprime*chi*costh*psibD1/(Radius**2*sinth) + Rprime*chi*psibD_dD11/Radius**2 + Rprime*chi*psibD_dD22/(Radius**2*sinth))/(2*Rprime - 1),
            *    const REAL_SIMD_ARRAY __RHS_exp_4 = psiuD0/2 + 5764607523034235*uu_dKOD0/288230376151711744 + 5764607523034235*uu_dKOD1/(288230376151711744*f0_of_xx0) + 5764607523034235*uu_dKOD2/(288230376151711744*f0_of_xx0*f1_of_xx1) - chip*uu/(2*chi) + up/(2*chi),
            *    const REAL_SIMD_ARRAY __RHS_exp_5 = 5764607523034235*up_dKOD0/288230376151711744 + (Rprime*chip*psiuD0 - Rprime*chipp*uu - Rprime*psiuD0 + Rprime*chip**2*uu/chi - Rprime*chip*up/chi + Rprime*chip*uu/chi - Rprime*up/chi - 2*chi*psiuD0*rcoord/(1 - rcoord**2) - chip*psiuD0/2 - 2*chip*rcoord*uu/(1 - rcoord**2) - chip*uu_dD0 + 2*rcoord*up/(1 - rcoord**2) + up_dD0 + chip**2*uu/(2*chi) - chip*up/(2*chi) + Rprime*chi*costh*psiuD1/(Radius**2*sinth) + Rprime*chi*psiuD_dD11/Radius**2 + Rprime*chi*psiuD_dD22/(Radius**2*sinth))/(2*Rprime - 1) + 5764607523034235*up_dKOD1/(288230376151711744*f0_of_xx0) + 5764607523034235*up_dKOD2/(288230376151711744*f0_of_xx0*f1_of_xx1),
            *    const REAL_SIMD_ARRAY __RHS_exp_6 = 5764607523034235*fgauge_dKOD0/288230376151711744 + psifD0/2 + 5764607523034235*fgauge_dKOD1/(288230376151711744*f0_of_xx0) + 5764607523034235*fgauge_dKOD2/(288230376151711744*f0_of_xx0*f1_of_xx1) - chip*fgauge/(2*chi) + fgaugep/(2*chi),
            *    const REAL_SIMD_ARRAY __RHS_exp_7 = 5764607523034235*fgaugep_dKOD0/288230376151711744 + (Rprime*chip*psifD0 - Rprime*chipp*fgauge - Rprime*psifD0 - Rprime*psigD0**2/2 + Rprime*chip**2*fgauge/chi + Rprime*chip*fgauge/chi - Rprime*chip*fgaugep/chi + Rprime*chip*gg*psigD0/chi - Rprime*fgaugep/chi - Rprime*gp*psigD0/chi - Rprime*chip**2*gg**2/(2*chi**2) + Rprime*chip*gg*gp/chi**2 - Rprime*gp**2/(2*chi**2) - 2*chi*psifD0*rcoord/(1 - rcoord**2) - 2*chip*fgauge*rcoord/(1 - rcoord**2) - chip*fgauge_dD0 - chip*psifD0/2 + 2*fgaugep*rcoord/(1 - rcoord**2) + fgaugep_dD0 + chip**2*fgauge/(2*chi) - chip*fgaugep/(2*chi) + Rprime*chi*costh*psifD1/(Radius**2*sinth) + Rprime*chi*psifD_dD11/Radius**2 + Rprime*chi*psifD_dD22/(Radius**2*sinth))/(2*Rprime - 1) + 5764607523034235*fgaugep_dKOD1/(288230376151711744*f0_of_xx0) + 5764607523034235*fgaugep_dKOD2/(288230376151711744*f0_of_xx0*f1_of_xx1),
            *    const REAL_SIMD_ARRAY __RHS_exp_8 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_9 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_10 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_11 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_12 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_13 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_14 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_15 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_16 = Rprime*chip*psigD0/chi - 2*psigD0*rcoord/(1 - rcoord**2) - psigD_dD00 + 5764607523034235*psigD_dKOD00/288230376151711744 + 5764607523034235*psigD_dKOD01/(288230376151711744*f0_of_xx0) + 5764607523034235*psigD_dKOD02/(288230376151711744*f0_of_xx0*f1_of_xx1) - 2*chip*gg*rcoord/(chi*(1 - rcoord**2)) + 2*gp*rcoord/(chi*(1 - rcoord**2)) + Rprime*costh*psigD1/(Radius**2*sinth) + Rprime*psigD_dD11/Radius**2 + Rprime*psigD_dD22/(Radius**2*sinth),
            *    const REAL_SIMD_ARRAY __RHS_exp_17 = psigD_dD01/2 + 5764607523034235*psigD_dKOD10/288230376151711744 + 5764607523034235*psigD_dKOD11/(288230376151711744*f0_of_xx0) + 5764607523034235*psigD_dKOD12/(288230376151711744*f0_of_xx0*f1_of_xx1) - chip*gg_dD1/(2*chi) + gp_dD1/(2*chi),
            *    const REAL_SIMD_ARRAY __RHS_exp_18 = psigD_dD02/(2*sinth) + 5764607523034235*psigD_dKOD20/288230376151711744 + 5764607523034235*psigD_dKOD21/(288230376151711744*f0_of_xx0) + 5764607523034235*psigD_dKOD22/(288230376151711744*f0_of_xx0*f1_of_xx1) - chip*gg_dD2/(2*chi*sinth) + gp_dD2/(2*chi*sinth),
            *    const REAL_SIMD_ARRAY __RHS_exp_19 = Rprime*chip*psibD0/chi - Rprime*psifD0/(2*chi) - Rprime*psigD0**2/(4*chi) + Rprime*chip*fgauge/(2*chi**2) + Rprime*chip*gg*psigD0/(2*chi**2) - Rprime*fgaugep/(2*chi**2) - Rprime*gp*psigD0/(2*chi**2) - Rprime*chip**2*gg**2/(4*chi**3) + Rprime*chip*gg*gp/(2*chi**3) - Rprime*gp**2/(4*chi**3) - 2*bb*chip*rcoord/(chi*(1 - rcoord**2)) + 2*bp*rcoord/(chi*(1 - rcoord**2)) - 2*psibD0*rcoord/(1 - rcoord**2) - psibD_dD00 + 5764607523034235*psibD_dKOD00/288230376151711744 + 5764607523034235*psibD_dKOD01/(288230376151711744*f0_of_xx0) + 5764607523034235*psibD_dKOD02/(288230376151711744*f0_of_xx0*f1_of_xx1) + Rprime*costh*psibD1/(Radius**2*sinth) + Rprime*psibD_dD11/Radius**2 + Rprime*psibD_dD22/(Radius**2*sinth),
            *    const REAL_SIMD_ARRAY __RHS_exp_20 = -bb_dD1*chip/(2*chi) + bp_dD1/(2*chi) + psibD_dD01/2 + 5764607523034235*psibD_dKOD10/288230376151711744 + 5764607523034235*psibD_dKOD11/(288230376151711744*f0_of_xx0) + 5764607523034235*psibD_dKOD12/(288230376151711744*f0_of_xx0*f1_of_xx1),
            *    const REAL_SIMD_ARRAY __RHS_exp_21 = -bb_dD2*chip/(2*chi*sinth) + bp_dD2/(2*chi*sinth) + psibD_dD02/(2*sinth) + 5764607523034235*psibD_dKOD20/288230376151711744 + 5764607523034235*psibD_dKOD21/(288230376151711744*f0_of_xx0) + 5764607523034235*psibD_dKOD22/(288230376151711744*f0_of_xx0*f1_of_xx1),
            *    const REAL_SIMD_ARRAY __RHS_exp_22 = Rprime*chip*psiuD0/chi - Rprime*psiuD0/chi + Rprime*chip*uu/chi**2 - Rprime*up/chi**2 - 2*psiuD0*rcoord/(1 - rcoord**2) - psiuD_dD00 + 5764607523034235*psiuD_dKOD00/288230376151711744 + 5764607523034235*psiuD_dKOD01/(288230376151711744*f0_of_xx0) + 5764607523034235*psiuD_dKOD02/(288230376151711744*f0_of_xx0*f1_of_xx1) - 2*chip*rcoord*uu/(chi*(1 - rcoord**2)) + 2*rcoord*up/(chi*(1 - rcoord**2)) + Rprime*costh*psiuD1/(Radius**2*sinth) + Rprime*psiuD_dD11/Radius**2 + Rprime*psiuD_dD22/(Radius**2*sinth),
            *    const REAL_SIMD_ARRAY __RHS_exp_23 = psiuD_dD01/2 + 5764607523034235*psiuD_dKOD10/288230376151711744 + 5764607523034235*psiuD_dKOD11/(288230376151711744*f0_of_xx0) + 5764607523034235*psiuD_dKOD12/(288230376151711744*f0_of_xx0*f1_of_xx1) - chip*uu_dD1/(2*chi) + up_dD1/(2*chi),
            *    const REAL_SIMD_ARRAY __RHS_exp_24 = psiuD_dD02/(2*sinth) + 5764607523034235*psiuD_dKOD20/288230376151711744 + 5764607523034235*psiuD_dKOD21/(288230376151711744*f0_of_xx0) + 5764607523034235*psiuD_dKOD22/(288230376151711744*f0_of_xx0*f1_of_xx1) - chip*uu_dD2/(2*chi*sinth) + up_dD2/(2*chi*sinth),
            *    const REAL_SIMD_ARRAY __RHS_exp_25 = Rprime*chip*psifD0/chi - Rprime*psifD0/chi - Rprime*psigD0**2/(2*chi) + Rprime*chip*fgauge/chi**2 + Rprime*chip*gg*psigD0/chi**2 - Rprime*fgaugep/chi**2 - Rprime*gp*psigD0/chi**2 - Rprime*chip**2*gg**2/(2*chi**3) + Rprime*chip*gg*gp/chi**3 - Rprime*gp**2/(2*chi**3) - 2*psifD0*rcoord/(1 - rcoord**2) - psifD_dD00 + 5764607523034235*psifD_dKOD00/288230376151711744 + 5764607523034235*psifD_dKOD01/(288230376151711744*f0_of_xx0) + 5764607523034235*psifD_dKOD02/(288230376151711744*f0_of_xx0*f1_of_xx1) - 2*chip*fgauge*rcoord/(chi*(1 - rcoord**2)) + 2*fgaugep*rcoord/(chi*(1 - rcoord**2)) + Rprime*costh*psifD1/(Radius**2*sinth) + Rprime*psifD_dD11/Radius**2 + Rprime*psifD_dD22/(Radius**2*sinth),
            *    const REAL_SIMD_ARRAY __RHS_exp_26 = psifD_dD01/2 + 5764607523034235*psifD_dKOD10/288230376151711744 + 5764607523034235*psifD_dKOD11/(288230376151711744*f0_of_xx0) + 5764607523034235*psifD_dKOD12/(288230376151711744*f0_of_xx0*f1_of_xx1) - chip*fgauge_dD1/(2*chi) + fgaugep_dD1/(2*chi),
            *    const REAL_SIMD_ARRAY __RHS_exp_27 = psifD_dD02/(2*sinth) + 5764607523034235*psifD_dKOD20/288230376151711744 + 5764607523034235*psifD_dKOD21/(288230376151711744*f0_of_xx0) + 5764607523034235*psifD_dKOD22/(288230376151711744*f0_of_xx0*f1_of_xx1) - chip*fgauge_dD2/(2*chi*sinth) + fgaugep_dD2/(2*chi*sinth)]"
            */
           const double tmpFDPart3_Integer_0 = 0.0;
           const REAL_SIMD_ARRAY FDPart3_Integer_0 = ConstSIMD(tmpFDPart3_Integer_0);

           const double tmpFDPart3_Integer_1 = 1.0;
           const REAL_SIMD_ARRAY FDPart3_Integer_1 = ConstSIMD(tmpFDPart3_Integer_1);

           const double tmpFDPart3_Integer_2 = 2.0;
           const REAL_SIMD_ARRAY FDPart3_Integer_2 = ConstSIMD(tmpFDPart3_Integer_2);

           const double tmpFDPart3_NegativeOne_ = -1.0;
           const REAL_SIMD_ARRAY FDPart3_NegativeOne_ = ConstSIMD(tmpFDPart3_NegativeOne_);

           const double tmpFDPart3_Rational_1_2 = 1.0/2.0;
           const REAL_SIMD_ARRAY FDPart3_Rational_1_2 = ConstSIMD(tmpFDPart3_Rational_1_2);

           const double tmpFDPart3_Rational_1_4 = 1.0/4.0;
           const REAL_SIMD_ARRAY FDPart3_Rational_1_4 = ConstSIMD(tmpFDPart3_Rational_1_4);

           const double tmpFDPart3_Rational_5764607523034235_288230376151711744 = 2.0*5764607523034235.0/2.8823037615171174e+16;
           const REAL_SIMD_ARRAY FDPart3_Rational_5764607523034235_288230376151711744 = ConstSIMD(tmpFDPart3_Rational_5764607523034235_288230376151711744);





//----------------------------------------------------------------------------------------------------
           // ADDED FUNCTIONS FOR EVANS METHOD
           const REAL_SIMD_ARRAY rcoord_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(RCOORDGF, i0+1,i1,i2)]);
           const REAL_SIMD_ARRAY rcoord_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(RCOORDGF, i0-1,i1,i2)]);
           const REAL_SIMD_ARRAY radiussquared = MulSIMD(rcoord,rcoord);
           const REAL_SIMD_ARRAY radiussquaredp1 = MulSIMD(rcoord_i0p1_i1_i2,rcoord_i0p1_i1_i2);
           const REAL_SIMD_ARRAY radiuscubedp1 = MulSIMD(rcoord_i0p1_i1_i2,radiussquaredp1);
           const REAL_SIMD_ARRAY radiussquaredm1 = MulSIMD(rcoord_i0m1_i1_i2,rcoord_i0m1_i1_i2);
           const REAL_SIMD_ARRAY radiuscubedm1 = MulSIMD(rcoord_i0m1_i1_i2,radiussquaredm1);
           const REAL_SIMD_ARRAY rcubeddif = SubSIMD(radiuscubedp1,radiuscubedm1);
           const double tmpFD_Integer_3 = 3.0;
           const REAL_SIMD_ARRAY FD_Integer_3 = ConstSIMD(tmpFD_Integer_3);

           // Evans finite difference of gg
           const REAL_SIMD_ARRAY rsggp1 = MulSIMD(radiussquaredp1,gg_i0p1_i1_i2);
           const REAL_SIMD_ARRAY rsggm1 = MulSIMD(radiussquaredm1,gg_i0m1_i1_i2);
           const REAL_SIMD_ARRAY rsggdif = SubSIMD(rsggp1,rsggm1);
           const REAL_SIMD_ARRAY gg_dED0 = MulSIMD(FD_Integer_3,DivSIMD(rsggdif,rcubeddif));
           //const REAL_SIMD_ARRAY gg_dED0 = MulSIMD(invdx0,DivSIMD(rsggdif,radiussquared));
           const REAL_SIMD_ARRAY ggor = MulSIMD(FDPart3_Rational_1_2,SubSIMD(gg_dED0,gg_dD0));
           const REAL_SIMD_ARRAY ggor_gp = MulSIMD(ggor, MulSIMD(FDPart3_NegativeOne_,chip) );
           const REAL_SIMD_ARRAY ggor_gm = MulSIMD(ggor, MulSIMD(FDPart3_NegativeOne_,DivSIMD(chip,chi)) );

           // Evans finite difference of gp
           const REAL_SIMD_ARRAY rsgpp1 = MulSIMD(radiussquaredp1,gp_i0p1_i1_i2);
           const REAL_SIMD_ARRAY rsgpm1 = MulSIMD(radiussquaredm1,gp_i0m1_i1_i2);
           const REAL_SIMD_ARRAY rsgpdif = SubSIMD(rsgpp1,rsgpm1);
           const REAL_SIMD_ARRAY gp_dED0 = MulSIMD(FD_Integer_3,DivSIMD(rsgpdif,rcubeddif));
           //const REAL_SIMD_ARRAY gp_dED0 = MulSIMD(invdx0,DivSIMD(rsgpdif,radiussquared));
           const REAL_SIMD_ARRAY gpor = MulSIMD(FDPart3_Rational_1_2,SubSIMD(gp_dED0,gp_dD0));
           const REAL_SIMD_ARRAY gpor_gm = DivSIMD(gpor,chi);

           // Evans finite difference of psig
           const REAL_SIMD_ARRAY rspsigD0p1 = MulSIMD(radiussquaredp1,psigD0_i0p1_i1_i2);
           const REAL_SIMD_ARRAY rspsigD0m1 = MulSIMD(radiussquaredm1,psigD0_i0m1_i1_i2);
           const REAL_SIMD_ARRAY rspsigD0dif = SubSIMD(rspsigD0p1,rspsigD0m1);
           const REAL_SIMD_ARRAY psigD_dED00 = MulSIMD(FD_Integer_3,DivSIMD(rspsigD0dif,rcubeddif));
           const REAL_SIMD_ARRAY psigor = MulSIMD(FDPart3_Rational_1_2,SubSIMD(psigD_dED00,psigD_dD00));
           const REAL_SIMD_ARRAY psigor_gp = MulSIMD(psigor, MulSIMD(FDPart3_NegativeOne_,chi) );
           const REAL_SIMD_ARRAY psigor_gm = MulSIMD(psigor,FDPart3_NegativeOne_);


           // Evans finite difference of bb
           const REAL_SIMD_ARRAY rsbbp1 = MulSIMD(radiussquaredp1,bb_i0p1_i1_i2);
           const REAL_SIMD_ARRAY rsbbm1 = MulSIMD(radiussquaredm1,bb_i0m1_i1_i2);
           const REAL_SIMD_ARRAY rsbbdif = SubSIMD(rsbbp1,rsbbm1);
           const REAL_SIMD_ARRAY bb_dED0 = MulSIMD(FD_Integer_3,DivSIMD(rsbbdif,rcubeddif));
           //const REAL_SIMD_ARRAY bb_dED0 = MulSIMD(invdx0,DivSIMD(rsbbdif,radiussquared));
           const REAL_SIMD_ARRAY bbor = MulSIMD(FDPart3_Rational_1_2,SubSIMD(bb_dED0,bb_dD0));
           const REAL_SIMD_ARRAY bbor_bp = MulSIMD(bbor, MulSIMD(FDPart3_NegativeOne_,chip) );
           const REAL_SIMD_ARRAY bbor_bm = MulSIMD(bbor, MulSIMD(FDPart3_NegativeOne_,DivSIMD(chip,chi)) );

           // Evans finite difference of bp
           const REAL_SIMD_ARRAY rsbpp1 = MulSIMD(radiussquaredp1,bp_i0p1_i1_i2);
           const REAL_SIMD_ARRAY rsbpm1 = MulSIMD(radiussquaredm1,bp_i0m1_i1_i2);
           const REAL_SIMD_ARRAY rsbpdif = SubSIMD(rsbpp1,rsbpm1);
           const REAL_SIMD_ARRAY bp_dED0 = MulSIMD(FD_Integer_3,DivSIMD(rsbpdif,rcubeddif));
           //const REAL_SIMD_ARRAY bp_dED0 = MulSIMD(invdx0,DivSIMD(rsbpdif,radiussquared));
           const REAL_SIMD_ARRAY bpor = MulSIMD(FDPart3_Rational_1_2,SubSIMD(bp_dED0,bp_dD0));
           const REAL_SIMD_ARRAY bpor_bm = DivSIMD(bpor,chi);

           // Evans finite difference of psib
           const REAL_SIMD_ARRAY rspsibD0p1 = MulSIMD(radiussquaredp1,psibD0_i0p1_i1_i2);
           const REAL_SIMD_ARRAY rspsibD0m1 = MulSIMD(radiussquaredm1,psibD0_i0m1_i1_i2);
           const REAL_SIMD_ARRAY rspsibD0dif = SubSIMD(rspsibD0p1,rspsibD0m1);
           const REAL_SIMD_ARRAY psibD_dED00 = MulSIMD(FD_Integer_3,DivSIMD(rspsibD0dif,rcubeddif));
           const REAL_SIMD_ARRAY psibor = MulSIMD(FDPart3_Rational_1_2,SubSIMD(psibD_dED00,psibD_dD00));
           const REAL_SIMD_ARRAY psibor_bp = MulSIMD(psibor, MulSIMD(FDPart3_NegativeOne_,chi) );
           const REAL_SIMD_ARRAY psibor_bm = MulSIMD(psibor,FDPart3_NegativeOne_);


           // Evans finite difference of uu
           const REAL_SIMD_ARRAY rsuup1 = MulSIMD(radiussquaredp1,uu_i0p1_i1_i2);
           const REAL_SIMD_ARRAY rsuum1 = MulSIMD(radiussquaredm1,uu_i0m1_i1_i2);
           const REAL_SIMD_ARRAY rsuudif = SubSIMD(rsuup1,rsuum1);
           const REAL_SIMD_ARRAY uu_dED0 = MulSIMD(FD_Integer_3,DivSIMD(rsuudif,rcubeddif));
           //const REAL_SIMD_ARRAY uu_dED0 = MulSIMD(invdx0,DivSIMD(rsuudif,radiussquared));
           const REAL_SIMD_ARRAY uuor = MulSIMD(FDPart3_Rational_1_2,SubSIMD(uu_dED0,uu_dD0));
           const REAL_SIMD_ARRAY uuor_up = MulSIMD(uuor, MulSIMD(FDPart3_NegativeOne_,chip) );
           const REAL_SIMD_ARRAY uuor_um = MulSIMD(uuor, MulSIMD(FDPart3_NegativeOne_,DivSIMD(chip,chi)) );

           // Evans finite difference of up
           const REAL_SIMD_ARRAY rsupp1 = MulSIMD(radiussquaredp1,up_i0p1_i1_i2);
           const REAL_SIMD_ARRAY rsupm1 = MulSIMD(radiussquaredm1,up_i0m1_i1_i2);
           const REAL_SIMD_ARRAY rsupdif = SubSIMD(rsupp1,rsupm1);
           const REAL_SIMD_ARRAY up_dED0 = MulSIMD(FD_Integer_3,DivSIMD(rsupdif,rcubeddif));
           //const REAL_SIMD_ARRAY up_dED0 = MulSIMD(invdx0,DivSIMD(rsupdif,radiussquared));
           const REAL_SIMD_ARRAY upor = MulSIMD(FDPart3_Rational_1_2,SubSIMD(up_dED0,up_dD0));
           const REAL_SIMD_ARRAY upor_um = DivSIMD(upor,chi);

           // Evans finite difference of psiu
           const REAL_SIMD_ARRAY rspsiuD0p1 = MulSIMD(radiussquaredp1,psiuD0_i0p1_i1_i2);
           const REAL_SIMD_ARRAY rspsiuD0m1 = MulSIMD(radiussquaredm1,psiuD0_i0m1_i1_i2);
           const REAL_SIMD_ARRAY rspsiuD0dif = SubSIMD(rspsiuD0p1,rspsiuD0m1);
           const REAL_SIMD_ARRAY psiuD_dED00 = MulSIMD(FD_Integer_3,DivSIMD(rspsiuD0dif,rcubeddif));
           const REAL_SIMD_ARRAY psiuor = MulSIMD(FDPart3_Rational_1_2,SubSIMD(psiuD_dED00,psiuD_dD00));
           const REAL_SIMD_ARRAY psiuor_up = MulSIMD(psiuor, MulSIMD(FDPart3_NegativeOne_,chi) );
           const REAL_SIMD_ARRAY psiuor_um = MulSIMD(psiuor,FDPart3_NegativeOne_);


           // Evans finite difference of fgauge
           const REAL_SIMD_ARRAY rsfgaugep1 = MulSIMD(radiussquaredp1,fgauge_i0p1_i1_i2);
           const REAL_SIMD_ARRAY rsfgaugem1 = MulSIMD(radiussquaredm1,fgauge_i0m1_i1_i2);
           const REAL_SIMD_ARRAY rsfgaugedif = SubSIMD(rsfgaugep1,rsfgaugem1);
           const REAL_SIMD_ARRAY fgauge_dED0 = MulSIMD(FD_Integer_3,DivSIMD(rsfgaugedif,rcubeddif));
           //const REAL_SIMD_ARRAY fgauge_dED0 = MulSIMD(invdx0,DivSIMD(rsfgaugedif,radiussquared));
           const REAL_SIMD_ARRAY fgaugeor = MulSIMD(FDPart3_Rational_1_2,SubSIMD(fgauge_dED0,fgauge_dD0));
           const REAL_SIMD_ARRAY fgaugeor_fp = MulSIMD(fgaugeor, MulSIMD(FDPart3_NegativeOne_,chip) );
           const REAL_SIMD_ARRAY fgaugeor_fm = MulSIMD(fgaugeor, MulSIMD(FDPart3_NegativeOne_,DivSIMD(chip,chi)) );

           // Evans finite difference of fgaugep
           const REAL_SIMD_ARRAY rsfgaugepp1 = MulSIMD(radiussquaredp1,fgaugep_i0p1_i1_i2);
           const REAL_SIMD_ARRAY rsfgaugepm1 = MulSIMD(radiussquaredm1,fgaugep_i0m1_i1_i2);
           const REAL_SIMD_ARRAY rsfgaugepdif = SubSIMD(rsfgaugepp1,rsfgaugepm1);
           const REAL_SIMD_ARRAY fgaugep_dED0 = MulSIMD(FD_Integer_3,DivSIMD(rsfgaugepdif,rcubeddif));
           //const REAL_SIMD_ARRAY fgaugep_dED0 = MulSIMD(invdx0,DivSIMD(rsfgaugepdif,radiussquared));
           const REAL_SIMD_ARRAY fgaugepor = MulSIMD(FDPart3_Rational_1_2,SubSIMD(fgaugep_dED0,fgaugep_dD0));
           const REAL_SIMD_ARRAY fgaugepor_fm = DivSIMD(fgaugepor,chi);

           // Evans finite difference of psif
           const REAL_SIMD_ARRAY rspsifD0p1 = MulSIMD(radiussquaredp1,psifD0_i0p1_i1_i2);
           const REAL_SIMD_ARRAY rspsifD0m1 = MulSIMD(radiussquaredm1,psifD0_i0m1_i1_i2);
           const REAL_SIMD_ARRAY rspsifD0dif = SubSIMD(rspsifD0p1,rspsifD0m1);
           const REAL_SIMD_ARRAY psifD_dED00 = MulSIMD(FD_Integer_3,DivSIMD(rspsifD0dif,rcubeddif));
           const REAL_SIMD_ARRAY psifor = MulSIMD(FDPart3_Rational_1_2,SubSIMD(psifD_dED00,psifD_dD00));
           const REAL_SIMD_ARRAY psifor_fp = MulSIMD(psifor, MulSIMD(FDPart3_NegativeOne_,chi) );
           const REAL_SIMD_ARRAY psifor_fm = MulSIMD(psifor,FDPart3_NegativeOne_);

//--------------------------------------------------------------------------------------------




           const REAL_SIMD_ARRAY FDPart3_1 = DivSIMD(FDPart3_Integer_1, chi);
           const REAL_SIMD_ARRAY FDPart3_2 = MulSIMD(FDPart3_1, FDPart3_Rational_1_2);
           const REAL_SIMD_ARRAY FDPart3_4 = MulSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(DivSIMD(FDPart3_Integer_1, f1_of_xx1), f0_of_xx0));
           const REAL_SIMD_ARRAY FDPart3_6 = DivSIMD(FDPart3_Integer_1, FusedMulAddSIMD(FDPart3_Integer_2, Rprime, FDPart3_NegativeOne_));
           const REAL_SIMD_ARRAY FDPart3_8 = MulSIMD(Rprime, chip);
           const REAL_SIMD_ARRAY FDPart3_9 = MulSIMD(FDPart3_8, psigD0);
           const REAL_SIMD_ARRAY FDPart3_11 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(Rprime, chipp));
           const REAL_SIMD_ARRAY FDPart3_12 = DivSIMD(Rprime, MulSIMD(Radius, Radius));
           const REAL_SIMD_ARRAY FDPart3_14 = MulSIMD(FDPart3_1, Rprime);
           const REAL_SIMD_ARRAY FDPart3_15 = MulSIMD(FDPart3_14, MulSIMD(FDPart3_NegativeOne_, chip));
           const REAL_SIMD_ARRAY FDPart3_16 = MulSIMD(chip, chip);
           const REAL_SIMD_ARRAY FDPart3_18 = DivSIMD(FDPart3_Integer_1, sinth);
           const REAL_SIMD_ARRAY FDPart3_19 = MulSIMD(FDPart3_12, FDPart3_18);
           const REAL_SIMD_ARRAY FDPart3_21 = MulSIMD(FDPart3_19, costh);
           const REAL_SIMD_ARRAY FDPart3_23 = DivSIMD(FDPart3_Integer_1, NegFusedMulAddSIMD(rcoord, rcoord, FDPart3_Integer_1));
           const REAL_SIMD_ARRAY FDPart3_24 = MulSIMD(FDPart3_23, MulSIMD(FDPart3_Integer_2, rcoord));
           const REAL_SIMD_ARRAY FDPart3_29 = MulSIMD(psigD0, psigD0);
           const REAL_SIMD_ARRAY FDPart3_32 = DivSIMD(FDPart3_Integer_1, MulSIMD(chi, chi));
           const REAL_SIMD_ARRAY FDPart3_33 = MulSIMD(FDPart3_32, FDPart3_8);
           const REAL_SIMD_ARRAY FDPart3_35 = MulSIMD(FDPart3_33, MulSIMD(gg, gp));
           const REAL_SIMD_ARRAY FDPart3_36 = MulSIMD(FDPart3_1, FDPart3_8);
           const REAL_SIMD_ARRAY FDPart3_37 = MulSIMD(FDPart3_Rational_1_2, MulSIMD(psigD0, gg));
           const REAL_SIMD_ARRAY FDPart3_39 = MulSIMD(gp, gp);
           const REAL_SIMD_ARRAY FDPart3_42 = MulSIMD(gg, gg);
           const REAL_SIMD_ARRAY FDPart3_43 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(psiuD0, Rprime));
           const REAL_SIMD_ARRAY FDPart3_45 = MulSIMD(FDPart3_1, MulSIMD(FDPart3_NegativeOne_, Rprime));
           const REAL_SIMD_ARRAY FDPart3_51 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(psifD0, Rprime));
           const REAL_SIMD_ARRAY FDPart3_60 = MulSIMD(FDPart3_18, FDPart3_Rational_1_2);
           const REAL_SIMD_ARRAY FDPart3_61 = MulSIMD(FDPart3_18, FDPart3_2);
           const REAL_SIMD_ARRAY FDPart3_63 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(chi, chi), chi));
           const REAL_SIMD_ARRAY FDPart3_64 = MulSIMD(MulSIMD(FDPart3_63, FDPart3_8), MulSIMD(gg, gp));
           const REAL_SIMD_ARRAY FDPart3_65 = MulSIMD(FDPart3_32, MulSIMD(FDPart3_NegativeOne_, Rprime));
           const REAL_SIMD_ARRAY __RHS_exp_0 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, gg_dKOD0, FusedMulAddSIMD(FDPart3_4, gg_dKOD2, FusedMulAddSIMD(FDPart3_Rational_1_2, psigD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(gg_dKOD1, f0_of_xx0), FusedMulSubSIMD(FDPart3_2, gp, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(chip, gg)))))));
           const REAL_SIMD_ARRAY __RHS_exp_1 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, gp_dKOD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(gp_dKOD1, f0_of_xx0), FusedMulAddSIMD(FDPart3_4, gp_dKOD2, MulSIMD( /*THIS MulSIMD IS THE 1/(2R'-1) FACTOR*/ FDPart3_6, AddSIMD(ggor_gp, AddSIMD(gpor, AddSIMD(psigor_gp,   FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(MulSIMD(FDPart3_1, FDPart3_NegativeOne_), MulSIMD(chip, gp)), NegFusedMulAddSIMD(psigD0, MulSIMD(MulSIMD(FDPart3_23, FDPart3_Integer_2), MulSIMD(chi, rcoord)), FusedMulAddSIMD(FDPart3_21, MulSIMD(psigD1, chi), NegFusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(psigD0, chip), FusedMulAddSIMD(FDPart3_16, MulSIMD(FDPart3_2, gg), FusedMulAddSIMD(FDPart3_19, MulSIMD(psigD_dD22, chi), FusedMulAddSIMD(FDPart3_12, MulSIMD(psigD_dD11, chi), FusedMulAddSIMD(FDPart3_14, MulSIMD(FDPart3_16, gg), FusedMulAddSIMD(FDPart3_15, gp, FusedMulAddSIMD(FDPart3_24, gp, AddSIMD(FusedMulAddSIMD(FDPart3_11, gg, gp_dD0), NegFusedMulAddSIMD(chip, MulSIMD(MulSIMD(FDPart3_23, FDPart3_Integer_2), MulSIMD(gg, rcoord)), NegFusedMulAddSIMD(gg_dD0, chip, FDPart3_9))))))))))))))))))));
           const REAL_SIMD_ARRAY __RHS_exp_2 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, bb_dKOD0, FusedMulAddSIMD(FDPart3_4, bb_dKOD2, FusedMulAddSIMD(FDPart3_Rational_1_2, psibD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(bb_dKOD1, f0_of_xx0), FusedMulSubSIMD(FDPart3_2, bp, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(bb, chip)))))));
           const REAL_SIMD_ARRAY __RHS_exp_3 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, bp_dKOD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(bp_dKOD1, f0_of_xx0), FusedMulAddSIMD(FDPart3_4, bp_dKOD2, MulSIMD( /*THIS MulSIMD IS THE 1/(2R'-1) FACTOR*/ FDPart3_6, AddSIMD(bbor_bp, AddSIMD(bpor, AddSIMD(psibor_bp,    NegFusedMulAddSIMD(psibD0, MulSIMD(MulSIMD(FDPart3_23, FDPart3_Integer_2), MulSIMD(chi, rcoord)), FusedMulAddSIMD(MulSIMD(FDPart3_Rational_1_2, psigD0), MulSIMD(MulSIMD(FDPart3_1, FDPart3_NegativeOne_), MulSIMD(Rprime, gp)), NegFusedMulAddSIMD(FDPart3_42, MulSIMD(MulSIMD(FDPart3_16, FDPart3_32), MulSIMD(FDPart3_Rational_1_4, Rprime)), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(MulSIMD(FDPart3_1, FDPart3_NegativeOne_), MulSIMD(bp, chip)), SubSIMD(NegFusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(psifD0, Rprime), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(MulSIMD(FDPart3_1, FDPart3_NegativeOne_), MulSIMD(Rprime, fgaugep)), NegFusedMulAddSIMD(FDPart3_29, MulSIMD(FDPart3_Rational_1_4, Rprime), NegFusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(psibD0, chip), FusedMulAddSIMD(FDPart3_2, MulSIMD(FDPart3_8, fgauge), FusedMulAddSIMD(FDPart3_21, MulSIMD(psibD1, chi), FusedMulAddSIMD(FDPart3_16, MulSIMD(FDPart3_2, bb), FusedMulAddSIMD(FDPart3_19, MulSIMD(psibD_dD22, chi), FusedMulAddSIMD(FDPart3_12, MulSIMD(psibD_dD11, chi), FusedMulAddSIMD(FDPart3_14, MulSIMD(FDPart3_16, bb), FusedMulAddSIMD(FDPart3_36, FDPart3_37, FusedMulAddSIMD(FDPart3_8, psibD0, FusedMulAddSIMD(FDPart3_24, bp, FusedMulAddSIMD(FDPart3_35, FDPart3_Rational_1_2, FusedMulAddSIMD(FDPart3_11, bb, FusedMulAddSIMD(FDPart3_15, bp, NegFusedMulAddSIMD(bb, MulSIMD(MulSIMD(FDPart3_23, FDPart3_Integer_2), MulSIMD(chip, rcoord)), NegFusedMulAddSIMD(bb_dD0, chip, bp_dD0)))))))))))))))))), MulSIMD(MulSIMD(FDPart3_32, FDPart3_39), MulSIMD(FDPart3_Rational_1_4, Rprime))))))))))))));
           const REAL_SIMD_ARRAY __RHS_exp_4 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, uu_dKOD0, FusedMulAddSIMD(FDPart3_4, uu_dKOD2, FusedMulAddSIMD(FDPart3_Rational_1_2, psiuD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(uu_dKOD1, f0_of_xx0), FusedMulSubSIMD(FDPart3_2, up, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(chip, uu)))))));
           const REAL_SIMD_ARRAY __RHS_exp_5 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, up_dKOD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(up_dKOD1, f0_of_xx0), FusedMulAddSIMD(FDPart3_4, up_dKOD2, MulSIMD( /*THIS MulSIMD IS THE 1/(2R'-1) FACTOR*/ FDPart3_6, AddSIMD(uuor_up, AddSIMD(upor, AddSIMD(psiuor_up,     NegFusedMulAddSIMD(psiuD0, MulSIMD(MulSIMD(FDPart3_23, FDPart3_Integer_2), MulSIMD(chi, rcoord)), NegFusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(psiuD0, chip), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(MulSIMD(FDPart3_1, FDPart3_NegativeOne_), MulSIMD(chip, up)), FusedMulAddSIMD(FDPart3_19, MulSIMD(psiuD_dD22, chi), FusedMulAddSIMD(FDPart3_21, MulSIMD(psiuD1, chi), FusedMulAddSIMD(FDPart3_14, MulSIMD(FDPart3_16, uu), FusedMulAddSIMD(FDPart3_16, MulSIMD(FDPart3_2, uu), FusedMulAddSIMD(FDPart3_8, psiuD0, FusedMulAddSIMD(FDPart3_12, MulSIMD(psiuD_dD11, chi), FusedMulAddSIMD(FDPart3_36, uu, FusedMulAddSIMD(FDPart3_45, up, FusedMulAddSIMD(FDPart3_15, up, FusedMulAddSIMD(FDPart3_24, up, AddSIMD(FusedMulAddSIMD(FDPart3_11, uu, up_dD0), NegFusedMulAddSIMD(chip, MulSIMD(MulSIMD(FDPart3_23, FDPart3_Integer_2), MulSIMD(rcoord, uu)), NegFusedMulAddSIMD(uu_dD0, chip, FDPart3_43)))))))))))))))))))))));
           const REAL_SIMD_ARRAY __RHS_exp_6 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, fgauge_dKOD0, FusedMulAddSIMD(FDPart3_4, fgauge_dKOD2, FusedMulAddSIMD(FDPart3_Rational_1_2, psifD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(fgauge_dKOD1, f0_of_xx0), FusedMulSubSIMD(FDPart3_2, fgaugep, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(chip, fgauge)))))));
           const REAL_SIMD_ARRAY __RHS_exp_7 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, fgaugep_dKOD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(fgaugep_dKOD1, f0_of_xx0), FusedMulAddSIMD(FDPart3_4, fgaugep_dKOD2, MulSIMD( /*THIS MulSIMD IS THE 1/(2R'-1) FACTOR*/ FDPart3_6, AddSIMD(fgaugeor_fp, AddSIMD(fgaugepor, AddSIMD(psifor_fp,     NegFusedMulAddSIMD(psifD0, MulSIMD(MulSIMD(FDPart3_23, FDPart3_Integer_2), MulSIMD(chi, rcoord)), SubSIMD(NegFusedMulAddSIMD(FDPart3_42, MulSIMD(MulSIMD(FDPart3_16, FDPart3_32), MulSIMD(FDPart3_Rational_1_2, Rprime)), NegFusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(psifD0, chip), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(MulSIMD(FDPart3_1, FDPart3_NegativeOne_), MulSIMD(chip, fgaugep)), FusedMulAddSIMD(FDPart3_45, MulSIMD(psigD0, gp), NegFusedMulAddSIMD(FDPart3_29, MulSIMD(FDPart3_Rational_1_2, Rprime), FusedMulAddSIMD(FDPart3_19, MulSIMD(psifD_dD22, chi), FusedMulAddSIMD(FDPart3_21, MulSIMD(psifD1, chi), FusedMulAddSIMD(FDPart3_14, MulSIMD(FDPart3_16, fgauge), FusedMulAddSIMD(FDPart3_16, MulSIMD(FDPart3_2, fgauge), FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_9, gg), FusedMulAddSIMD(FDPart3_12, MulSIMD(psifD_dD11, chi), FusedMulAddSIMD(FDPart3_45, fgaugep, FusedMulAddSIMD(FDPart3_8, psifD0, FusedMulAddSIMD(FDPart3_24, fgaugep, FusedMulAddSIMD(FDPart3_36, fgauge, FusedMulAddSIMD(FDPart3_11, fgauge, FusedMulAddSIMD(FDPart3_15, fgaugep, AddSIMD(AddSIMD(FDPart3_51, fgaugep_dD0), NegFusedMulAddSIMD(chip, MulSIMD(MulSIMD(FDPart3_23, FDPart3_Integer_2), MulSIMD(fgauge, rcoord)), NegFusedMulAddSIMD(fgauge_dD0, chip, FDPart3_35)))))))))))))))))))), MulSIMD(MulSIMD(FDPart3_32, FDPart3_39), MulSIMD(FDPart3_Rational_1_2, Rprime)))))))))));
           const REAL_SIMD_ARRAY __RHS_exp_8 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_9 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_10 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_11 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_12 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_13 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_14 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_15 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_16 = AddSIMD(ggor_gm, AddSIMD(gpor_gm, AddSIMD(psigor_gm,   FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(psigD_dKOD01, f0_of_xx0), SubSIMD(FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psigD_dKOD00, FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_24, gp), FusedMulAddSIMD(FDPart3_21, psigD1, FusedMulAddSIMD(FDPart3_4, psigD_dKOD02, FusedMulAddSIMD(FDPart3_12, psigD_dD11, FusedMulAddSIMD(FDPart3_19, psigD_dD22, FusedMulAddSIMD(chip, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(MulSIMD(FDPart3_1, FDPart3_23), MulSIMD(gg, rcoord))), FusedMulSubSIMD(FDPart3_1, FDPart3_9, psigD_dD00)))))))), MulSIMD(MulSIMD(FDPart3_23, FDPart3_Integer_2), MulSIMD(psigD0, rcoord)))))));
           const REAL_SIMD_ARRAY __RHS_exp_17 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psigD_dKOD10, FusedMulAddSIMD(FDPart3_4, psigD_dKOD12, FusedMulAddSIMD(FDPart3_Rational_1_2, psigD_dD01, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(psigD_dKOD11, f0_of_xx0), FusedMulSubSIMD(FDPart3_2, gp_dD1, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psigD1, chip))))))); //Changing from this and on 
           const REAL_SIMD_ARRAY __RHS_exp_18 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psigD_dKOD20, FusedMulAddSIMD(FDPart3_60, psigD_dD02, FusedMulAddSIMD(FDPart3_61, gp_dD2, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(psigD_dKOD21, f0_of_xx0), FusedMulSubSIMD(FDPart3_4, psigD_dKOD22, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psigD2, chip)))))));
           const REAL_SIMD_ARRAY __RHS_exp_19 = AddSIMD(bbor_bm, AddSIMD(bpor_bm, AddSIMD(psibor_bm,   FusedMulAddSIMD(MulSIMD(FDPart3_Rational_1_2, psigD0), MulSIMD(MulSIMD(FDPart3_32, FDPart3_NegativeOne_), MulSIMD(Rprime, gp)), SubSIMD(NegFusedMulAddSIMD(FDPart3_63, MulSIMD(MulSIMD(FDPart3_16, FDPart3_42), MulSIMD(FDPart3_Rational_1_4, Rprime)), SubSIMD(FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(MulSIMD(FDPart3_32, FDPart3_NegativeOne_), MulSIMD(Rprime, fgaugep)), SubSIMD(FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(MulSIMD(FDPart3_1, FDPart3_NegativeOne_), MulSIMD(psifD0, Rprime)), FusedMulAddSIMD(FDPart3_33, MulSIMD(FDPart3_Rational_1_2, fgauge), FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(psibD_dKOD01, f0_of_xx0), FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_24, bp), FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_8, psibD0), FusedMulAddSIMD(FDPart3_64, FDPart3_Rational_1_2, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psibD_dKOD00, FusedMulAddSIMD(FDPart3_33, FDPart3_37, FusedMulAddSIMD(FDPart3_4, psibD_dKOD02, FusedMulAddSIMD(FDPart3_19, psibD_dD22, FusedMulAddSIMD(FDPart3_21, psibD1, FusedMulAddSIMD(bb, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(MulSIMD(FDPart3_1, FDPart3_23), MulSIMD(chip, rcoord))), FusedMulSubSIMD(FDPart3_12, psibD_dD11, psibD_dD00))))))))))))), MulSIMD(MulSIMD(FDPart3_1, FDPart3_29), MulSIMD(FDPart3_Rational_1_4, Rprime)))), MulSIMD(MulSIMD(FDPart3_23, FDPart3_Integer_2), MulSIMD(psibD0, rcoord)))), MulSIMD(MulSIMD(FDPart3_39, FDPart3_63), MulSIMD(FDPart3_Rational_1_4, Rprime)))))));
           const REAL_SIMD_ARRAY __RHS_exp_20 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psibD_dKOD10, FusedMulAddSIMD(FDPart3_4, psibD_dKOD12, FusedMulAddSIMD(FDPart3_Rational_1_2, psibD_dD01, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(psibD_dKOD11, f0_of_xx0), FusedMulSubSIMD(FDPart3_2, bp_dD1, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psibD1, chip)))))));
           const REAL_SIMD_ARRAY __RHS_exp_21 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psibD_dKOD20, FusedMulAddSIMD(FDPart3_60, psibD_dD02, FusedMulAddSIMD(FDPart3_61, bp_dD2, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(psibD_dKOD21, f0_of_xx0), FusedMulSubSIMD(FDPart3_4, psibD_dKOD22, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psibD2, chip)))))));
           const REAL_SIMD_ARRAY __RHS_exp_22 = AddSIMD(uuor_um, AddSIMD(upor_um, AddSIMD(psiuor_um,   SubSIMD(FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_8, psiuD0), FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(psiuD_dKOD01, f0_of_xx0), FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psiuD_dKOD00, FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_24, up), FusedMulAddSIMD(FDPart3_4, psiuD_dKOD02, FusedMulAddSIMD(FDPart3_65, up, FusedMulAddSIMD(FDPart3_21, psiuD1, FusedMulAddSIMD(FDPart3_33, uu, FusedMulAddSIMD(FDPart3_12, psiuD_dD11, FusedMulAddSIMD(FDPart3_19, psiuD_dD22, FusedMulAddSIMD(chip, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(MulSIMD(FDPart3_1, FDPart3_23), MulSIMD(rcoord, uu))), FusedMulSubSIMD(FDPart3_1, FDPart3_43, psiuD_dD00)))))))))))), MulSIMD(MulSIMD(FDPart3_23, FDPart3_Integer_2), MulSIMD(psiuD0, rcoord))))));
           const REAL_SIMD_ARRAY __RHS_exp_23 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psiuD_dKOD10, FusedMulAddSIMD(FDPart3_4, psiuD_dKOD12, FusedMulAddSIMD(FDPart3_Rational_1_2, psiuD_dD01, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(psiuD_dKOD11, f0_of_xx0), FusedMulSubSIMD(FDPart3_2, up_dD1, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psiuD1, chip)))))));
           const REAL_SIMD_ARRAY __RHS_exp_24 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psiuD_dKOD20, FusedMulAddSIMD(FDPart3_60, psiuD_dD02, FusedMulAddSIMD(FDPart3_61, up_dD2, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(psiuD_dKOD21, f0_of_xx0), FusedMulSubSIMD(FDPart3_4, psiuD_dKOD22, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psiuD2, chip)))))));
           const REAL_SIMD_ARRAY __RHS_exp_25 = AddSIMD(fgaugeor_fm, AddSIMD(fgaugepor_fm, AddSIMD(psifor_fm,   NegFusedMulAddSIMD(FDPart3_63, MulSIMD(MulSIMD(FDPart3_16, FDPart3_42), MulSIMD(FDPart3_Rational_1_2, Rprime)), SubSIMD(SubSIMD(FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(psifD_dKOD01, f0_of_xx0), SubSIMD(FusedMulAddSIMD(FDPart3_32, MulSIMD(FDPart3_9, gg), FusedMulAddSIMD(FDPart3_65, MulSIMD(psigD0, gp), FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_24, fgaugep), FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_8, psifD0), FusedMulAddSIMD(FDPart3_65, fgaugep, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psifD_dKOD00, FusedMulAddSIMD(FDPart3_33, fgauge, FusedMulAddSIMD(FDPart3_4, psifD_dKOD02, FusedMulAddSIMD(FDPart3_19, psifD_dD22, FusedMulAddSIMD(FDPart3_21, psifD1, FusedMulAddSIMD(FDPart3_1, FDPart3_51, FusedMulAddSIMD(FDPart3_12, psifD_dD11, FusedMulAddSIMD(chip, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(MulSIMD(FDPart3_1, FDPart3_23), MulSIMD(fgauge, rcoord))), SubSIMD(FDPart3_64, psifD_dD00)))))))))))))), MulSIMD(MulSIMD(FDPart3_1, FDPart3_29), MulSIMD(FDPart3_Rational_1_2, Rprime)))), MulSIMD(MulSIMD(FDPart3_39, FDPart3_63), MulSIMD(FDPart3_Rational_1_2, Rprime))), MulSIMD(MulSIMD(FDPart3_23, FDPart3_Integer_2), MulSIMD(psifD0, rcoord)))))));
           const REAL_SIMD_ARRAY __RHS_exp_26 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psifD_dKOD10, FusedMulAddSIMD(FDPart3_4, psifD_dKOD12, FusedMulAddSIMD(FDPart3_Rational_1_2, psifD_dD01, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(psifD_dKOD11, f0_of_xx0), FusedMulSubSIMD(FDPart3_2, fgaugep_dD1, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psifD1, chip)))))));
           const REAL_SIMD_ARRAY __RHS_exp_27 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psifD_dKOD20, FusedMulAddSIMD(FDPart3_60, psifD_dD02, FusedMulAddSIMD(FDPart3_61, fgaugep_dD2, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, DivSIMD(psifD_dKOD21, f0_of_xx0), FusedMulSubSIMD(FDPart3_4, psifD_dKOD22, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psifD2, chip)))))));
           WriteSIMD(&rhs_gfs[IDX4S(GGGF, i0, i1, i2)], __RHS_exp_0);
           WriteSIMD(&rhs_gfs[IDX4S(GPGF, i0, i1, i2)], __RHS_exp_1);
           WriteSIMD(&rhs_gfs[IDX4S(BBGF, i0, i1, i2)], __RHS_exp_2);
           WriteSIMD(&rhs_gfs[IDX4S(BPGF, i0, i1, i2)], __RHS_exp_3);
           WriteSIMD(&rhs_gfs[IDX4S(UUGF, i0, i1, i2)], __RHS_exp_4);
           WriteSIMD(&rhs_gfs[IDX4S(UPGF, i0, i1, i2)], __RHS_exp_5);
           WriteSIMD(&rhs_gfs[IDX4S(FGAUGEGF, i0, i1, i2)], __RHS_exp_6);
           WriteSIMD(&rhs_gfs[IDX4S(FGAUGEPGF, i0, i1, i2)], __RHS_exp_7);
           WriteSIMD(&rhs_gfs[IDX4S(CHIGF, i0, i1, i2)], __RHS_exp_8);
           WriteSIMD(&rhs_gfs[IDX4S(CHIPGF, i0, i1, i2)], __RHS_exp_9);
           WriteSIMD(&rhs_gfs[IDX4S(CHIPPGF, i0, i1, i2)], __RHS_exp_10);
           WriteSIMD(&rhs_gfs[IDX4S(RCOORDGF, i0, i1, i2)], __RHS_exp_11);
           WriteSIMD(&rhs_gfs[IDX4S(RADIUSGF, i0, i1, i2)], __RHS_exp_12);
           WriteSIMD(&rhs_gfs[IDX4S(RPRIMEGF, i0, i1, i2)], __RHS_exp_13);
           WriteSIMD(&rhs_gfs[IDX4S(SINTHGF, i0, i1, i2)], __RHS_exp_14);
           WriteSIMD(&rhs_gfs[IDX4S(COSTHGF, i0, i1, i2)], __RHS_exp_15);
           WriteSIMD(&rhs_gfs[IDX4S(PSIGD0GF, i0, i1, i2)], __RHS_exp_16);
           WriteSIMD(&rhs_gfs[IDX4S(PSIGD1GF, i0, i1, i2)], __RHS_exp_17);
           WriteSIMD(&rhs_gfs[IDX4S(PSIGD2GF, i0, i1, i2)], __RHS_exp_18);
           WriteSIMD(&rhs_gfs[IDX4S(PSIBD0GF, i0, i1, i2)], __RHS_exp_19);
           WriteSIMD(&rhs_gfs[IDX4S(PSIBD1GF, i0, i1, i2)], __RHS_exp_20);
           WriteSIMD(&rhs_gfs[IDX4S(PSIBD2GF, i0, i1, i2)], __RHS_exp_21);
           WriteSIMD(&rhs_gfs[IDX4S(PSIUD0GF, i0, i1, i2)], __RHS_exp_22);
           WriteSIMD(&rhs_gfs[IDX4S(PSIUD1GF, i0, i1, i2)], __RHS_exp_23);
           WriteSIMD(&rhs_gfs[IDX4S(PSIUD2GF, i0, i1, i2)], __RHS_exp_24);
           WriteSIMD(&rhs_gfs[IDX4S(PSIFD0GF, i0, i1, i2)], __RHS_exp_25);
           WriteSIMD(&rhs_gfs[IDX4S(PSIFD1GF, i0, i1, i2)], __RHS_exp_26);
           WriteSIMD(&rhs_gfs[IDX4S(PSIFD2GF, i0, i1, i2)], __RHS_exp_27);
        }
      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0 += SIMD_width)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
