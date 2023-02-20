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
           const REAL_SIMD_ARRAY FDPart1_1 = MulSIMD(FDPart1_Rational_3_8, bb);
           const REAL_SIMD_ARRAY FDPart1_2 = MulSIMD(FDPart1_Rational_1_2, invdx1);
           const REAL_SIMD_ARRAY FDPart1_3 = MulSIMD(FDPart1_Rational_1_2, invdx2);
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
           const REAL_SIMD_ARRAY bb_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(bb_i0m1_i1_i2, bb_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(bb_i0m2_i1_i2, bb_i0p2_i1_i2), FDPart1_1)));
           const REAL_SIMD_ARRAY bb_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(bb_i0_i1m1_i2, bb_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(bb_i0_i1m2_i2, bb_i0_i1p2_i2), FDPart1_1)));
           const REAL_SIMD_ARRAY bb_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(bb_i0_i1_i2m1, bb_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(bb_i0_i1_i2m2, bb_i0_i1_i2p2), FDPart1_1)));
           const REAL_SIMD_ARRAY bp_dD0 = MulSIMD(FDPart1_0, SubSIMD(bp_i0p1_i1_i2, bp_i0m1_i1_i2));
           const REAL_SIMD_ARRAY bp_dD1 = MulSIMD(FDPart1_2, SubSIMD(bp_i0_i1p1_i2, bp_i0_i1m1_i2));
           const REAL_SIMD_ARRAY bp_dD2 = MulSIMD(FDPart1_3, SubSIMD(bp_i0_i1_i2p1, bp_i0_i1_i2m1));
           const REAL_SIMD_ARRAY bp_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(bp_i0m1_i1_i2, bp_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(bp_i0m2_i1_i2, bp_i0p2_i1_i2), FDPart1_4)));
           const REAL_SIMD_ARRAY bp_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(bp_i0_i1m1_i2, bp_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(bp_i0_i1m2_i2, bp_i0_i1p2_i2), FDPart1_4)));
           const REAL_SIMD_ARRAY bp_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(bp_i0_i1_i2m1, bp_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(bp_i0_i1_i2m2, bp_i0_i1_i2p2), FDPart1_4)));
           const REAL_SIMD_ARRAY fgauge_dD0 = MulSIMD(FDPart1_0, SubSIMD(fgauge_i0p1_i1_i2, fgauge_i0m1_i1_i2));
           const REAL_SIMD_ARRAY fgauge_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(fgauge_i0m1_i1_i2, fgauge_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(fgauge_i0m2_i1_i2, fgauge_i0p2_i1_i2), FDPart1_5)));
           const REAL_SIMD_ARRAY fgauge_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(fgauge_i0_i1m1_i2, fgauge_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(fgauge_i0_i1m2_i2, fgauge_i0_i1p2_i2), FDPart1_5)));
           const REAL_SIMD_ARRAY fgauge_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(fgauge_i0_i1_i2m1, fgauge_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(fgauge_i0_i1_i2m2, fgauge_i0_i1_i2p2), FDPart1_5)));
           const REAL_SIMD_ARRAY fgaugep_dD0 = MulSIMD(FDPart1_0, SubSIMD(fgaugep_i0p1_i1_i2, fgaugep_i0m1_i1_i2));
           const REAL_SIMD_ARRAY fgaugep_dD1 = MulSIMD(FDPart1_2, SubSIMD(fgaugep_i0_i1p1_i2, fgaugep_i0_i1m1_i2));
           const REAL_SIMD_ARRAY fgaugep_dD2 = MulSIMD(FDPart1_3, SubSIMD(fgaugep_i0_i1_i2p1, fgaugep_i0_i1_i2m1));
           const REAL_SIMD_ARRAY fgaugep_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(fgaugep_i0m1_i1_i2, fgaugep_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(fgaugep_i0m2_i1_i2, fgaugep_i0p2_i1_i2), FDPart1_6)));
           const REAL_SIMD_ARRAY fgaugep_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(fgaugep_i0_i1m1_i2, fgaugep_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(fgaugep_i0_i1m2_i2, fgaugep_i0_i1p2_i2), FDPart1_6)));
           const REAL_SIMD_ARRAY fgaugep_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(fgaugep_i0_i1_i2m1, fgaugep_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(fgaugep_i0_i1_i2m2, fgaugep_i0_i1_i2p2), FDPart1_6)));
           const REAL_SIMD_ARRAY gg_dD0 = MulSIMD(FDPart1_0, SubSIMD(gg_i0p1_i1_i2, gg_i0m1_i1_i2));
           const REAL_SIMD_ARRAY gg_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(gg_i0m1_i1_i2, gg_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(gg_i0m2_i1_i2, gg_i0p2_i1_i2), FDPart1_7)));
           const REAL_SIMD_ARRAY gg_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(gg_i0_i1m1_i2, gg_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(gg_i0_i1m2_i2, gg_i0_i1p2_i2), FDPart1_7)));
           const REAL_SIMD_ARRAY gg_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(gg_i0_i1_i2m1, gg_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(gg_i0_i1_i2m2, gg_i0_i1_i2p2), FDPart1_7)));
           const REAL_SIMD_ARRAY gp_dD0 = MulSIMD(FDPart1_0, SubSIMD(gp_i0p1_i1_i2, gp_i0m1_i1_i2));
           const REAL_SIMD_ARRAY gp_dD1 = MulSIMD(FDPart1_2, SubSIMD(gp_i0_i1p1_i2, gp_i0_i1m1_i2));
           const REAL_SIMD_ARRAY gp_dD2 = MulSIMD(FDPart1_3, SubSIMD(gp_i0_i1_i2p1, gp_i0_i1_i2m1));
           const REAL_SIMD_ARRAY gp_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(gp_i0m1_i1_i2, gp_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(gp_i0m2_i1_i2, gp_i0p2_i1_i2), FDPart1_8)));
           const REAL_SIMD_ARRAY gp_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(gp_i0_i1m1_i2, gp_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(gp_i0_i1m2_i2, gp_i0_i1p2_i2), FDPart1_8)));
           const REAL_SIMD_ARRAY gp_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(gp_i0_i1_i2m1, gp_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(gp_i0_i1_i2m2, gp_i0_i1_i2p2), FDPart1_8)));
           const REAL_SIMD_ARRAY psibD_dD00 = MulSIMD(FDPart1_0, SubSIMD(psibD0_i0p1_i1_i2, psibD0_i0m1_i1_i2));
           const REAL_SIMD_ARRAY psibD_dD01 = MulSIMD(FDPart1_2, SubSIMD(psibD0_i0_i1p1_i2, psibD0_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psibD_dD02 = MulSIMD(FDPart1_3, SubSIMD(psibD0_i0_i1_i2p1, psibD0_i0_i1_i2m1));
           const REAL_SIMD_ARRAY psibD_dD11 = MulSIMD(FDPart1_2, SubSIMD(psibD1_i0_i1p1_i2, psibD1_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psibD_dD22 = MulSIMD(FDPart1_3, SubSIMD(psibD2_i0_i1_i2p1, psibD2_i0_i1_i2m1));
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
           const REAL_SIMD_ARRAY psifD_dD01 = MulSIMD(FDPart1_2, SubSIMD(psifD0_i0_i1p1_i2, psifD0_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psifD_dD02 = MulSIMD(FDPart1_3, SubSIMD(psifD0_i0_i1_i2p1, psifD0_i0_i1_i2m1));
           const REAL_SIMD_ARRAY psifD_dD11 = MulSIMD(FDPart1_2, SubSIMD(psifD1_i0_i1p1_i2, psifD1_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psifD_dD22 = MulSIMD(FDPart1_3, SubSIMD(psifD2_i0_i1_i2p1, psifD2_i0_i1_i2m1));
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
           const REAL_SIMD_ARRAY psigD_dD01 = MulSIMD(FDPart1_2, SubSIMD(psigD0_i0_i1p1_i2, psigD0_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psigD_dD02 = MulSIMD(FDPart1_3, SubSIMD(psigD0_i0_i1_i2p1, psigD0_i0_i1_i2m1));
           const REAL_SIMD_ARRAY psigD_dD11 = MulSIMD(FDPart1_2, SubSIMD(psigD1_i0_i1p1_i2, psigD1_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psigD_dD22 = MulSIMD(FDPart1_3, SubSIMD(psigD2_i0_i1_i2p1, psigD2_i0_i1_i2m1));
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
           const REAL_SIMD_ARRAY psiuD_dD01 = MulSIMD(FDPart1_2, SubSIMD(psiuD0_i0_i1p1_i2, psiuD0_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psiuD_dD02 = MulSIMD(FDPart1_3, SubSIMD(psiuD0_i0_i1_i2p1, psiuD0_i0_i1_i2m1));
           const REAL_SIMD_ARRAY psiuD_dD11 = MulSIMD(FDPart1_2, SubSIMD(psiuD1_i0_i1p1_i2, psiuD1_i0_i1m1_i2));
           const REAL_SIMD_ARRAY psiuD_dD22 = MulSIMD(FDPart1_3, SubSIMD(psiuD2_i0_i1_i2p1, psiuD2_i0_i1_i2m1));
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
           const REAL_SIMD_ARRAY up_dD1 = MulSIMD(FDPart1_2, SubSIMD(up_i0_i1p1_i2, up_i0_i1m1_i2));
           const REAL_SIMD_ARRAY up_dD2 = MulSIMD(FDPart1_3, SubSIMD(up_i0_i1_i2p1, up_i0_i1_i2m1));
           const REAL_SIMD_ARRAY up_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(up_i0m1_i1_i2, up_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(up_i0m2_i1_i2, up_i0p2_i1_i2), FDPart1_21)));
           const REAL_SIMD_ARRAY up_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(up_i0_i1m1_i2, up_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(up_i0_i1m2_i2, up_i0_i1p2_i2), FDPart1_21)));
           const REAL_SIMD_ARRAY up_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(up_i0_i1_i2m1, up_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(up_i0_i1_i2m2, up_i0_i1_i2p2), FDPart1_21)));
           const REAL_SIMD_ARRAY uu_dD0 = MulSIMD(FDPart1_0, SubSIMD(uu_i0p1_i1_i2, uu_i0m1_i1_i2));
           const REAL_SIMD_ARRAY uu_dKOD0 = MulSIMD(invdx0, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(uu_i0m1_i1_i2, uu_i0p1_i1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(uu_i0m2_i1_i2, uu_i0p2_i1_i2), FDPart1_22)));
           const REAL_SIMD_ARRAY uu_dKOD1 = MulSIMD(invdx1, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(uu_i0_i1m1_i2, uu_i0_i1p1_i2), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(uu_i0_i1m2_i2, uu_i0_i1p2_i2), FDPart1_22)));
           const REAL_SIMD_ARRAY uu_dKOD2 = MulSIMD(invdx2, FusedMulSubSIMD(FDPart1_Rational_1_4, AddSIMD(uu_i0_i1_i2m1, uu_i0_i1_i2p1), FusedMulAddSIMD(FDPart1_Rational_1_16, AddSIMD(uu_i0_i1_i2m2, uu_i0_i1_i2p2), FDPart1_22)));
           /*
            * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
            */
           /*
            *  Original SymPy expressions:
            *  "[const REAL_SIMD_ARRAY __RHS_exp_0 = 5764607523034235*gg_dKOD0/288230376151711744 + 5764607523034235*gg_dKOD1/288230376151711744 + 5764607523034235*gg_dKOD2/288230376151711744 + psigD0/2 - chip*gg/(2*chi) + gp/(2*chi),
            *    const REAL_SIMD_ARRAY __RHS_exp_1 = 5764607523034235*gp_dKOD0/288230376151711744 + 5764607523034235*gp_dKOD1/288230376151711744 + 5764607523034235*gp_dKOD2/288230376151711744 + (Rprime*chip*psigD0 - Rprime*chipp*gg + Rprime*chip**2*gg/chi - Rprime*chip*gp/chi - 2*chi*psigD0*rcoord/(1 - rcoord**2) - 2*chip*gg*rcoord/(1 - rcoord**2) - chip*gg_dD0 - chip*psigD0/2 + 2*gp*rcoord/(1 - rcoord**2) + gp_dD0 + chip**2*gg/(2*chi) - chip*gp/(2*chi) + Rprime*chi*costh*psigD1/(Radius**2*sinth) + Rprime*chi*psigD_dD11/Radius**2 + Rprime*chi*psigD_dD22/(Radius**2*sinth))/(2*Rprime - 1),
            *    const REAL_SIMD_ARRAY __RHS_exp_2 = -bb*chip/(2*chi) + 5764607523034235*bb_dKOD0/288230376151711744 + 5764607523034235*bb_dKOD1/288230376151711744 + 5764607523034235*bb_dKOD2/288230376151711744 + bp/(2*chi) + psibD0/2,
            *    const REAL_SIMD_ARRAY __RHS_exp_3 = 5764607523034235*bp_dKOD0/288230376151711744 + 5764607523034235*bp_dKOD1/288230376151711744 + 5764607523034235*bp_dKOD2/288230376151711744 + (-Rprime*bb*chipp + Rprime*bb*chip**2/chi - Rprime*bp*chip/chi + Rprime*chip*psibD0 - Rprime*psifD0/2 - Rprime*psigD0**2/4 + Rprime*chip*fgauge/(2*chi) + Rprime*chip*gg*psigD0/(2*chi) - Rprime*fgaugep/(2*chi) - Rprime*gp*psigD0/(2*chi) - Rprime*chip**2*gg**2/(4*chi**2) + Rprime*chip*gg*gp/(2*chi**2) - Rprime*gp**2/(4*chi**2) - 2*bb*chip*rcoord/(1 - rcoord**2) + bb*chip**2/(2*chi) - bb_dD0*chip + 2*bp*rcoord/(1 - rcoord**2) - bp*chip/(2*chi) + bp_dD0 - 2*chi*psibD0*rcoord/(1 - rcoord**2) - chip*psibD0/2 + Rprime*chi*costh*psibD1/(Radius**2*sinth) + Rprime*chi*psibD_dD11/Radius**2 + Rprime*chi*psibD_dD22/(Radius**2*sinth))/(2*Rprime - 1),
            *    const REAL_SIMD_ARRAY __RHS_exp_4 = psiuD0/2 + 5764607523034235*uu_dKOD0/288230376151711744 + 5764607523034235*uu_dKOD1/288230376151711744 + 5764607523034235*uu_dKOD2/288230376151711744 - chip*uu/(2*chi) + up/(2*chi),
            *    const REAL_SIMD_ARRAY __RHS_exp_5 = 5764607523034235*up_dKOD0/288230376151711744 + 5764607523034235*up_dKOD1/288230376151711744 + 5764607523034235*up_dKOD2/288230376151711744 + (Rprime*chip*psiuD0 - Rprime*chipp*uu - Rprime*psiuD0 + Rprime*chip**2*uu/chi - Rprime*chip*up/chi + Rprime*chip*uu/chi - Rprime*up/chi - 2*chi*psiuD0*rcoord/(1 - rcoord**2) - chip*psiuD0/2 - 2*chip*rcoord*uu/(1 - rcoord**2) - chip*uu_dD0 + 2*rcoord*up/(1 - rcoord**2) + up_dD0 + chip**2*uu/(2*chi) - chip*up/(2*chi) + Rprime*chi*costh*psiuD1/(Radius**2*sinth) + Rprime*chi*psiuD_dD11/Radius**2 + Rprime*chi*psiuD_dD22/(Radius**2*sinth))/(2*Rprime - 1),
            *    const REAL_SIMD_ARRAY __RHS_exp_6 = 5764607523034235*fgauge_dKOD0/288230376151711744 + 5764607523034235*fgauge_dKOD1/288230376151711744 + 5764607523034235*fgauge_dKOD2/288230376151711744 + psifD0/2 - chip*fgauge/(2*chi) + fgaugep/(2*chi),
            *    const REAL_SIMD_ARRAY __RHS_exp_7 = 5764607523034235*fgaugep_dKOD0/288230376151711744 + 5764607523034235*fgaugep_dKOD1/288230376151711744 + 5764607523034235*fgaugep_dKOD2/288230376151711744 + (Rprime*chip*psifD0 - Rprime*chipp*fgauge - Rprime*psifD0 - Rprime*psigD0**2/2 + Rprime*chip**2*fgauge/chi + Rprime*chip*fgauge/chi - Rprime*chip*fgaugep/chi + Rprime*chip*gg*psigD0/chi - Rprime*fgaugep/chi - Rprime*gp*psigD0/chi - Rprime*chip**2*gg**2/(2*chi**2) + Rprime*chip*gg*gp/chi**2 - Rprime*gp**2/(2*chi**2) - 2*chi*psifD0*rcoord/(1 - rcoord**2) - 2*chip*fgauge*rcoord/(1 - rcoord**2) - chip*fgauge_dD0 - chip*psifD0/2 + 2*fgaugep*rcoord/(1 - rcoord**2) + fgaugep_dD0 + chip**2*fgauge/(2*chi) - chip*fgaugep/(2*chi) + Rprime*chi*costh*psifD1/(Radius**2*sinth) + Rprime*chi*psifD_dD11/Radius**2 + Rprime*chi*psifD_dD22/(Radius**2*sinth))/(2*Rprime - 1),
            *    const REAL_SIMD_ARRAY __RHS_exp_8 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_9 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_10 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_11 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_12 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_13 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_14 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_15 = 0,
            *    const REAL_SIMD_ARRAY __RHS_exp_16 = Rprime*chip*psigD0/chi - 2*psigD0*rcoord/(1 - rcoord**2) - psigD_dD00 + 5764607523034235*psigD_dKOD00/288230376151711744 + 5764607523034235*psigD_dKOD01/288230376151711744 + 5764607523034235*psigD_dKOD02/288230376151711744 - 2*chip*gg*rcoord/(chi*(1 - rcoord**2)) + 2*gp*rcoord/(chi*(1 - rcoord**2)) + Rprime*costh*psigD1/(Radius**2*sinth) + Rprime*psigD_dD11/Radius**2 + Rprime*psigD_dD22/(Radius**2*sinth),
            *    const REAL_SIMD_ARRAY __RHS_exp_17 = psigD_dD01/2 + 5764607523034235*psigD_dKOD10/288230376151711744 + 5764607523034235*psigD_dKOD11/288230376151711744 + 5764607523034235*psigD_dKOD12/288230376151711744 - chip*psigD1/(2*chi) + gp_dD1/(2*chi),
            *    const REAL_SIMD_ARRAY __RHS_exp_18 = psigD_dD02/(2*sinth) + 5764607523034235*psigD_dKOD20/288230376151711744 + 5764607523034235*psigD_dKOD21/288230376151711744 + 5764607523034235*psigD_dKOD22/288230376151711744 - chip*psigD2/(2*chi) + gp_dD2/(2*chi*sinth),
            *    const REAL_SIMD_ARRAY __RHS_exp_19 = Rprime*chip*psibD0/chi - Rprime*psifD0/(2*chi) - Rprime*psigD0**2/(4*chi) + Rprime*chip*fgauge/(2*chi**2) + Rprime*chip*gg*psigD0/(2*chi**2) - Rprime*fgaugep/(2*chi**2) - Rprime*gp*psigD0/(2*chi**2) - Rprime*chip**2*gg**2/(4*chi**3) + Rprime*chip*gg*gp/(2*chi**3) - Rprime*gp**2/(4*chi**3) - 2*bb*chip*rcoord/(chi*(1 - rcoord**2)) + 2*bp*rcoord/(chi*(1 - rcoord**2)) - 2*psibD0*rcoord/(1 - rcoord**2) - psibD_dD00 + 5764607523034235*psibD_dKOD00/288230376151711744 + 5764607523034235*psibD_dKOD01/288230376151711744 + 5764607523034235*psibD_dKOD02/288230376151711744 + Rprime*costh*psibD1/(Radius**2*sinth) + Rprime*psibD_dD11/Radius**2 + Rprime*psibD_dD22/(Radius**2*sinth),
            *    const REAL_SIMD_ARRAY __RHS_exp_20 = bp_dD1/(2*chi) + psibD_dD01/2 + 5764607523034235*psibD_dKOD10/288230376151711744 + 5764607523034235*psibD_dKOD11/288230376151711744 + 5764607523034235*psibD_dKOD12/288230376151711744 - chip*psibD1/(2*chi),
            *    const REAL_SIMD_ARRAY __RHS_exp_21 = bp_dD2/(2*chi*sinth) + psibD_dD02/(2*sinth) + 5764607523034235*psibD_dKOD20/288230376151711744 + 5764607523034235*psibD_dKOD21/288230376151711744 + 5764607523034235*psibD_dKOD22/288230376151711744 - chip*psibD2/(2*chi),
            *    const REAL_SIMD_ARRAY __RHS_exp_22 = Rprime*chip*psiuD0/chi - Rprime*psiuD0/chi + Rprime*chip*uu/chi**2 - Rprime*up/chi**2 - 2*psiuD0*rcoord/(1 - rcoord**2) - psiuD_dD00 + 5764607523034235*psiuD_dKOD00/288230376151711744 + 5764607523034235*psiuD_dKOD01/288230376151711744 + 5764607523034235*psiuD_dKOD02/288230376151711744 - 2*chip*rcoord*uu/(chi*(1 - rcoord**2)) + 2*rcoord*up/(chi*(1 - rcoord**2)) + Rprime*costh*psiuD1/(Radius**2*sinth) + Rprime*psiuD_dD11/Radius**2 + Rprime*psiuD_dD22/(Radius**2*sinth),
            *    const REAL_SIMD_ARRAY __RHS_exp_23 = psiuD_dD01/2 + 5764607523034235*psiuD_dKOD10/288230376151711744 + 5764607523034235*psiuD_dKOD11/288230376151711744 + 5764607523034235*psiuD_dKOD12/288230376151711744 - chip*psiuD1/(2*chi) + up_dD1/(2*chi),
            *    const REAL_SIMD_ARRAY __RHS_exp_24 = psiuD_dD02/(2*sinth) + 5764607523034235*psiuD_dKOD20/288230376151711744 + 5764607523034235*psiuD_dKOD21/288230376151711744 + 5764607523034235*psiuD_dKOD22/288230376151711744 - chip*psiuD2/(2*chi) + up_dD2/(2*chi*sinth),
            *    const REAL_SIMD_ARRAY __RHS_exp_25 = Rprime*chip*psifD0/chi - Rprime*psifD0/chi - Rprime*psigD0**2/(2*chi) + Rprime*chip*fgauge/chi**2 + Rprime*chip*gg*psigD0/chi**2 - Rprime*fgaugep/chi**2 - Rprime*gp*psigD0/chi**2 - Rprime*chip**2*gg**2/(2*chi**3) + Rprime*chip*gg*gp/chi**3 - Rprime*gp**2/(2*chi**3) - 2*psifD0*rcoord/(1 - rcoord**2) - psifD_dD00 + 5764607523034235*psifD_dKOD00/288230376151711744 + 5764607523034235*psifD_dKOD01/288230376151711744 + 5764607523034235*psifD_dKOD02/288230376151711744 - 2*chip*fgauge*rcoord/(chi*(1 - rcoord**2)) + 2*fgaugep*rcoord/(chi*(1 - rcoord**2)) + Rprime*costh*psifD1/(Radius**2*sinth) + Rprime*psifD_dD11/Radius**2 + Rprime*psifD_dD22/(Radius**2*sinth),
            *    const REAL_SIMD_ARRAY __RHS_exp_26 = psifD_dD01/2 + 5764607523034235*psifD_dKOD10/288230376151711744 + 5764607523034235*psifD_dKOD11/288230376151711744 + 5764607523034235*psifD_dKOD12/288230376151711744 - chip*psifD1/(2*chi) + fgaugep_dD1/(2*chi),
            *    const REAL_SIMD_ARRAY __RHS_exp_27 = psifD_dD02/(2*sinth) + 5764607523034235*psifD_dKOD20/288230376151711744 + 5764607523034235*psifD_dKOD21/288230376151711744 + 5764607523034235*psifD_dKOD22/288230376151711744 - chip*psifD2/(2*chi) + fgaugep_dD2/(2*chi*sinth)]"
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
           
           const double tmpFDPart3_Rational_5764607523034235_288230376151711744 = 5764607523034235.0/2.8823037615171174e+17;
           const REAL_SIMD_ARRAY FDPart3_Rational_5764607523034235_288230376151711744 = ConstSIMD(tmpFDPart3_Rational_5764607523034235_288230376151711744);
           
           const REAL_SIMD_ARRAY FDPart3_1 = DivSIMD(FDPart3_Integer_1, chi);
           const REAL_SIMD_ARRAY FDPart3_2 = MulSIMD(FDPart3_1, FDPart3_Rational_1_2);
           const REAL_SIMD_ARRAY FDPart3_4 = DivSIMD(FDPart3_Integer_1, FusedMulAddSIMD(FDPart3_Integer_2, Rprime, FDPart3_NegativeOne_));
           const REAL_SIMD_ARRAY FDPart3_6 = MulSIMD(Rprime, chip);
           const REAL_SIMD_ARRAY FDPart3_7 = MulSIMD(FDPart3_6, psigD0);
           const REAL_SIMD_ARRAY FDPart3_9 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(Rprime, chipp));
           const REAL_SIMD_ARRAY FDPart3_10 = DivSIMD(Rprime, MulSIMD(Radius, Radius));
           const REAL_SIMD_ARRAY FDPart3_12 = MulSIMD(FDPart3_1, Rprime);
           const REAL_SIMD_ARRAY FDPart3_13 = MulSIMD(FDPart3_12, MulSIMD(FDPart3_NegativeOne_, chip));
           const REAL_SIMD_ARRAY FDPart3_14 = MulSIMD(chip, chip);
           const REAL_SIMD_ARRAY FDPart3_16 = DivSIMD(FDPart3_Integer_1, sinth);
           const REAL_SIMD_ARRAY FDPart3_17 = MulSIMD(FDPart3_10, FDPart3_16);
           const REAL_SIMD_ARRAY FDPart3_19 = MulSIMD(FDPart3_17, costh);
           const REAL_SIMD_ARRAY FDPart3_21 = DivSIMD(FDPart3_Integer_1, NegFusedMulAddSIMD(rcoord, rcoord, FDPart3_Integer_1));
           const REAL_SIMD_ARRAY FDPart3_22 = MulSIMD(FDPart3_21, MulSIMD(FDPart3_Integer_2, rcoord));
           const REAL_SIMD_ARRAY FDPart3_27 = MulSIMD(psigD0, psigD0);
           const REAL_SIMD_ARRAY FDPart3_30 = DivSIMD(FDPart3_Integer_1, MulSIMD(chi, chi));
           const REAL_SIMD_ARRAY FDPart3_31 = MulSIMD(FDPart3_30, FDPart3_6);
           const REAL_SIMD_ARRAY FDPart3_33 = MulSIMD(FDPart3_31, MulSIMD(gg, gp));
           const REAL_SIMD_ARRAY FDPart3_34 = MulSIMD(FDPart3_1, FDPart3_6);
           const REAL_SIMD_ARRAY FDPart3_35 = MulSIMD(FDPart3_Rational_1_2, MulSIMD(psigD0, gg));
           const REAL_SIMD_ARRAY FDPart3_37 = MulSIMD(gp, gp);
           const REAL_SIMD_ARRAY FDPart3_40 = MulSIMD(gg, gg);
           const REAL_SIMD_ARRAY FDPart3_41 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(psiuD0, Rprime));
           const REAL_SIMD_ARRAY FDPart3_43 = MulSIMD(FDPart3_1, MulSIMD(FDPart3_NegativeOne_, Rprime));
           const REAL_SIMD_ARRAY FDPart3_49 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(psifD0, Rprime));
           const REAL_SIMD_ARRAY FDPart3_58 = MulSIMD(FDPart3_16, FDPart3_Rational_1_2);
           const REAL_SIMD_ARRAY FDPart3_59 = MulSIMD(FDPart3_16, FDPart3_2);
           const REAL_SIMD_ARRAY FDPart3_61 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(chi, chi), chi));
           const REAL_SIMD_ARRAY FDPart3_62 = MulSIMD(MulSIMD(FDPart3_6, FDPart3_61), MulSIMD(gg, gp));
           const REAL_SIMD_ARRAY FDPart3_63 = MulSIMD(FDPart3_30, MulSIMD(FDPart3_NegativeOne_, Rprime));
           const REAL_SIMD_ARRAY __RHS_exp_0 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, gg_dKOD1, FusedMulAddSIMD(FDPart3_Rational_1_2, psigD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, gg_dKOD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, gg_dKOD2, FusedMulSubSIMD(FDPart3_2, gp, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(chip, gg)))))));
           const REAL_SIMD_ARRAY __RHS_exp_1 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, gp_dKOD1, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, gp_dKOD2, FusedMulAddSIMD(FDPart3_4, FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(MulSIMD(FDPart3_1, FDPart3_NegativeOne_), MulSIMD(chip, gp)), NegFusedMulAddSIMD(psigD0, MulSIMD(MulSIMD(FDPart3_21, FDPart3_Integer_2), MulSIMD(chi, rcoord)), FusedMulAddSIMD(FDPart3_19, MulSIMD(psigD1, chi), NegFusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(psigD0, chip), FusedMulAddSIMD(FDPart3_14, MulSIMD(FDPart3_2, gg), FusedMulAddSIMD(FDPart3_17, MulSIMD(psigD_dD22, chi), FusedMulAddSIMD(FDPart3_10, MulSIMD(psigD_dD11, chi), FusedMulAddSIMD(FDPart3_12, MulSIMD(FDPart3_14, gg), FusedMulAddSIMD(FDPart3_22, gp, FusedMulAddSIMD(FDPart3_9, gg, AddSIMD(FusedMulAddSIMD(FDPart3_13, gp, gp_dD0), NegFusedMulAddSIMD(chip, MulSIMD(MulSIMD(FDPart3_21, FDPart3_Integer_2), MulSIMD(gg, rcoord)), NegFusedMulAddSIMD(gg_dD0, chip, FDPart3_7))))))))))))), MulSIMD(FDPart3_Rational_5764607523034235_288230376151711744, gp_dKOD0))));
           const REAL_SIMD_ARRAY __RHS_exp_2 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, bb_dKOD1, FusedMulAddSIMD(FDPart3_Rational_1_2, psibD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, bb_dKOD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, bb_dKOD2, FusedMulSubSIMD(FDPart3_2, bp, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(bb, chip)))))));
           const REAL_SIMD_ARRAY __RHS_exp_3 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, bp_dKOD1, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, bp_dKOD2, FusedMulAddSIMD(FDPart3_4, NegFusedMulAddSIMD(psibD0, MulSIMD(MulSIMD(FDPart3_21, FDPart3_Integer_2), MulSIMD(chi, rcoord)), FusedMulAddSIMD(MulSIMD(FDPart3_Rational_1_2, psigD0), MulSIMD(MulSIMD(FDPart3_1, FDPart3_NegativeOne_), MulSIMD(Rprime, gp)), NegFusedMulAddSIMD(FDPart3_40, MulSIMD(MulSIMD(FDPart3_14, FDPart3_30), MulSIMD(FDPart3_Rational_1_4, Rprime)), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(MulSIMD(FDPart3_1, FDPart3_NegativeOne_), MulSIMD(bp, chip)), SubSIMD(NegFusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(psifD0, Rprime), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(MulSIMD(FDPart3_1, FDPart3_NegativeOne_), MulSIMD(Rprime, fgaugep)), NegFusedMulAddSIMD(FDPart3_27, MulSIMD(FDPart3_Rational_1_4, Rprime), NegFusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(psibD0, chip), FusedMulAddSIMD(FDPart3_19, MulSIMD(psibD1, chi), FusedMulAddSIMD(FDPart3_2, MulSIMD(FDPart3_6, fgauge), FusedMulAddSIMD(FDPart3_14, MulSIMD(FDPart3_2, bb), FusedMulAddSIMD(FDPart3_17, MulSIMD(psibD_dD22, chi), FusedMulAddSIMD(FDPart3_10, MulSIMD(psibD_dD11, chi), FusedMulAddSIMD(FDPart3_12, MulSIMD(FDPart3_14, bb), FusedMulAddSIMD(FDPart3_6, psibD0, FusedMulAddSIMD(FDPart3_9, bb, FusedMulAddSIMD(FDPart3_33, FDPart3_Rational_1_2, FusedMulAddSIMD(FDPart3_34, FDPart3_35, FusedMulAddSIMD(FDPart3_13, bp, FusedMulAddSIMD(FDPart3_22, bp, NegFusedMulAddSIMD(bb, MulSIMD(MulSIMD(FDPart3_21, FDPart3_Integer_2), MulSIMD(chip, rcoord)), NegFusedMulAddSIMD(bb_dD0, chip, bp_dD0)))))))))))))))))), MulSIMD(MulSIMD(FDPart3_30, FDPart3_37), MulSIMD(FDPart3_Rational_1_4, Rprime))))))), MulSIMD(FDPart3_Rational_5764607523034235_288230376151711744, bp_dKOD0))));
           const REAL_SIMD_ARRAY __RHS_exp_4 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, uu_dKOD1, FusedMulAddSIMD(FDPart3_Rational_1_2, psiuD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, uu_dKOD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, uu_dKOD2, FusedMulSubSIMD(FDPart3_2, up, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(chip, uu)))))));
           const REAL_SIMD_ARRAY __RHS_exp_5 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, up_dKOD1, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, up_dKOD2, FusedMulAddSIMD(FDPart3_4, NegFusedMulAddSIMD(psiuD0, MulSIMD(MulSIMD(FDPart3_21, FDPart3_Integer_2), MulSIMD(chi, rcoord)), NegFusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(psiuD0, chip), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(MulSIMD(FDPart3_1, FDPart3_NegativeOne_), MulSIMD(chip, up)), FusedMulAddSIMD(FDPart3_17, MulSIMD(psiuD_dD22, chi), FusedMulAddSIMD(FDPart3_19, MulSIMD(psiuD1, chi), FusedMulAddSIMD(FDPart3_12, MulSIMD(FDPart3_14, uu), FusedMulAddSIMD(FDPart3_14, MulSIMD(FDPart3_2, uu), FusedMulAddSIMD(FDPart3_9, uu, FusedMulAddSIMD(FDPart3_10, MulSIMD(psiuD_dD11, chi), FusedMulAddSIMD(FDPart3_43, up, FusedMulAddSIMD(FDPart3_6, psiuD0, FusedMulAddSIMD(FDPart3_22, up, FusedMulAddSIMD(FDPart3_34, uu, AddSIMD(FusedMulAddSIMD(FDPart3_13, up, up_dD0), NegFusedMulAddSIMD(chip, MulSIMD(MulSIMD(FDPart3_21, FDPart3_Integer_2), MulSIMD(rcoord, uu)), NegFusedMulAddSIMD(uu_dD0, chip, FDPart3_41)))))))))))))))), MulSIMD(FDPart3_Rational_5764607523034235_288230376151711744, up_dKOD0))));
           const REAL_SIMD_ARRAY __RHS_exp_6 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, fgauge_dKOD1, FusedMulAddSIMD(FDPart3_Rational_1_2, psifD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, fgauge_dKOD0, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, fgauge_dKOD2, FusedMulSubSIMD(FDPart3_2, fgaugep, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(chip, fgauge)))))));
           const REAL_SIMD_ARRAY __RHS_exp_7 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, fgaugep_dKOD1, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, fgaugep_dKOD2, FusedMulAddSIMD(FDPart3_4, NegFusedMulAddSIMD(psifD0, MulSIMD(MulSIMD(FDPart3_21, FDPart3_Integer_2), MulSIMD(chi, rcoord)), SubSIMD(NegFusedMulAddSIMD(FDPart3_40, MulSIMD(MulSIMD(FDPart3_14, FDPart3_30), MulSIMD(FDPart3_Rational_1_2, Rprime)), NegFusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(psifD0, chip), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(MulSIMD(FDPart3_1, FDPart3_NegativeOne_), MulSIMD(chip, fgaugep)), FusedMulAddSIMD(FDPart3_43, MulSIMD(psigD0, gp), NegFusedMulAddSIMD(FDPart3_27, MulSIMD(FDPart3_Rational_1_2, Rprime), FusedMulAddSIMD(FDPart3_17, MulSIMD(psifD_dD22, chi), FusedMulAddSIMD(FDPart3_19, MulSIMD(psifD1, chi), FusedMulAddSIMD(FDPart3_12, MulSIMD(FDPart3_14, fgauge), FusedMulAddSIMD(FDPart3_14, MulSIMD(FDPart3_2, fgauge), FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_7, gg), FusedMulAddSIMD(FDPart3_10, MulSIMD(psifD_dD11, chi), FusedMulAddSIMD(FDPart3_6, psifD0, FusedMulAddSIMD(FDPart3_9, fgauge, FusedMulAddSIMD(FDPart3_34, fgauge, FusedMulAddSIMD(FDPart3_43, fgaugep, FusedMulAddSIMD(FDPart3_13, fgaugep, FusedMulAddSIMD(FDPart3_22, fgaugep, AddSIMD(AddSIMD(FDPart3_49, fgaugep_dD0), NegFusedMulAddSIMD(chip, MulSIMD(MulSIMD(FDPart3_21, FDPart3_Integer_2), MulSIMD(fgauge, rcoord)), NegFusedMulAddSIMD(fgauge_dD0, chip, FDPart3_33)))))))))))))))))))), MulSIMD(MulSIMD(FDPart3_30, FDPart3_37), MulSIMD(FDPart3_Rational_1_2, Rprime)))), MulSIMD(FDPart3_Rational_5764607523034235_288230376151711744, fgaugep_dKOD0))));
           const REAL_SIMD_ARRAY __RHS_exp_8 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_9 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_10 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_11 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_12 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_13 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_14 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_15 = FDPart3_Integer_0;
           const REAL_SIMD_ARRAY __RHS_exp_16 = FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_22, gp), SubSIMD(FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psigD_dKOD01, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psigD_dKOD02, FusedMulAddSIMD(FDPart3_19, psigD1, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psigD_dKOD00, FusedMulAddSIMD(FDPart3_10, psigD_dD11, FusedMulAddSIMD(FDPart3_17, psigD_dD22, FusedMulAddSIMD(chip, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(MulSIMD(FDPart3_1, FDPart3_21), MulSIMD(gg, rcoord))), FusedMulSubSIMD(FDPart3_1, FDPart3_7, psigD_dD00)))))))), MulSIMD(MulSIMD(FDPart3_21, FDPart3_Integer_2), MulSIMD(psigD0, rcoord))));
           const REAL_SIMD_ARRAY __RHS_exp_17 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psigD_dKOD11, FusedMulAddSIMD(FDPart3_Rational_1_2, psigD_dD01, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psigD_dKOD10, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psigD_dKOD12, FusedMulSubSIMD(FDPart3_2, gp_dD1, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psigD1, chip)))))));
           const REAL_SIMD_ARRAY __RHS_exp_18 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psigD_dKOD21, FusedMulAddSIMD(FDPart3_59, gp_dD2, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psigD_dKOD20, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psigD_dKOD22, FusedMulSubSIMD(FDPart3_58, psigD_dD02, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psigD2, chip)))))));
           const REAL_SIMD_ARRAY __RHS_exp_19 = FusedMulAddSIMD(MulSIMD(FDPart3_Rational_1_2, psigD0), MulSIMD(MulSIMD(FDPart3_30, FDPart3_NegativeOne_), MulSIMD(Rprime, gp)), SubSIMD(NegFusedMulAddSIMD(FDPart3_61, MulSIMD(MulSIMD(FDPart3_14, FDPart3_40), MulSIMD(FDPart3_Rational_1_4, Rprime)), SubSIMD(FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(MulSIMD(FDPart3_30, FDPart3_NegativeOne_), MulSIMD(Rprime, fgaugep)), SubSIMD(FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(MulSIMD(FDPart3_1, FDPart3_NegativeOne_), MulSIMD(psifD0, Rprime)), FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_6, psibD0), FusedMulAddSIMD(FDPart3_31, MulSIMD(FDPart3_Rational_1_2, fgauge), FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psibD_dKOD02, FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_22, bp), FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psibD_dKOD00, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psibD_dKOD01, FusedMulAddSIMD(FDPart3_31, FDPart3_35, FusedMulAddSIMD(FDPart3_62, FDPart3_Rational_1_2, FusedMulAddSIMD(FDPart3_17, psibD_dD22, FusedMulAddSIMD(FDPart3_19, psibD1, FusedMulAddSIMD(bb, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(MulSIMD(FDPart3_1, FDPart3_21), MulSIMD(chip, rcoord))), FusedMulSubSIMD(FDPart3_10, psibD_dD11, psibD_dD00))))))))))))), MulSIMD(MulSIMD(FDPart3_1, FDPart3_27), MulSIMD(FDPart3_Rational_1_4, Rprime)))), MulSIMD(MulSIMD(FDPart3_21, FDPart3_Integer_2), MulSIMD(psibD0, rcoord)))), MulSIMD(MulSIMD(FDPart3_37, FDPart3_61), MulSIMD(FDPart3_Rational_1_4, Rprime))));
           const REAL_SIMD_ARRAY __RHS_exp_20 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psibD_dKOD11, FusedMulAddSIMD(FDPart3_Rational_1_2, psibD_dD01, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psibD_dKOD10, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psibD_dKOD12, FusedMulSubSIMD(FDPart3_2, bp_dD1, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psibD1, chip)))))));
           const REAL_SIMD_ARRAY __RHS_exp_21 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psibD_dKOD21, FusedMulAddSIMD(FDPart3_59, bp_dD2, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psibD_dKOD20, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psibD_dKOD22, FusedMulSubSIMD(FDPart3_58, psibD_dD02, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psibD2, chip)))))));
           const REAL_SIMD_ARRAY __RHS_exp_22 = SubSIMD(FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_22, up), FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_6, psiuD0), FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psiuD_dKOD01, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psiuD_dKOD02, FusedMulAddSIMD(FDPart3_63, up, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psiuD_dKOD00, FusedMulAddSIMD(FDPart3_19, psiuD1, FusedMulAddSIMD(FDPart3_31, uu, FusedMulAddSIMD(FDPart3_10, psiuD_dD11, FusedMulAddSIMD(FDPart3_17, psiuD_dD22, FusedMulAddSIMD(chip, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(MulSIMD(FDPart3_1, FDPart3_21), MulSIMD(rcoord, uu))), FusedMulSubSIMD(FDPart3_1, FDPart3_41, psiuD_dD00)))))))))))), MulSIMD(MulSIMD(FDPart3_21, FDPart3_Integer_2), MulSIMD(psiuD0, rcoord)));
           const REAL_SIMD_ARRAY __RHS_exp_23 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psiuD_dKOD11, FusedMulAddSIMD(FDPart3_Rational_1_2, psiuD_dD01, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psiuD_dKOD10, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psiuD_dKOD12, FusedMulSubSIMD(FDPart3_2, up_dD1, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psiuD1, chip)))))));
           const REAL_SIMD_ARRAY __RHS_exp_24 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psiuD_dKOD21, FusedMulAddSIMD(FDPart3_59, up_dD2, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psiuD_dKOD20, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psiuD_dKOD22, FusedMulSubSIMD(FDPart3_58, psiuD_dD02, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psiuD2, chip)))))));
           const REAL_SIMD_ARRAY __RHS_exp_25 = NegFusedMulAddSIMD(FDPart3_61, MulSIMD(MulSIMD(FDPart3_14, FDPart3_40), MulSIMD(FDPart3_Rational_1_2, Rprime)), SubSIMD(SubSIMD(FusedMulAddSIMD(FDPart3_63, MulSIMD(psigD0, gp), SubSIMD(FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_6, psifD0), FusedMulAddSIMD(FDPart3_30, MulSIMD(FDPart3_7, gg), FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psifD_dKOD02, FusedMulAddSIMD(FDPart3_1, MulSIMD(FDPart3_22, fgaugep), FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psifD_dKOD00, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psifD_dKOD01, FusedMulAddSIMD(FDPart3_31, fgauge, FusedMulAddSIMD(FDPart3_63, fgaugep, FusedMulAddSIMD(FDPart3_17, psifD_dD22, FusedMulAddSIMD(FDPart3_19, psifD1, FusedMulAddSIMD(FDPart3_1, FDPart3_49, FusedMulAddSIMD(FDPart3_10, psifD_dD11, FusedMulAddSIMD(chip, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(MulSIMD(FDPart3_1, FDPart3_21), MulSIMD(fgauge, rcoord))), SubSIMD(FDPart3_62, psifD_dD00)))))))))))))), MulSIMD(MulSIMD(FDPart3_1, FDPart3_27), MulSIMD(FDPart3_Rational_1_2, Rprime)))), MulSIMD(MulSIMD(FDPart3_37, FDPart3_61), MulSIMD(FDPart3_Rational_1_2, Rprime))), MulSIMD(MulSIMD(FDPart3_21, FDPart3_Integer_2), MulSIMD(psifD0, rcoord))));
           const REAL_SIMD_ARRAY __RHS_exp_26 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psifD_dKOD11, FusedMulAddSIMD(FDPart3_Rational_1_2, psifD_dD01, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psifD_dKOD10, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psifD_dKOD12, FusedMulSubSIMD(FDPart3_2, fgaugep_dD1, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psifD1, chip)))))));
           const REAL_SIMD_ARRAY __RHS_exp_27 = FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psifD_dKOD21, FusedMulAddSIMD(FDPart3_59, fgaugep_dD2, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psifD_dKOD20, FusedMulAddSIMD(FDPart3_Rational_5764607523034235_288230376151711744, psifD_dKOD22, FusedMulSubSIMD(FDPart3_58, psifD_dD02, MulSIMD(MulSIMD(FDPart3_1, FDPart3_Rational_1_2), MulSIMD(psifD2, chip)))))));
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
