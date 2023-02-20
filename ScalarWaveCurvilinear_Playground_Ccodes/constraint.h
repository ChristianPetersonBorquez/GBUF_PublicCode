/*
 * Evaluate the constraint
 */
void constraint(rfm_struct *restrict rfmstruct, const paramstruct *restrict params,
                  REAL *restrict in_gfs, REAL *restrict aux_gfs) {
#include "./set_Cparameters.h"

  #pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++) {
    #include "rfm_files/rfm_struct__read2.h"
    for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++) {
      #include "rfm_files/rfm_struct__read1.h"
      for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++) {
        #include "rfm_files/rfm_struct__read0.h"
        {
           /*
            * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
            */
           const double gg_i0_i1_i2m1 = in_gfs[IDX4S(GGGF, i0,i1,i2-1)];
           const double gg_i0_i1m1_i2 = in_gfs[IDX4S(GGGF, i0,i1-1,i2)];
           const double gg_i0m1_i1_i2 = in_gfs[IDX4S(GGGF, i0-1,i1,i2)];
           const double gg = in_gfs[IDX4S(GGGF, i0,i1,i2)];
           const double gg_i0p1_i1_i2 = in_gfs[IDX4S(GGGF, i0+1,i1,i2)];
           const double gg_i0_i1p1_i2 = in_gfs[IDX4S(GGGF, i0,i1+1,i2)];
           const double gg_i0_i1_i2p1 = in_gfs[IDX4S(GGGF, i0,i1,i2+1)];
           const double gp = in_gfs[IDX4S(GPGF, i0,i1,i2)];
           const double bb_i0_i1_i2m1 = in_gfs[IDX4S(BBGF, i0,i1,i2-1)];
           const double bb_i0_i1m1_i2 = in_gfs[IDX4S(BBGF, i0,i1-1,i2)];
           const double bb_i0m1_i1_i2 = in_gfs[IDX4S(BBGF, i0-1,i1,i2)];
           const double bb = in_gfs[IDX4S(BBGF, i0,i1,i2)];
           const double bb_i0p1_i1_i2 = in_gfs[IDX4S(BBGF, i0+1,i1,i2)];
           const double bb_i0_i1p1_i2 = in_gfs[IDX4S(BBGF, i0,i1+1,i2)];
           const double bb_i0_i1_i2p1 = in_gfs[IDX4S(BBGF, i0,i1,i2+1)];
           const double bp = in_gfs[IDX4S(BPGF, i0,i1,i2)];
           const double uu_i0_i1_i2m1 = in_gfs[IDX4S(UUGF, i0,i1,i2-1)];
           const double uu_i0_i1m1_i2 = in_gfs[IDX4S(UUGF, i0,i1-1,i2)];
           const double uu_i0m1_i1_i2 = in_gfs[IDX4S(UUGF, i0-1,i1,i2)];
           const double uu = in_gfs[IDX4S(UUGF, i0,i1,i2)];
           const double uu_i0p1_i1_i2 = in_gfs[IDX4S(UUGF, i0+1,i1,i2)];
           const double uu_i0_i1p1_i2 = in_gfs[IDX4S(UUGF, i0,i1+1,i2)];
           const double uu_i0_i1_i2p1 = in_gfs[IDX4S(UUGF, i0,i1,i2+1)];
           const double up = in_gfs[IDX4S(UPGF, i0,i1,i2)];
           const double fgauge_i0_i1_i2m1 = in_gfs[IDX4S(FGAUGEGF, i0,i1,i2-1)];
           const double fgauge_i0_i1m1_i2 = in_gfs[IDX4S(FGAUGEGF, i0,i1-1,i2)];
           const double fgauge_i0m1_i1_i2 = in_gfs[IDX4S(FGAUGEGF, i0-1,i1,i2)];
           const double fgauge = in_gfs[IDX4S(FGAUGEGF, i0,i1,i2)];
           const double fgauge_i0p1_i1_i2 = in_gfs[IDX4S(FGAUGEGF, i0+1,i1,i2)];
           const double fgauge_i0_i1p1_i2 = in_gfs[IDX4S(FGAUGEGF, i0,i1+1,i2)];
           const double fgauge_i0_i1_i2p1 = in_gfs[IDX4S(FGAUGEGF, i0,i1,i2+1)];
           const double fgaugep = in_gfs[IDX4S(FGAUGEPGF, i0,i1,i2)];
           const double chi = in_gfs[IDX4S(CHIGF, i0,i1,i2)];
           const double chip = in_gfs[IDX4S(CHIPGF, i0,i1,i2)];
           const double psigD0 = in_gfs[IDX4S(PSIGD0GF, i0,i1,i2)];
           const double psigD1 = in_gfs[IDX4S(PSIGD1GF, i0,i1,i2)];
           const double psigD2 = in_gfs[IDX4S(PSIGD2GF, i0,i1,i2)];
           const double psibD0 = in_gfs[IDX4S(PSIBD0GF, i0,i1,i2)];
           const double psibD1 = in_gfs[IDX4S(PSIBD1GF, i0,i1,i2)];
           const double psibD2 = in_gfs[IDX4S(PSIBD2GF, i0,i1,i2)];
           const double psiuD0 = in_gfs[IDX4S(PSIUD0GF, i0,i1,i2)];
           const double psiuD1 = in_gfs[IDX4S(PSIUD1GF, i0,i1,i2)];
           const double psiuD2 = in_gfs[IDX4S(PSIUD2GF, i0,i1,i2)];
           const double psifD0 = in_gfs[IDX4S(PSIFD0GF, i0,i1,i2)];
           const double psifD1 = in_gfs[IDX4S(PSIFD1GF, i0,i1,i2)];
           const double psifD2 = in_gfs[IDX4S(PSIFD2GF, i0,i1,i2)];
           const double Rprime = in_gfs[IDX4S(RPRIMEGF, i0,i1,i2)];
           const double sinth = in_gfs[IDX4S(SINTHGF, i0,i1,i2)];
           const double FDPart1_Rational_1_2 = 1.0/2.0;
           const double FDPart1_0 = FDPart1_Rational_1_2*invdx0;
           const double FDPart1_1 = FDPart1_Rational_1_2*invdx1;
           const double FDPart1_2 = FDPart1_Rational_1_2*invdx2;
           const double bb_dD0 = FDPart1_0*(-bb_i0m1_i1_i2 + bb_i0p1_i1_i2);
           const double bb_dD1 = FDPart1_1*(-bb_i0_i1m1_i2 + bb_i0_i1p1_i2);
           const double bb_dD2 = FDPart1_2*(-bb_i0_i1_i2m1 + bb_i0_i1_i2p1);
           const double fgauge_dD0 = FDPart1_0*(-fgauge_i0m1_i1_i2 + fgauge_i0p1_i1_i2);
           const double fgauge_dD1 = FDPart1_1*(-fgauge_i0_i1m1_i2 + fgauge_i0_i1p1_i2);
           const double fgauge_dD2 = FDPart1_2*(-fgauge_i0_i1_i2m1 + fgauge_i0_i1_i2p1);
           const double gg_dD0 = FDPart1_0*(-gg_i0m1_i1_i2 + gg_i0p1_i1_i2);
           const double gg_dD1 = FDPart1_1*(-gg_i0_i1m1_i2 + gg_i0_i1p1_i2);
           const double gg_dD2 = FDPart1_2*(-gg_i0_i1_i2m1 + gg_i0_i1_i2p1);
           const double uu_dD0 = FDPart1_0*(-uu_i0m1_i1_i2 + uu_i0p1_i1_i2);
           const double uu_dD1 = FDPart1_1*(-uu_i0_i1m1_i2 + uu_i0_i1p1_i2);
           const double uu_dD2 = FDPart1_2*(-uu_i0_i1_i2m1 + uu_i0_i1_i2p1);
           /*
            * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
            */
           const double FDPart3_0 = (1.0/(chi));
           const double FDPart3_1 = (1.0/(2*Rprime - 1));
           const double FDPart3_2 = FDPart3_0*chip;
           aux_gfs[IDX4S(CDG0GF, i0, i1, i2)] = -FDPart3_0*gp + FDPart3_1*(-FDPart3_2*gg + 2*gg_dD0 + psigD0);
           aux_gfs[IDX4S(CDG1GF, i0, i1, i2)] = -gg_dD1 + psigD1;
           aux_gfs[IDX4S(CDG2GF, i0, i1, i2)] = -gg_dD2 + psigD2*sinth;
           aux_gfs[IDX4S(CDB0GF, i0, i1, i2)] = -FDPart3_0*bp + FDPart3_1*(-FDPart3_2*bb + 2*bb_dD0 + psibD0);
           aux_gfs[IDX4S(CDB1GF, i0, i1, i2)] = -bb_dD1 + psibD1;
           aux_gfs[IDX4S(CDB2GF, i0, i1, i2)] = -bb_dD2 + psibD2*sinth;
           aux_gfs[IDX4S(CDU0GF, i0, i1, i2)] = -FDPart3_0*up + FDPart3_1*(-FDPart3_2*uu + psiuD0 + 2*uu_dD0);
           aux_gfs[IDX4S(CDU1GF, i0, i1, i2)] = psiuD1 - uu_dD1;
           aux_gfs[IDX4S(CDU2GF, i0, i1, i2)] = psiuD2*sinth - uu_dD2;
           aux_gfs[IDX4S(CDF0GF, i0, i1, i2)] = -FDPart3_0*fgaugep + FDPart3_1*(-FDPart3_2*fgauge + 2*fgauge_dD0 + psifD0);
           aux_gfs[IDX4S(CDF1GF, i0, i1, i2)] = -fgauge_dD1 + psifD1;
           aux_gfs[IDX4S(CDF2GF, i0, i1, i2)] = -fgauge_dD2 + psifD2*sinth;
        }
      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
