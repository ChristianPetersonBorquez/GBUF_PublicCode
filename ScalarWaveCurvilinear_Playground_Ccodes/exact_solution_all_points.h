/*
 * Part P4: Declare the function for the exact solution at all points. time==0 corresponds to the initial data.
 */
void exact_solution_all_points(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict in_gfs) {
#include "./set_Cparameters.h"

  #pragma omp parallel for
  for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
    for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
      for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
        exact_solution_single_point(xx[0][i0],xx[1][i1],xx[2][i2],params,
                                    &in_gfs[IDX4S(GGGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(GPGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(BBGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(BPGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(UUGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(UPGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(FGAUGEGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(FGAUGEPGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(CHIGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(CHIPGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(CHIPPGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(RCOORDGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(RADIUSGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(RPRIMEGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(SINTHGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(COSTHGF,i0,i1,i2)],
                                    &in_gfs[IDX4S(PSIGD0GF,i0,i1,i2)],
                                    &in_gfs[IDX4S(PSIGD1GF,i0,i1,i2)],
                                    &in_gfs[IDX4S(PSIGD2GF,i0,i1,i2)],
                                    &in_gfs[IDX4S(PSIBD0GF,i0,i1,i2)],
                                    &in_gfs[IDX4S(PSIBD1GF,i0,i1,i2)],
                                    &in_gfs[IDX4S(PSIBD2GF,i0,i1,i2)],
                                    &in_gfs[IDX4S(PSIUD0GF,i0,i1,i2)],
                                    &in_gfs[IDX4S(PSIUD1GF,i0,i1,i2)],
                                    &in_gfs[IDX4S(PSIUD2GF,i0,i1,i2)],
                                    &in_gfs[IDX4S(PSIFD0GF,i0,i1,i2)],
                                    &in_gfs[IDX4S(PSIFD1GF,i0,i1,i2)],
                                    &in_gfs[IDX4S(PSIFD2GF,i0,i1,i2)]);
      } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
    } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
  } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
