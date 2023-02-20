
// Declare boundary condition BC_UPDATE_OUTER macro,
//          which updates a single outer boundary face
//          of the 3D grid cube using quadratic polynomial
//          extrapolation.

#define BC_UPDATE_OUTER(which_gf, i0,i1,i2, FACEX0,FACEX1,FACEX2) {     \
    const int idx3 = IDX3S(i0,i1,i2);                                   \
    gfs[IDX4S(which_gf,i0,i1,i2)] =                                     \
        +3.0*gfs[IDX4S(which_gf,i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2)]   \
        -3.0*gfs[IDX4S(which_gf,i0+2*FACEX0,i1+2*FACEX1,i2+2*FACEX2)]   \
        +1.0*gfs[IDX4S(which_gf,i0+3*FACEX0,i1+3*FACEX1,i2+3*FACEX2)];  \
    }

// Curvilinear boundary condition driver routine: Apply BCs to all six
//          boundary faces of the 3D numerical domain, filling in the
//          innermost ghost zone layer first, and moving outward.

void apply_bcs_curvilinear_single_gf(const paramstruct *restrict params, const bc_struct *restrict bcstruct, const int8_t *restrict gfs_parity, const int which_gf, REAL *restrict gfs) {

#include "RELATIVE_PATH__set_Cparameters.h" /* Header file containing correct #include for set_Cparameters.h;
                                             * accounting for the relative path */

  for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
    // First apply OUTER boundary conditions,
    //   in case an INNER (parity) boundary point
    //   needs data at the outer boundary:
    // After updating each face, adjust imin[] and imax[]
    //   to reflect the newly-updated face extents.
    for(int pt=0;pt<bcstruct->num_ob_gz_pts[which_gz];pt++) {
      BC_UPDATE_OUTER(which_gf,
                      bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i0,
                      bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i1,
                      bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i2,
                      bcstruct->outer[which_gz][pt].FACEi0,
                      bcstruct->outer[which_gz][pt].FACEi1,
                      bcstruct->outer[which_gz][pt].FACEi2);
    }
    //printf("%d\n", which_gf);

    // Then apply INNER (parity) boundary conditions:
    for(int pt=0;pt<bcstruct->num_ib_gz_pts[which_gz];pt++) {
      const int i0dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i0;
      const int i1dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i1;
      const int i2dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i2;
      const int i0src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i0;
      const int i1src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i1;
      const int i2src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i2;
      const int8_t *prty= bcstruct->inner[which_gz][pt].parity;
      //printf("%d\n", i0dest);
      //                printf("%d\n",bcstruct->inner_bc_parity[which_gz][pt].parity[gfs_parity[which_gf]]);
      gfs[IDX4S(which_gf,i0dest,i1dest,i2dest)] =
        bcstruct->inner[which_gz][pt].parity[gfs_parity[which_gf]] * gfs[IDX4S(which_gf, i0src,i1src,i2src)];
    } // END for(int pt=0;pt<num_ib_gz_pts[which_gz];pt++)
  } // END for(int which_gz = 0; which_gz < NGHOSTS; which_gz++)
} // END function

void apply_bcs_curvilinear(const paramstruct *restrict params, const bc_struct *restrict bcstruct,
                           const int NUM_GFS, const int8_t *restrict gfs_parity, REAL *restrict gfs) {
#include "RELATIVE_PATH__set_Cparameters.h"
#pragma omp parallel for
    for(int which_gf=0;which_gf<NUM_GFS;which_gf++) {

        apply_bcs_curvilinear_single_gf(params, bcstruct, gfs_parity, which_gf, gfs);

    } // END for(int which_gf=0;which_gf<NUM_GFS;which_gf++)



// THIS IS A SECOND ATTEMPT TO GET THE CORRECT PARITY CONDITIONS AT THE ORIGIN. STARTING FROM SCRATCH. BASED ON CHRISTIAN'S NOTES
// ALL: 4 STANDS FOR CHI (chi) GRIDFUNCTION, 5 STANDS FOR CHIP (chip) GRIDFUNCTION,


    // GOOD FUNCTIONS: 10 STANDS FOR GZERO (gg) GRIDFUNCTION, 11 STANDS FOR GPLUS (gp) GRIDFUNCTION, 18 STANDS FOR GMINUS (psig[0]) GRIDFUNCTION

    // THIS FIRST PART IS THETA LEFT GHOSTPOINTS, PHI LEFT GHOST POINTS  // THIS IS WRONG IN PHI! //THINK THE ERROR IS CORRECTED NOW
    //FIRST -R POINTS
    gfs[IDX4S(11,1,1,1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(18,1,1,1)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(11,1,1,0)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(18,1,1,0)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(11,1,0,1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(18,1,0,1)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(11,1,0,0)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(18,1,0,0)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    //SECOND -R POINTS
    gfs[IDX4S(11,0,1,1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(18,0,1,1)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(11,0,1,0)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(18,0,1,0)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(11,0,0,1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(18,0,0,1)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(11,0,0,0)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(18,0,0,0)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA RIGHT GHOSTPOINTS, PHI LEFT GHOST POINTS
    //FIRST -R POINTS
    gfs[IDX4S(11,1,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(18,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(18,1,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(11,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(11,1,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(18,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(18,1,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(11,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(11,1,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(18,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(18,1,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(11,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(11,1,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(18,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(18,1,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(11,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    //SECOND -R POINTS
    gfs[IDX4S(11,0,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(18,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(18,0,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(11,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(11,0,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(18,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(18,0,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(11,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(11,0,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(18,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(18,0,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(11,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(10,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(11,0,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(18,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(18,0,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(11,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(10,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];



// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA LEFT GHOSTPOINTS, PHI RIGHT GHOST POINTS
    //FIRST -R POINTS
    gfs[IDX4S(11,1,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];
    gfs[IDX4S(18,1,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];

    gfs[IDX4S(11,1,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];
    gfs[IDX4S(18,1,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];

    gfs[IDX4S(11,1,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];
    gfs[IDX4S(18,1,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];

    gfs[IDX4S(11,1,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];
    gfs[IDX4S(18,1,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];

    //SECOND -R POINTS
    gfs[IDX4S(11,0,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];
    gfs[IDX4S(18,0,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];

    gfs[IDX4S(11,0,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];
    gfs[IDX4S(18,0,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];

    gfs[IDX4S(11,0,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];
    gfs[IDX4S(18,0,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];

    gfs[IDX4S(11,0,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];
    gfs[IDX4S(18,0,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA RIGHT GHOSTPOINTS, PHI RIGHT GHOST POINTS
    //FIRST -R POINTS
    gfs[IDX4S(11,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,2,2)]*gfs[IDX4S(18,2,2,2)]-gfs[IDX4S(5,2,2,2)]*gfs[IDX4S(10,2,2,2)];
    gfs[IDX4S(18,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(11,2,2,2)]/gfs[IDX4S(4,2,2,2)]-gfs[IDX4S(5,2,2,2)]*gfs[IDX4S(10,2,2,2)]/gfs[IDX4S(4,2,2,2)];

    gfs[IDX4S(11,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,2,3)]*gfs[IDX4S(18,2,2,3)]-gfs[IDX4S(5,2,2,3)]*gfs[IDX4S(10,2,2,3)];
    gfs[IDX4S(18,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(11,2,2,3)]/gfs[IDX4S(4,2,2,3)]-gfs[IDX4S(5,2,2,3)]*gfs[IDX4S(10,2,2,3)]/gfs[IDX4S(4,2,2,3)];

    gfs[IDX4S(11,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,3,2)]*gfs[IDX4S(18,2,3,2)]-gfs[IDX4S(5,2,3,2)]*gfs[IDX4S(10,2,3,2)];
    gfs[IDX4S(18,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(11,2,3,2)]/gfs[IDX4S(4,2,3,2)]-gfs[IDX4S(5,2,3,2)]*gfs[IDX4S(10,2,3,2)]/gfs[IDX4S(4,2,3,2)];

    gfs[IDX4S(11,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,3,3)]*gfs[IDX4S(18,2,3,3)]-gfs[IDX4S(5,2,3,3)]*gfs[IDX4S(10,2,3,3)];
    gfs[IDX4S(18,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(11,2,3,3)]/gfs[IDX4S(4,2,3,3)]-gfs[IDX4S(5,2,3,3)]*gfs[IDX4S(10,2,3,3)]/gfs[IDX4S(4,2,3,3)];

    //SECOND -R POINTS
    gfs[IDX4S(11,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,2,2)]*gfs[IDX4S(18,3,2,2)]-gfs[IDX4S(5,3,2,2)]*gfs[IDX4S(10,3,2,2)];
    gfs[IDX4S(18,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(11,3,2,2)]/gfs[IDX4S(4,3,2,2)]-gfs[IDX4S(5,3,2,2)]*gfs[IDX4S(10,3,2,2)]/gfs[IDX4S(4,3,2,2)];

    gfs[IDX4S(11,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,2,3)]*gfs[IDX4S(18,3,2,3)]-gfs[IDX4S(5,3,2,3)]*gfs[IDX4S(10,3,2,3)];
    gfs[IDX4S(18,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(11,3,2,3)]/gfs[IDX4S(4,3,2,3)]-gfs[IDX4S(5,3,2,3)]*gfs[IDX4S(10,3,2,3)]/gfs[IDX4S(4,3,2,3)];

    gfs[IDX4S(11,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,3,2)]*gfs[IDX4S(18,3,3,2)]-gfs[IDX4S(5,3,3,2)]*gfs[IDX4S(10,3,3,2)];
    gfs[IDX4S(18,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(11,3,3,2)]/gfs[IDX4S(4,3,3,2)]-gfs[IDX4S(5,3,3,2)]*gfs[IDX4S(10,3,3,2)]/gfs[IDX4S(4,3,3,2)];

    gfs[IDX4S(11,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,3,3)]*gfs[IDX4S(18,3,3,3)]-gfs[IDX4S(5,3,3,3)]*gfs[IDX4S(10,3,3,3)];
    gfs[IDX4S(18,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(11,3,3,3)]/gfs[IDX4S(4,3,3,3)]-gfs[IDX4S(5,3,3,3)]*gfs[IDX4S(10,3,3,3)]/gfs[IDX4S(4,3,3,3)];

// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA LEFT GHOSTPOINTS, PHI IN THE GRID
    for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++){
    //FIRST -R POINTS
     gfs[IDX4S(11,1,1,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];
     gfs[IDX4S(18,1,1,i2d)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];

     gfs[IDX4S(11,1,0,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];
     gfs[IDX4S(18,1,0,i2d)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];


    //SECOND -R POINTS
     gfs[IDX4S(11,0,1,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];
     gfs[IDX4S(18,0,1,i2d)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];

     gfs[IDX4S(11,0,0,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];
     gfs[IDX4S(18,0,0,i2d)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];

    } // END for i2d


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA RIGHT GHOSTPOINTS, PHI IN THE GRID
    for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++){
    //FIRST -R POINTS
     gfs[IDX4S(11,1,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(4,2,2,i2d)]*gfs[IDX4S(18,2,2,i2d)]-gfs[IDX4S(5,2,2,i2d)]*gfs[IDX4S(10,2,2,i2d)];
     gfs[IDX4S(18,1,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(11,2,2,i2d)]/gfs[IDX4S(4,2,2,i2d)]-gfs[IDX4S(5,2,2,i2d)]*gfs[IDX4S(10,2,2,i2d)]/gfs[IDX4S(4,2,2,i2d)];

     gfs[IDX4S(11,1,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(4,2,3,i2d)]*gfs[IDX4S(18,2,3,i2d)]-gfs[IDX4S(5,2,3,i2d)]*gfs[IDX4S(10,2,3,i2d)];
     gfs[IDX4S(18,1,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(11,2,3,i2d)]/gfs[IDX4S(4,2,3,i2d)]-gfs[IDX4S(5,2,3,i2d)]*gfs[IDX4S(10,2,3,i2d)]/gfs[IDX4S(4,2,3,i2d)];


    //SECOND -R POINTS
     gfs[IDX4S(11,0,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(4,3,2,i2d)]*gfs[IDX4S(18,3,2,i2d)]-gfs[IDX4S(5,3,2,i2d)]*gfs[IDX4S(10,3,2,i2d)];
     gfs[IDX4S(18,0,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(11,3,2,i2d)]/gfs[IDX4S(4,3,2,i2d)]-gfs[IDX4S(5,3,2,i2d)]*gfs[IDX4S(10,3,2,i2d)]/gfs[IDX4S(4,3,2,i2d)];

     gfs[IDX4S(11,0,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(4,3,3,i2d)]*gfs[IDX4S(18,3,3,i2d)]-gfs[IDX4S(5,3,3,i2d)]*gfs[IDX4S(10,3,3,i2d)];
     gfs[IDX4S(18,0,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(11,3,3,i2d)]/gfs[IDX4S(4,3,3,i2d)]-gfs[IDX4S(5,3,3,i2d)]*gfs[IDX4S(10,3,3,i2d)]/gfs[IDX4S(4,3,3,i2d)];

    }



// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA IN THE GRID, PHI LEFT GHOST POINTS
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++) {
    //FIRST -R POINTS
     gfs[IDX4S(11,1,i1d,1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];
     gfs[IDX4S(18,1,i1d,1)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];

     gfs[IDX4S(11,1,i1d,0)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];
     gfs[IDX4S(18,1,i1d,0)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];


    //SECOND -R POINTS
     gfs[IDX4S(11,0,i1d,1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];
     gfs[IDX4S(18,0,i1d,1)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];

     gfs[IDX4S(11,0,i1d,0)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];
     gfs[IDX4S(18,0,i1d,0)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];

    }


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA IN THE GRID, PHI RIGHT GHOST POINTS
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++){
    //FIRST -R POINTS
     gfs[IDX4S(11,1,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];
     gfs[IDX4S(18,1,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];

     gfs[IDX4S(11,1,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];
     gfs[IDX4S(18,1,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];


    //SECOND -R POINTS
     gfs[IDX4S(11,0,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];
     gfs[IDX4S(18,0,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];

     gfs[IDX4S(11,0,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];
     gfs[IDX4S(18,0,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];

    }


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA IN THE GRID, PHI IN THE GRID
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++){
        for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++){
            if( i2d<Nxx_plus_2NGHOSTS2/2 ){
     //FIRST -R POINTS
      gfs[IDX4S(11,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(18,1,i1d,i2d)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];

      gfs[IDX4S(11,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(18,1,i1d,i2d)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];

     //SECOND -R POINTS
      gfs[IDX4S(11,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(18,0,i1d,i2d)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];

      gfs[IDX4S(11,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(18,0,i1d,i2d)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
            }
            else{
         //FIRST -R POINTS
      gfs[IDX4S(11,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(18,1,i1d,i2d)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];

      gfs[IDX4S(11,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(18,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(18,1,i1d,i2d)] = gfs[IDX4S(11,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(10,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];

     //SECOND -R POINTS
      gfs[IDX4S(11,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(18,0,i1d,i2d)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];

      gfs[IDX4S(11,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(18,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(18,0,i1d,i2d)] = gfs[IDX4S(11,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(10,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
            }
        } // END for i2d
    } // END for i1d


//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------



    // BAD FUNCTIONS: 2 STANDS FOR BZERO (bb) GRIDFUNCTION, 3 STANDS FOR BPLUS (bp) GRIDFUNCTION, 12 STANDS FOR BMINUS (psib[0]) GRIDFUNCTION

     // THIS FIRST PART IS THETA LEFT GHOSTPOINTS, PHI LEFT GHOST POINTS  // THIS IS WRONG IN PHI! //THINK THE ERROR IS CORRECTED NOW
    //FIRST -R POINTS
    gfs[IDX4S(3,1,1,1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(12,1,1,1)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(3,1,1,0)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(12,1,1,0)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(3,1,0,1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(12,1,0,1)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(3,1,0,0)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(12,1,0,0)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    //SECOND -R POINTS
    gfs[IDX4S(3,0,1,1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(12,0,1,1)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(3,0,1,0)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(12,0,1,0)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(3,0,0,1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(12,0,0,1)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(3,0,0,0)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(12,0,0,0)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA RIGHT GHOSTPOINTS, PHI LEFT GHOST POINTS
    //FIRST -R POINTS
    gfs[IDX4S(3,1,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(12,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(12,1,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(3,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(3,1,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(12,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(12,1,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(3,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(3,1,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(12,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(12,1,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(3,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(3,1,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(12,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(12,1,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(3,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    //SECOND -R POINTS
    gfs[IDX4S(3,0,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(12,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(12,0,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(3,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(3,0,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(12,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(12,0,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(3,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(3,0,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(12,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(12,0,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(3,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(2,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(3,0,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(12,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(12,0,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(3,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(2,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];



// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA LEFT GHOSTPOINTS, PHI RIGHT GHOST POINTS
    //FIRST -R POINTS
    gfs[IDX4S(3,1,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];
    gfs[IDX4S(12,1,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];

    gfs[IDX4S(3,1,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];
    gfs[IDX4S(12,1,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];

    gfs[IDX4S(3,1,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];
    gfs[IDX4S(12,1,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];

    gfs[IDX4S(3,1,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];
    gfs[IDX4S(12,1,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];

    //SECOND -R POINTS
    gfs[IDX4S(3,0,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];
    gfs[IDX4S(12,0,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];

    gfs[IDX4S(3,0,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];
    gfs[IDX4S(12,0,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];

    gfs[IDX4S(3,0,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];
    gfs[IDX4S(12,0,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];

    gfs[IDX4S(3,0,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];
    gfs[IDX4S(12,0,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA RIGHT GHOSTPOINTS, PHI RIGHT GHOST POINTS
    //FIRST -R POINTS
    gfs[IDX4S(3,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,2,2)]*gfs[IDX4S(12,2,2,2)]-gfs[IDX4S(5,2,2,2)]*gfs[IDX4S(2,2,2,2)];
    gfs[IDX4S(12,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(3,2,2,2)]/gfs[IDX4S(4,2,2,2)]-gfs[IDX4S(5,2,2,2)]*gfs[IDX4S(2,2,2,2)]/gfs[IDX4S(4,2,2,2)];

    gfs[IDX4S(3,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,2,3)]*gfs[IDX4S(12,2,2,3)]-gfs[IDX4S(5,2,2,3)]*gfs[IDX4S(2,2,2,3)];
    gfs[IDX4S(12,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(3,2,2,3)]/gfs[IDX4S(4,2,2,3)]-gfs[IDX4S(5,2,2,3)]*gfs[IDX4S(2,2,2,3)]/gfs[IDX4S(4,2,2,3)];

    gfs[IDX4S(3,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,3,2)]*gfs[IDX4S(12,2,3,2)]-gfs[IDX4S(5,2,3,2)]*gfs[IDX4S(2,2,3,2)];
    gfs[IDX4S(12,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(3,2,3,2)]/gfs[IDX4S(4,2,3,2)]-gfs[IDX4S(5,2,3,2)]*gfs[IDX4S(2,2,3,2)]/gfs[IDX4S(4,2,3,2)];

    gfs[IDX4S(3,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,3,3)]*gfs[IDX4S(12,2,3,3)]-gfs[IDX4S(5,2,3,3)]*gfs[IDX4S(2,2,3,3)];
    gfs[IDX4S(12,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(3,2,3,3)]/gfs[IDX4S(4,2,3,3)]-gfs[IDX4S(5,2,3,3)]*gfs[IDX4S(2,2,3,3)]/gfs[IDX4S(4,2,3,3)];

    //SECOND -R POINTS
    gfs[IDX4S(3,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,2,2)]*gfs[IDX4S(12,3,2,2)]-gfs[IDX4S(5,3,2,2)]*gfs[IDX4S(2,3,2,2)];
    gfs[IDX4S(12,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(3,3,2,2)]/gfs[IDX4S(4,3,2,2)]-gfs[IDX4S(5,3,2,2)]*gfs[IDX4S(2,3,2,2)]/gfs[IDX4S(4,3,2,2)];

    gfs[IDX4S(3,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,2,3)]*gfs[IDX4S(12,3,2,3)]-gfs[IDX4S(5,3,2,3)]*gfs[IDX4S(2,3,2,3)];
    gfs[IDX4S(12,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(3,3,2,3)]/gfs[IDX4S(4,3,2,3)]-gfs[IDX4S(5,3,2,3)]*gfs[IDX4S(2,3,2,3)]/gfs[IDX4S(4,3,2,3)];

    gfs[IDX4S(3,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,3,2)]*gfs[IDX4S(12,3,3,2)]-gfs[IDX4S(5,3,3,2)]*gfs[IDX4S(2,3,3,2)];
    gfs[IDX4S(12,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(3,3,3,2)]/gfs[IDX4S(4,3,3,2)]-gfs[IDX4S(5,3,3,2)]*gfs[IDX4S(2,3,3,2)]/gfs[IDX4S(4,3,3,2)];

    gfs[IDX4S(3,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,3,3)]*gfs[IDX4S(12,3,3,3)]-gfs[IDX4S(5,3,3,3)]*gfs[IDX4S(2,3,3,3)];
    gfs[IDX4S(12,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(3,3,3,3)]/gfs[IDX4S(4,3,3,3)]-gfs[IDX4S(5,3,3,3)]*gfs[IDX4S(2,3,3,3)]/gfs[IDX4S(4,3,3,3)];

// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA LEFT GHOSTPOINTS, PHI IN THE GRID
    for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++){
    //FIRST -R POINTS
     gfs[IDX4S(3,1,1,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];
     gfs[IDX4S(12,1,1,i2d)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];

     gfs[IDX4S(3,1,0,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];
     gfs[IDX4S(12,1,0,i2d)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];


    //SECOND -R POINTS
     gfs[IDX4S(3,0,1,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];
     gfs[IDX4S(12,0,1,i2d)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];

     gfs[IDX4S(3,0,0,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];
     gfs[IDX4S(12,0,0,i2d)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];

    } // END for i2d


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA RIGHT GHOSTPOINTS, PHI IN THE GRID
    for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++){
    //FIRST -R POINTS
     gfs[IDX4S(3,1,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(4,2,2,i2d)]*gfs[IDX4S(12,2,2,i2d)]-gfs[IDX4S(5,2,2,i2d)]*gfs[IDX4S(2,2,2,i2d)];
     gfs[IDX4S(12,1,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(3,2,2,i2d)]/gfs[IDX4S(4,2,2,i2d)]-gfs[IDX4S(5,2,2,i2d)]*gfs[IDX4S(2,2,2,i2d)]/gfs[IDX4S(4,2,2,i2d)];

     gfs[IDX4S(3,1,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(4,2,3,i2d)]*gfs[IDX4S(12,2,3,i2d)]-gfs[IDX4S(5,2,3,i2d)]*gfs[IDX4S(2,2,3,i2d)];
     gfs[IDX4S(12,1,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(3,2,3,i2d)]/gfs[IDX4S(4,2,3,i2d)]-gfs[IDX4S(5,2,3,i2d)]*gfs[IDX4S(2,2,3,i2d)]/gfs[IDX4S(4,2,3,i2d)];


    //SECOND -R POINTS
     gfs[IDX4S(3,0,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(4,3,2,i2d)]*gfs[IDX4S(12,3,2,i2d)]-gfs[IDX4S(5,3,2,i2d)]*gfs[IDX4S(2,3,2,i2d)];
     gfs[IDX4S(12,0,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(3,3,2,i2d)]/gfs[IDX4S(4,3,2,i2d)]-gfs[IDX4S(5,3,2,i2d)]*gfs[IDX4S(2,3,2,i2d)]/gfs[IDX4S(4,3,2,i2d)];

     gfs[IDX4S(3,0,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(4,3,3,i2d)]*gfs[IDX4S(12,3,3,i2d)]-gfs[IDX4S(5,3,3,i2d)]*gfs[IDX4S(2,3,3,i2d)];
     gfs[IDX4S(12,0,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(3,3,3,i2d)]/gfs[IDX4S(4,3,3,i2d)]-gfs[IDX4S(5,3,3,i2d)]*gfs[IDX4S(2,3,3,i2d)]/gfs[IDX4S(4,3,3,i2d)];

    }



// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA IN THE GRID, PHI LEFT GHOST POINTS
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++) {
    //FIRST -R POINTS
     gfs[IDX4S(3,1,i1d,1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];
     gfs[IDX4S(12,1,i1d,1)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];

     gfs[IDX4S(3,1,i1d,0)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];
     gfs[IDX4S(12,1,i1d,0)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];


    //SECOND -R POINTS
     gfs[IDX4S(3,0,i1d,1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];
     gfs[IDX4S(12,0,i1d,1)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];

     gfs[IDX4S(3,0,i1d,0)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];
     gfs[IDX4S(12,0,i1d,0)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];

    }


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA IN THE GRID, PHI RIGHT GHOST POINTS
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++){
    //FIRST -R POINTS
     gfs[IDX4S(3,1,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];
     gfs[IDX4S(12,1,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];

     gfs[IDX4S(3,1,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];
     gfs[IDX4S(12,1,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];


    //SECOND -R POINTS
     gfs[IDX4S(3,0,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];
     gfs[IDX4S(12,0,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];

     gfs[IDX4S(3,0,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];
     gfs[IDX4S(12,0,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];

    }


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA IN THE GRID, PHI IN THE GRID
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++){
        for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++){
            if( i2d<Nxx_plus_2NGHOSTS2/2 ){
     //FIRST -R POINTS
      gfs[IDX4S(3,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(12,1,i1d,i2d)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];

      gfs[IDX4S(3,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(12,1,i1d,i2d)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];

     //SECOND -R POINTS
      gfs[IDX4S(3,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(12,0,i1d,i2d)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];

      gfs[IDX4S(3,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(12,0,i1d,i2d)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
            }
            else{
         //FIRST -R POINTS
      gfs[IDX4S(3,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(12,1,i1d,i2d)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];

      gfs[IDX4S(3,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(12,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(12,1,i1d,i2d)] = gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];

     //SECOND -R POINTS
      gfs[IDX4S(3,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(12,0,i1d,i2d)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];

      gfs[IDX4S(3,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(12,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(12,0,i1d,i2d)] = gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
            }
        } // END for i2d
    } // END for i1d


//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------

    // UGLY FUNCTIONS: 27 STANDS FOR UZERO (uu) GRIDFUNCTION, 26 STANDS FOR UPLUS (up) GRIDFUNCTION, 21 STANDS FOR UMINUS (psiu[0]) GRIDFUNCTION

     // THIS FIRST PART IS THETA LEFT GHOSTPOINTS, PHI LEFT GHOST POINTS  // THIS IS WRONG IN PHI! //THINK THE ERROR IS CORRECTED NOW
    //FIRST -R POINTS
    gfs[IDX4S(26,1,1,1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(21,1,1,1)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(26,1,1,0)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(21,1,1,0)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(26,1,0,1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(21,1,0,1)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(26,1,0,0)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(21,1,0,0)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    //SECOND -R POINTS
    gfs[IDX4S(26,0,1,1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(21,0,1,1)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(26,0,1,0)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(21,0,1,0)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(26,0,0,1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(21,0,0,1)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(26,0,0,0)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(21,0,0,0)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA RIGHT GHOSTPOINTS, PHI LEFT GHOST POINTS
    //FIRST -R POINTS
    gfs[IDX4S(26,1,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(21,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(21,1,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(26,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(26,1,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(21,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(21,1,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(26,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(26,1,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(21,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(21,1,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(26,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(26,1,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(21,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(21,1,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(26,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    //SECOND -R POINTS
    gfs[IDX4S(26,0,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(21,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(21,0,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(26,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(26,0,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(21,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(21,0,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(26,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(26,0,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(21,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(21,0,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(26,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(27,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(26,0,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(21,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(21,0,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(26,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(27,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];



// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA LEFT GHOSTPOINTS, PHI RIGHT GHOST POINTS
    //FIRST -R POINTS
    gfs[IDX4S(26,1,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];
    gfs[IDX4S(21,1,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];

    gfs[IDX4S(26,1,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];
    gfs[IDX4S(21,1,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];

    gfs[IDX4S(26,1,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];
    gfs[IDX4S(21,1,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];

    gfs[IDX4S(26,1,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];
    gfs[IDX4S(21,1,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];

    //SECOND -R POINTS
    gfs[IDX4S(26,0,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];
    gfs[IDX4S(21,0,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];

    gfs[IDX4S(26,0,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];
    gfs[IDX4S(21,0,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];

    gfs[IDX4S(26,0,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];
    gfs[IDX4S(21,0,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];

    gfs[IDX4S(26,0,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];
    gfs[IDX4S(21,0,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA RIGHT GHOSTPOINTS, PHI RIGHT GHOST POINTS
    //FIRST -R POINTS
    gfs[IDX4S(26,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,2,2)]*gfs[IDX4S(21,2,2,2)]-gfs[IDX4S(5,2,2,2)]*gfs[IDX4S(27,2,2,2)];
    gfs[IDX4S(21,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(26,2,2,2)]/gfs[IDX4S(4,2,2,2)]-gfs[IDX4S(5,2,2,2)]*gfs[IDX4S(27,2,2,2)]/gfs[IDX4S(4,2,2,2)];

    gfs[IDX4S(26,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,2,3)]*gfs[IDX4S(21,2,2,3)]-gfs[IDX4S(5,2,2,3)]*gfs[IDX4S(27,2,2,3)];
    gfs[IDX4S(21,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(26,2,2,3)]/gfs[IDX4S(4,2,2,3)]-gfs[IDX4S(5,2,2,3)]*gfs[IDX4S(27,2,2,3)]/gfs[IDX4S(4,2,2,3)];

    gfs[IDX4S(26,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,3,2)]*gfs[IDX4S(21,2,3,2)]-gfs[IDX4S(5,2,3,2)]*gfs[IDX4S(27,2,3,2)];
    gfs[IDX4S(21,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(26,2,3,2)]/gfs[IDX4S(4,2,3,2)]-gfs[IDX4S(5,2,3,2)]*gfs[IDX4S(27,2,3,2)]/gfs[IDX4S(4,2,3,2)];

    gfs[IDX4S(26,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,3,3)]*gfs[IDX4S(21,2,3,3)]-gfs[IDX4S(5,2,3,3)]*gfs[IDX4S(27,2,3,3)];
    gfs[IDX4S(21,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(26,2,3,3)]/gfs[IDX4S(4,2,3,3)]-gfs[IDX4S(5,2,3,3)]*gfs[IDX4S(27,2,3,3)]/gfs[IDX4S(4,2,3,3)];

    //SECOND -R POINTS
    gfs[IDX4S(26,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,2,2)]*gfs[IDX4S(21,3,2,2)]-gfs[IDX4S(5,3,2,2)]*gfs[IDX4S(27,3,2,2)];
    gfs[IDX4S(21,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(26,3,2,2)]/gfs[IDX4S(4,3,2,2)]-gfs[IDX4S(5,3,2,2)]*gfs[IDX4S(27,3,2,2)]/gfs[IDX4S(4,3,2,2)];

    gfs[IDX4S(26,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,2,3)]*gfs[IDX4S(21,3,2,3)]-gfs[IDX4S(5,3,2,3)]*gfs[IDX4S(27,3,2,3)];
    gfs[IDX4S(21,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(26,3,2,3)]/gfs[IDX4S(4,3,2,3)]-gfs[IDX4S(5,3,2,3)]*gfs[IDX4S(27,3,2,3)]/gfs[IDX4S(4,3,2,3)];

    gfs[IDX4S(26,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,3,2)]*gfs[IDX4S(21,3,3,2)]-gfs[IDX4S(5,3,3,2)]*gfs[IDX4S(27,3,3,2)];
    gfs[IDX4S(21,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(26,3,3,2)]/gfs[IDX4S(4,3,3,2)]-gfs[IDX4S(5,3,3,2)]*gfs[IDX4S(27,3,3,2)]/gfs[IDX4S(4,3,3,2)];

    gfs[IDX4S(26,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,3,3)]*gfs[IDX4S(21,3,3,3)]-gfs[IDX4S(5,3,3,3)]*gfs[IDX4S(27,3,3,3)];
    gfs[IDX4S(21,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(26,3,3,3)]/gfs[IDX4S(4,3,3,3)]-gfs[IDX4S(5,3,3,3)]*gfs[IDX4S(27,3,3,3)]/gfs[IDX4S(4,3,3,3)];

// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA LEFT GHOSTPOINTS, PHI IN THE GRID
    for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++){
    //FIRST -R POINTS
     gfs[IDX4S(26,1,1,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];
     gfs[IDX4S(21,1,1,i2d)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];

     gfs[IDX4S(26,1,0,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];
     gfs[IDX4S(21,1,0,i2d)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];


    //SECOND -R POINTS
     gfs[IDX4S(26,0,1,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];
     gfs[IDX4S(21,0,1,i2d)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];

     gfs[IDX4S(26,0,0,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];
     gfs[IDX4S(21,0,0,i2d)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];

    } // END for i2d


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA RIGHT GHOSTPOINTS, PHI IN THE GRID
    for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++){
    //FIRST -R POINTS
     gfs[IDX4S(26,1,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(4,2,2,i2d)]*gfs[IDX4S(21,2,2,i2d)]-gfs[IDX4S(5,2,2,i2d)]*gfs[IDX4S(27,2,2,i2d)];
     gfs[IDX4S(21,1,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(26,2,2,i2d)]/gfs[IDX4S(4,2,2,i2d)]-gfs[IDX4S(5,2,2,i2d)]*gfs[IDX4S(27,2,2,i2d)]/gfs[IDX4S(4,2,2,i2d)];

     gfs[IDX4S(26,1,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(4,2,3,i2d)]*gfs[IDX4S(21,2,3,i2d)]-gfs[IDX4S(5,2,3,i2d)]*gfs[IDX4S(27,2,3,i2d)];
     gfs[IDX4S(21,1,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(26,2,3,i2d)]/gfs[IDX4S(4,2,3,i2d)]-gfs[IDX4S(5,2,3,i2d)]*gfs[IDX4S(27,2,3,i2d)]/gfs[IDX4S(4,2,3,i2d)];


    //SECOND -R POINTS
     gfs[IDX4S(26,0,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(4,3,2,i2d)]*gfs[IDX4S(21,3,2,i2d)]-gfs[IDX4S(5,3,2,i2d)]*gfs[IDX4S(27,3,2,i2d)];
     gfs[IDX4S(21,0,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(26,3,2,i2d)]/gfs[IDX4S(4,3,2,i2d)]-gfs[IDX4S(5,3,2,i2d)]*gfs[IDX4S(27,3,2,i2d)]/gfs[IDX4S(4,3,2,i2d)];

     gfs[IDX4S(26,0,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(4,3,3,i2d)]*gfs[IDX4S(21,3,3,i2d)]-gfs[IDX4S(5,3,3,i2d)]*gfs[IDX4S(27,3,3,i2d)];
     gfs[IDX4S(21,0,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(26,3,3,i2d)]/gfs[IDX4S(4,3,3,i2d)]-gfs[IDX4S(5,3,3,i2d)]*gfs[IDX4S(27,3,3,i2d)]/gfs[IDX4S(4,3,3,i2d)];

    }



// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA IN THE GRID, PHI LEFT GHOST POINTS
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++) {
    //FIRST -R POINTS
     gfs[IDX4S(26,1,i1d,1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];
     gfs[IDX4S(21,1,i1d,1)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];

     gfs[IDX4S(26,1,i1d,0)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];
     gfs[IDX4S(21,1,i1d,0)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];


    //SECOND -R POINTS
     gfs[IDX4S(26,0,i1d,1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];
     gfs[IDX4S(21,0,i1d,1)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];

     gfs[IDX4S(26,0,i1d,0)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];
     gfs[IDX4S(21,0,i1d,0)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];

    }


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA IN THE GRID, PHI RIGHT GHOST POINTS
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++){
    //FIRST -R POINTS
     gfs[IDX4S(26,1,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];
     gfs[IDX4S(21,1,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];

     gfs[IDX4S(26,1,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];
     gfs[IDX4S(21,1,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];


    //SECOND -R POINTS
     gfs[IDX4S(26,0,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];
     gfs[IDX4S(21,0,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];

     gfs[IDX4S(26,0,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];
     gfs[IDX4S(21,0,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];

    }


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA IN THE GRID, PHI IN THE GRID
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++){
        for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++){
            if( i2d<Nxx_plus_2NGHOSTS2/2 ){
     //FIRST -R POINTS
      gfs[IDX4S(26,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(21,1,i1d,i2d)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];

      gfs[IDX4S(26,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(21,1,i1d,i2d)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];

     //SECOND -R POINTS
      gfs[IDX4S(26,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(21,0,i1d,i2d)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];

      gfs[IDX4S(26,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(21,0,i1d,i2d)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
            }
            else{
         //FIRST -R POINTS
      gfs[IDX4S(26,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(21,1,i1d,i2d)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];

      gfs[IDX4S(26,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(21,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(21,1,i1d,i2d)] = gfs[IDX4S(26,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(27,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];

     //SECOND -R POINTS
      gfs[IDX4S(26,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(21,0,i1d,i2d)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];

      gfs[IDX4S(26,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(21,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(21,0,i1d,i2d)] = gfs[IDX4S(26,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(27,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
            }
        } // END for i2d
    } // END for i1d


//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------


    // F FUNCTIONS: 8 STANDS FOR FZERO (fgauge) GRIDFUNCTION, 9 STANDS FOR FPLUS (fgaugep) GRIDFUNCTION, 15 STANDS FOR FMINUS (psif[0]) GRIDFUNCTION

     // THIS FIRST PART IS THETA LEFT GHOSTPOINTS, PHI LEFT GHOST POINTS  // THIS IS WRONG IN PHI! //THINK THE ERROR IS CORRECTED NOW
    //FIRST -R POINTS
    gfs[IDX4S(9,1,1,1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(15,1,1,1)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(9,1,1,0)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(15,1,1,0)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(9,1,0,1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(15,1,0,1)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(9,1,0,0)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(15,1,0,0)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    //SECOND -R POINTS
    gfs[IDX4S(9,0,1,1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(15,0,1,1)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(9,0,1,0)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(15,0,1,0)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(9,0,0,1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(15,0,0,1)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(9,0,0,0)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(15,0,0,0)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA RIGHT GHOSTPOINTS, PHI LEFT GHOST POINTS
    //FIRST -R POINTS
    gfs[IDX4S(9,1,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(15,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(15,1,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(9,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(9,1,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(15,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(15,1,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(9,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(9,1,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(15,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(15,1,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(9,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(9,1,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(15,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(15,1,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(9,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,2,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    //SECOND -R POINTS
    gfs[IDX4S(9,0,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(15,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(15,0,Nxx_plus_2NGHOSTS1-2,1)] = gfs[IDX4S(9,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(9,0,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(15,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(15,0,Nxx_plus_2NGHOSTS1-2,0)] = gfs[IDX4S(9,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,2,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];

    gfs[IDX4S(9,0,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(15,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];
    gfs[IDX4S(15,0,Nxx_plus_2NGHOSTS1-1,1)] = gfs[IDX4S(9,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]*gfs[IDX4S(8,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-1)];

    gfs[IDX4S(9,0,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(15,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];
    gfs[IDX4S(15,0,Nxx_plus_2NGHOSTS1-1,0)] = gfs[IDX4S(9,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]-gfs[IDX4S(5,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]*gfs[IDX4S(8,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)]/gfs[IDX4S(4,3,3,Nxx_plus_2NGHOSTS2-NGHOSTS-2)];



// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA LEFT GHOSTPOINTS, PHI RIGHT GHOST POINTS
    //FIRST -R POINTS
    gfs[IDX4S(9,1,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];
    gfs[IDX4S(15,1,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];

    gfs[IDX4S(9,1,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];
    gfs[IDX4S(15,1,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];

    gfs[IDX4S(9,1,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];
    gfs[IDX4S(15,1,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];

    gfs[IDX4S(9,1,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];
    gfs[IDX4S(15,1,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];

    //SECOND -R POINTS
    gfs[IDX4S(9,0,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];
    gfs[IDX4S(15,0,1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,2)];

    gfs[IDX4S(9,0,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];
    gfs[IDX4S(15,0,1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,3)];

    gfs[IDX4S(9,0,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];
    gfs[IDX4S(15,0,0,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,2)];

    gfs[IDX4S(9,0,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];
    gfs[IDX4S(15,0,0,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,3)];


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA RIGHT GHOSTPOINTS, PHI RIGHT GHOST POINTS
    //FIRST -R POINTS
    gfs[IDX4S(9,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,2,2)]*gfs[IDX4S(15,2,2,2)]-gfs[IDX4S(5,2,2,2)]*gfs[IDX4S(8,2,2,2)];
    gfs[IDX4S(15,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(9,2,2,2)]/gfs[IDX4S(4,2,2,2)]-gfs[IDX4S(5,2,2,2)]*gfs[IDX4S(8,2,2,2)]/gfs[IDX4S(4,2,2,2)];

    gfs[IDX4S(9,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,2,3)]*gfs[IDX4S(15,2,2,3)]-gfs[IDX4S(5,2,2,3)]*gfs[IDX4S(8,2,2,3)];
    gfs[IDX4S(15,1,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(9,2,2,3)]/gfs[IDX4S(4,2,2,3)]-gfs[IDX4S(5,2,2,3)]*gfs[IDX4S(8,2,2,3)]/gfs[IDX4S(4,2,2,3)];

    gfs[IDX4S(9,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,3,2)]*gfs[IDX4S(15,2,3,2)]-gfs[IDX4S(5,2,3,2)]*gfs[IDX4S(8,2,3,2)];
    gfs[IDX4S(15,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(9,2,3,2)]/gfs[IDX4S(4,2,3,2)]-gfs[IDX4S(5,2,3,2)]*gfs[IDX4S(8,2,3,2)]/gfs[IDX4S(4,2,3,2)];

    gfs[IDX4S(9,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,3,3)]*gfs[IDX4S(15,2,3,3)]-gfs[IDX4S(5,2,3,3)]*gfs[IDX4S(8,2,3,3)];
    gfs[IDX4S(15,1,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(9,2,3,3)]/gfs[IDX4S(4,2,3,3)]-gfs[IDX4S(5,2,3,3)]*gfs[IDX4S(8,2,3,3)]/gfs[IDX4S(4,2,3,3)];

    //SECOND -R POINTS
    gfs[IDX4S(9,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,2,2)]*gfs[IDX4S(15,3,2,2)]-gfs[IDX4S(5,3,2,2)]*gfs[IDX4S(8,3,2,2)];
    gfs[IDX4S(15,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(9,3,2,2)]/gfs[IDX4S(4,3,2,2)]-gfs[IDX4S(5,3,2,2)]*gfs[IDX4S(8,3,2,2)]/gfs[IDX4S(4,3,2,2)];

    gfs[IDX4S(9,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,2,3)]*gfs[IDX4S(15,3,2,3)]-gfs[IDX4S(5,3,2,3)]*gfs[IDX4S(8,3,2,3)];
    gfs[IDX4S(15,0,Nxx_plus_2NGHOSTS1-2,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(9,3,2,3)]/gfs[IDX4S(4,3,2,3)]-gfs[IDX4S(5,3,2,3)]*gfs[IDX4S(8,3,2,3)]/gfs[IDX4S(4,3,2,3)];

    gfs[IDX4S(9,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,3,2)]*gfs[IDX4S(15,3,3,2)]-gfs[IDX4S(5,3,3,2)]*gfs[IDX4S(8,3,3,2)];
    gfs[IDX4S(15,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(9,3,3,2)]/gfs[IDX4S(4,3,3,2)]-gfs[IDX4S(5,3,3,2)]*gfs[IDX4S(8,3,3,2)]/gfs[IDX4S(4,3,3,2)];

    gfs[IDX4S(9,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,3,3)]*gfs[IDX4S(15,3,3,3)]-gfs[IDX4S(5,3,3,3)]*gfs[IDX4S(8,3,3,3)];
    gfs[IDX4S(15,0,Nxx_plus_2NGHOSTS1-1,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(9,3,3,3)]/gfs[IDX4S(4,3,3,3)]-gfs[IDX4S(5,3,3,3)]*gfs[IDX4S(8,3,3,3)]/gfs[IDX4S(4,3,3,3)];

// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA LEFT GHOSTPOINTS, PHI IN THE GRID
    for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++){
    //FIRST -R POINTS
     gfs[IDX4S(9,1,1,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];
     gfs[IDX4S(15,1,1,i2d)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];

     gfs[IDX4S(9,1,0,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];
     gfs[IDX4S(15,1,0,i2d)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];


    //SECOND -R POINTS
     gfs[IDX4S(9,0,1,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];
     gfs[IDX4S(15,0,1,i2d)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];

     gfs[IDX4S(9,0,0,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];
     gfs[IDX4S(15,0,0,i2d)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)];

    } // END for i2d


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA RIGHT GHOSTPOINTS, PHI IN THE GRID
    for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++){
    //FIRST -R POINTS
     gfs[IDX4S(9,1,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(4,2,2,i2d)]*gfs[IDX4S(15,2,2,i2d)]-gfs[IDX4S(5,2,2,i2d)]*gfs[IDX4S(8,2,2,i2d)];
     gfs[IDX4S(15,1,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(9,2,2,i2d)]/gfs[IDX4S(4,2,2,i2d)]-gfs[IDX4S(5,2,2,i2d)]*gfs[IDX4S(8,2,2,i2d)]/gfs[IDX4S(4,2,2,i2d)];

     gfs[IDX4S(9,1,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(4,2,3,i2d)]*gfs[IDX4S(15,2,3,i2d)]-gfs[IDX4S(5,2,3,i2d)]*gfs[IDX4S(8,2,3,i2d)];
     gfs[IDX4S(15,1,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(9,2,3,i2d)]/gfs[IDX4S(4,2,3,i2d)]-gfs[IDX4S(5,2,3,i2d)]*gfs[IDX4S(8,2,3,i2d)]/gfs[IDX4S(4,2,3,i2d)];


    //SECOND -R POINTS
     gfs[IDX4S(9,0,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(4,3,2,i2d)]*gfs[IDX4S(15,3,2,i2d)]-gfs[IDX4S(5,3,2,i2d)]*gfs[IDX4S(8,3,2,i2d)];
     gfs[IDX4S(15,0,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(9,3,2,i2d)]/gfs[IDX4S(4,3,2,i2d)]-gfs[IDX4S(5,3,2,i2d)]*gfs[IDX4S(8,3,2,i2d)]/gfs[IDX4S(4,3,2,i2d)];

     gfs[IDX4S(9,0,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(4,3,3,i2d)]*gfs[IDX4S(15,3,3,i2d)]-gfs[IDX4S(5,3,3,i2d)]*gfs[IDX4S(8,3,3,i2d)];
     gfs[IDX4S(15,0,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(9,3,3,i2d)]/gfs[IDX4S(4,3,3,i2d)]-gfs[IDX4S(5,3,3,i2d)]*gfs[IDX4S(8,3,3,i2d)]/gfs[IDX4S(4,3,3,i2d)];

    }



// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA IN THE GRID, PHI LEFT GHOST POINTS
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++) {
    //FIRST -R POINTS
     gfs[IDX4S(9,1,i1d,1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];
     gfs[IDX4S(15,1,i1d,1)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];

     gfs[IDX4S(9,1,i1d,0)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];
     gfs[IDX4S(15,1,i1d,0)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];


    //SECOND -R POINTS
     gfs[IDX4S(9,0,i1d,1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];
     gfs[IDX4S(15,0,i1d,1)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-1)];

     gfs[IDX4S(9,0,i1d,0)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];
     gfs[IDX4S(15,0,i1d,0)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2-2)];

    }


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA IN THE GRID, PHI RIGHT GHOST POINTS
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++){
    //FIRST -R POINTS
     gfs[IDX4S(9,1,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];
     gfs[IDX4S(15,1,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];

     gfs[IDX4S(9,1,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];
     gfs[IDX4S(15,1,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];


    //SECOND -R POINTS
     gfs[IDX4S(9,0,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];
     gfs[IDX4S(15,0,i1d,Nxx_plus_2NGHOSTS2-2)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2)];

     gfs[IDX4S(9,0,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];
     gfs[IDX4S(15,0,i1d,Nxx_plus_2NGHOSTS2-1)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,Nxx_plus_2NGHOSTS2/2+1)];

    }


// ---------------------------------------------------------------------------------------------------------------------

    // THIS PART IS THETA IN THE GRID, PHI IN THE GRID
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++){
        for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++){
            if( i2d<Nxx_plus_2NGHOSTS2/2 ){
     //FIRST -R POINTS
      gfs[IDX4S(9,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(15,1,i1d,i2d)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];

      gfs[IDX4S(9,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(15,1,i1d,i2d)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];

     //SECOND -R POINTS
      gfs[IDX4S(9,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(15,0,i1d,i2d)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];

      gfs[IDX4S(9,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
      gfs[IDX4S(15,0,i1d,i2d)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d+Nxx_plus_2NGHOSTS2/2-2)];
            }
            else{
         //FIRST -R POINTS
      gfs[IDX4S(9,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(15,1,i1d,i2d)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];

      gfs[IDX4S(9,1,i1d,i2d)] = gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(15,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(15,1,i1d,i2d)] = gfs[IDX4S(9,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,2,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];

     //SECOND -R POINTS
      gfs[IDX4S(9,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(15,0,i1d,i2d)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];

      gfs[IDX4S(9,0,i1d,i2d)] = gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(15,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
      gfs[IDX4S(15,0,i1d,i2d)] = gfs[IDX4S(9,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]-gfs[IDX4S(5,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)]/gfs[IDX4S(4,3,Nxx_plus_2NGHOSTS1-i1d-1,i2d-Nxx_plus_2NGHOSTS2/2+2)];
            }
        } // END for i2d
    } // END for i1d


//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------














//--------------------------------------PREVIOUS-ATTEMPTS----------------------------------------------------------

    //double temp;

    //IN THE FOLLOWING INDEX 7 STANDS FOR GPLUS (GP) GRIDFUNCTION, 8 FOR GMINUS (PSI[0])
    //APPARENTLY THE FIRST TWO METHODS WORK, AT LEAST IN SPHERICAL SYMMETRY
    //THIS IS HARDCODED FOR FD_ORDER=2 -> NGHOSTS=2

    // FIRST LAYER OF GHOST POINTS

    /*for(int pt=0;pt<bcstruct->num_ib_gz_pts[0];pt++) {

    const int i0d = bcstruct->inner[0][pt].inner_bc_dest_pt.i0;
    const int i1d = bcstruct->inner[0][pt].inner_bc_dest_pt.i1;
    const int i2d = bcstruct->inner[0][pt].inner_bc_dest_pt.i2;
    const int i0s  = bcstruct->inner[0][pt].inner_bc_src_pt.i0;
    const int i1s  = bcstruct->inner[0][pt].inner_bc_src_pt.i1;
    const int i2s  = bcstruct->inner[0][pt].inner_bc_src_pt.i2;

    // THE FOLLOWING ARE THE BCs FOR THE CHARACTERISTIC VARIABLES
    gfs[IDX4S(7,i0d,i1d,i2d)] = gfs[IDX4S(2,i0s,i1s,i2s)]*gfs[IDX4S(8,i0s,i1s,i2s)]-gfs[IDX4S(3,i0s,i1s,i2s)]*gfs[IDX4S(6,i0s,i1s,i2s)];
    gfs[IDX4S(8,i0d,i1d,i2d)] = gfs[IDX4S(7,i0s,i1s,i2s)]/gfs[IDX4S(2,i0s,i1s,i2s)]-gfs[IDX4S(3,i0s,i1s,i2s)]*gfs[IDX4S(6,i0s,i1s,i2s)]/gfs[IDX4S(2,i0s,i1s,i2s)];

    // THE FOLLOWING ARE THE BCs FOR THE GFs INTRODUCED FOR THE EVANS TRICK
    gfs[IDX4S(3,i0d,i1d,i2d)] = -gfs[IDX4S(3,i0s,i1s,i2s)]; // CHIP.. SHOULD THIS BE BEFORE FILLING THE OTHER GFs GHOSTPOINTS?
    gfs[IDX4S(0,i0d,i1d,i2d)] = -gfs[IDX4S(0,i0s,i1s,i2s)]; // RADIUS
    gfs[IDX4S(11,i0d,i1d,i2d)] = -gfs[IDX4S(11,i0s,i1s,i2s)]; // SINTH .. CHECK

    } */ // END for(int pt=0;pt<bcstruct->num_ib_gz_pts[which_gz];pt++)

    // THE FOLLOWING FILLS R NEGATIVE BUT THETA,PHI IN THE GRID (NOT THETA OR PHI GHOSTS)
    /*for(int i1d=0; i1d<Nxx_plus_2NGHOSTS1; i1d++){
        for(int i2d=0; i2d<Nxx_plus_2NGHOSTS2; i2d++){
            //THE FOLLOWING BCs SEEM TO WORK IN AXISYMMETRY BUT NOT IN FULL 3D
            gfs[IDX4S(11,1,i1d,i2d)] = gfs[IDX4S(4,2,i1d,i2d)]*gfs[IDX4S(18,2,i1d,i2d)]-gfs[IDX4S(5,2,i1d,i2d)]*gfs[IDX4S(10,2,i1d,i2d)];  //G+(-r)=chi(r)G-(r)-chi'(r)G(r)
            gfs[IDX4S(18,1,i1d,i2d)] = gfs[IDX4S(11,2,i1d,i2d)]/gfs[IDX4S(4,2,i1d,i2d)]-gfs[IDX4S(5,2,i1d,i2d)]*gfs[IDX4S(10,2,i1d,i2d)]/gfs[IDX4S(4,2,i1d,i2d)]; //G-(-r)=G+(r)/chi(r)-chi'(r)G(r)/chi(r)

            gfs[IDX4S(3,1,i1d,i2d)] = gfs[IDX4S(4,2,i1d,i2d)]*gfs[IDX4S(12,2,i1d,i2d)]-gfs[IDX4S(5,2,i1d,i2d)]*gfs[IDX4S(2,2,i1d,i2d)];
            gfs[IDX4S(12,1,i1d,i2d)] = gfs[IDX4S(3,2,i1d,i2d)]/gfs[IDX4S(4,2,i1d,i2d)]-gfs[IDX4S(5,2,i1d,i2d)]*gfs[IDX4S(2,2,i1d,i2d)]/gfs[IDX4S(4,2,i1d,i2d)];

            gfs[IDX4S(26,1,i1d,i2d)] = gfs[IDX4S(4,2,i1d,i2d)]*gfs[IDX4S(21,2,i1d,i2d)]-gfs[IDX4S(5,2,i1d,i2d)]*gfs[IDX4S(27,2,i1d,i2d)];
            gfs[IDX4S(21,1,i1d,i2d)] = gfs[IDX4S(26,2,i1d,i2d)]/gfs[IDX4S(4,2,i1d,i2d)]-gfs[IDX4S(5,2,i1d,i2d)]*gfs[IDX4S(27,2,i1d,i2d)]/gfs[IDX4S(4,2,i1d,i2d)];

            gfs[IDX4S(9,1,i1d,i2d)] = gfs[IDX4S(4,2,i1d,i2d)]*gfs[IDX4S(15,2,i1d,i2d)]-gfs[IDX4S(5,2,i1d,i2d)]*gfs[IDX4S(8,2,i1d,i2d)];
            gfs[IDX4S(15,1,i1d,i2d)] = gfs[IDX4S(9,2,i1d,i2d)]/gfs[IDX4S(4,2,i1d,i2d)]-gfs[IDX4S(5,2,i1d,i2d)]*gfs[IDX4S(8,2,i1d,i2d)]/gfs[IDX4S(4,2,i1d,i2d)];


            /*temp = gfs[IDX4S(2,2,i1d,i2d)]*gfs[IDX4S(8,1,i1d,i2d)]-gfs[IDX4S(3,2,i1d,i2d)]*gfs[IDX4S(6,1,i1d,i2d)];
            gfs[IDX4S(8,1,i1d,i2d)] = gfs[IDX4S(7,1,i1d,i2d)]/gfs[IDX4S(2,2,i1d,i2d)]-gfs[IDX4S(3,2,i1d,i2d)]*gfs[IDX4S(6,1,i1d,i2d)]/gfs[IDX4S(2,2,i1d,i2d)];
            gfs[IDX4S(7,1,i1d,i2d)] = temp;*/


            /*if(i2d<=(Nxx_plus_2NGHOSTS2/2)){
                gfs[IDX4S(7,1,i1d,i2d)] = gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]; // theta --> pi-theta, phi --> phi+pi
                gfs[IDX4S(8,1,i1d,i2d)] = gfs[IDX4S(7,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))];
            }
            if(i2d>(Nxx_plus_2NGHOSTS2/2)){
                gfs[IDX4S(7,1,i1d,i2d)] = gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]; // theta --> pi-theta, phi --> phi-pi
                gfs[IDX4S(8,1,i1d,i2d)] = gfs[IDX4S(7,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))];
            }*/
      //  } // END for i1d
   // }  // END for i2d

    /*
    //NEGATIVE THETA, PHI IN THE GRID ... not sure if its OK (check phi)
    for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++) {
        gfs[IDX4S(7,1,1,i2d)] = gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]; //-THETA -> PI-|THETA|.. Nxx_plus_2NGHOSTS1-NGHOSTS-1 OR Nxx_plus_2NGHOSTS1-NGHOSTS-2???
        gfs[IDX4S(8,1,1,i2d)] = gfs[IDX4S(7,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];
    } //end for phi loop (i2d)

    // THETA > PI, PHI IN THE GRID ... not sure if its OK (check phi)
    for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++) {
        gfs[IDX4S(7,1,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(2,2,2,i2d)]*gfs[IDX4S(8,2,2,i2d)]-gfs[IDX4S(3,2,2,i2d)]*gfs[IDX4S(6,2,2,i2d)]; //-THETA -> PI-|THETA|.. Nxx_plus_2NGHOSTS1-NGHOSTS-1 OR Nxx_plus_2NGHOSTS1-NGHOSTS-2???
        gfs[IDX4S(8,1,Nxx_plus_2NGHOSTS1-2,i2d)] = gfs[IDX4S(7,2,2,i2d)]/gfs[IDX4S(2,2,2,i2d)]-gfs[IDX4S(3,2,2,i2d)]*gfs[IDX4S(6,2,2,i2d)]/gfs[IDX4S(2,2,2,i2d)];
    } //end for phi loop (i2d)

    // PHI < -PI, THETA IN THE GRID
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++) {
        gfs[IDX4S(7,1,i1d,1)] = gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-1)]; //-THETA -> PI-|THETA|.. Nxx_plus_2NGHOSTS1-NGHOSTS-1 OR Nxx_plus_2NGHOSTS1-NGHOSTS-2???
        gfs[IDX4S(8,1,i1d,1)] = gfs[IDX4S(7,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-1)]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-1)]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-1)]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-1)];
    } //end for phi loop (i1d)

    // PHI > PI, THETA IN THE GRID
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++) {
        gfs[IDX4S(7,1,i1d,Nxx_plus_2NGHOSTS1-2)] = gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2)]; //-THETA -> PI-|THETA|.. Nxx_plus_2NGHOSTS1-NGHOSTS-1 OR Nxx_plus_2NGHOSTS1-NGHOSTS-2???
        gfs[IDX4S(8,1,i1d,Nxx_plus_2NGHOSTS1-2)] = gfs[IDX4S(7,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2)]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2)]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2)]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2)];
    } //end for phi loop (i1d)

    // MISSING {PHI<-PI,THETA<0}, {PHI<-PI,THETA>PI}, {PHI>PI,THETA<0}, {PHI>PI,THETA>0}
    */



    // SECOND LAYER OF GHOST POINTS
    /*for(int pt=0;pt<bcstruct->num_ib_gz_pts[1];pt++) {

    const int i0d2 = bcstruct->inner[1][pt].inner_bc_dest_pt.i0;
    const int i1d2 = bcstruct->inner[1][pt].inner_bc_dest_pt.i1;
    const int i2d2 = bcstruct->inner[1][pt].inner_bc_dest_pt.i2;
    const int i0s2  = bcstruct->inner[1][pt].inner_bc_src_pt.i0;
    const int i1s2  = bcstruct->inner[1][pt].inner_bc_src_pt.i1;
    const int i2s2  = bcstruct->inner[1][pt].inner_bc_src_pt.i2;

    // THE FOLLOWING ARE THE BCs FOR THE CHARACTERISTIC VARIABLES
    gfs[IDX4S(7,i0d2,i1d2,i2d2)] = gfs[IDX4S(2,i0s2,i1s2,i2s2)]*gfs[IDX4S(8,i0s2,i1s2,i2s2)]-gfs[IDX4S(3,i0s2,i1s2,i2s2)]*gfs[IDX4S(6,i0s2,i1s2,i2s2)];
    gfs[IDX4S(8,i0d2,i1d2,i2d2)] = gfs[IDX4S(7,i0s2,i1s2,i2s2)]/gfs[IDX4S(2,i0s2,i1s2,i2s2)]-gfs[IDX4S(3,i0s2,i1s2,i2s2)]*gfs[IDX4S(6,i0s2,i1s2,i2s2)]/gfs[IDX4S(2,i0s2,i1s2,i2s2)];


    // THE FOLLOWING ARE THE BCs FOR THE GFs INTRODUCED FOR THE EVANS TRICK
    gfs[IDX4S(3,i0d2,i1d2,i2d2)] = -gfs[IDX4S(3,i0s2,i1s2,i2s2)]; // CHIP
    gfs[IDX4S(0,i0d2,i1d2,i2d2)] = -gfs[IDX4S(0,i0s2,i1s2,i2s2)]; // RADIUS
    gfs[IDX4S(11,i0d2,i1d2,i2d2)] = -gfs[IDX4S(11,i0s2,i1s2,i2s2)]; // SINTH .. CHECK

    } */ // END for(int pt=0;pt<bcstruct->num_ib_gz_pts[which_gz];pt++)

    // THE FOLLOWING FILLS R NEGATIVE BUT THETA,PHI IN THE GRID (NOT THETA OR PHI GHOSTS)
    //for(int i1d=0; i1d<Nxx_plus_2NGHOSTS1; i1d++){
        //for(int i2d=0; i2d<Nxx_plus_2NGHOSTS2; i2d++){

            //THE FOLLOWING BCs SEEM TO WORK IN AXISYMMETRY BUT NOT IN FULL 3D
           /* gfs[IDX4S(11,0,i1d,i2d)] = gfs[IDX4S(4,3,i1d,i2d)]*gfs[IDX4S(18,3,i1d,i2d)]-gfs[IDX4S(5,3,i1d,i2d)]*gfs[IDX4S(10,3,i1d,i2d)];  //G+(-r)=chi(r)G-(r)-chi'(r)G(r)
            gfs[IDX4S(18,0,i1d,i2d)] = gfs[IDX4S(11,3,i1d,i2d)]/gfs[IDX4S(4,3,i1d,i2d)]-gfs[IDX4S(5,3,i1d,i2d)]*gfs[IDX4S(10,3,i1d,i2d)]/gfs[IDX4S(4,3,i1d,i2d)]; //G-(-r)=G+(r)/chi(r)-chi'(r)G(r)/chi(r)

            gfs[IDX4S(3,0,i1d,i2d)] = gfs[IDX4S(4,3,i1d,i2d)]*gfs[IDX4S(12,3,i1d,i2d)]-gfs[IDX4S(5,3,i1d,i2d)]*gfs[IDX4S(2,3,i1d,i2d)];
            gfs[IDX4S(12,0,i1d,i2d)] = gfs[IDX4S(3,3,i1d,i2d)]/gfs[IDX4S(4,3,i1d,i2d)]-gfs[IDX4S(5,3,i1d,i2d)]*gfs[IDX4S(2,3,i1d,i2d)]/gfs[IDX4S(4,3,i1d,i2d)];

            gfs[IDX4S(26,0,i1d,i2d)] = gfs[IDX4S(4,3,i1d,i2d)]*gfs[IDX4S(21,3,i1d,i2d)]-gfs[IDX4S(5,3,i1d,i2d)]*gfs[IDX4S(27,3,i1d,i2d)];
            gfs[IDX4S(21,0,i1d,i2d)] = gfs[IDX4S(26,3,i1d,i2d)]/gfs[IDX4S(4,3,i1d,i2d)]-gfs[IDX4S(5,3,i1d,i2d)]*gfs[IDX4S(27,3,i1d,i2d)]/gfs[IDX4S(4,3,i1d,i2d)];

            gfs[IDX4S(9,0,i1d,i2d)] = gfs[IDX4S(4,3,i1d,i2d)]*gfs[IDX4S(15,3,i1d,i2d)]-gfs[IDX4S(5,3,i1d,i2d)]*gfs[IDX4S(8,3,i1d,i2d)];
            gfs[IDX4S(15,0,i1d,i2d)] = gfs[IDX4S(9,3,i1d,i2d)]/gfs[IDX4S(4,3,i1d,i2d)]-gfs[IDX4S(5,3,i1d,i2d)]*gfs[IDX4S(8,3,i1d,i2d)]/gfs[IDX4S(4,3,i1d,i2d)]; */


            /*temp = gfs[IDX4S(2,3,i1d,i2d)]*gfs[IDX4S(8,0,i1d,i2d)]-gfs[IDX4S(3,3,i1d,i2d)]*gfs[IDX4S(6,0,i1d,i2d)];
            gfs[IDX4S(8,0,i1d,i2d)] = gfs[IDX4S(7,0,i1d,i2d)]/gfs[IDX4S(2,3,i1d,i2d)]-gfs[IDX4S(3,3,i1d,i2d)]*gfs[IDX4S(6,0,i1d,i2d)]/gfs[IDX4S(2,3,i1d,i2d)];
            gfs[IDX4S(7,0,i1d,i2d)] = temp;*/

            /*if(i2d<=(Nxx_plus_2NGHOSTS2/2)){
                gfs[IDX4S(7,0,i1d,i2d)] = gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]-gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]*gfs[IDX4S(6,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))];     // theta --> pi-theta, phi --> phi+theta
                gfs[IDX4S(8,0,i1d,i2d)] = gfs[IDX4S(7,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]/gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]-gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]*gfs[IDX4S(6,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]/gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d+((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))];
            }
            if(i2d>(Nxx_plus_2NGHOSTS2/2)){
                gfs[IDX4S(7,0,i1d,i2d)] = gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]*gfs[IDX4S(8,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]-gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]*gfs[IDX4S(6,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))];     // theta --> pi-theta, phi --> phi+theta
                gfs[IDX4S(8,0,i1d,i2d)] = gfs[IDX4S(7,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]/gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]-gfs[IDX4S(3,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]*gfs[IDX4S(6,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))]/gfs[IDX4S(2,3,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,i2d-((Nxx_plus_2NGHOSTS2-2*NGHOSTS)/2))];
            }*/
     //   } // END for i2d
   // }   // END for i1d

    /*
    //NEGATIVE THETA, PHI IN THE GRID ... not sure if its OK (check phi)
    for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++) {
        gfs[IDX4S(7,1,0,i2d)] = gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]; //-THETA -> PI-|THETA|.. Nxx_plus_2NGHOSTS1-NGHOSTS-1 OR Nxx_plus_2NGHOSTS1-NGHOSTS-2???
        gfs[IDX4S(8,1,0,i2d)] = gfs[IDX4S(7,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-2,i2d)]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-1,i2d)];
    } //end for phi loop (i2d)

    // THETA > PI, PHI IN THE GRID ... not sure if its OK (check phi)
    for(int i2d=NGHOSTS; i2d<Nxx_plus_2NGHOSTS2-NGHOSTS; i2d++) {
        gfs[IDX4S(7,1,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(2,2,3,i2d)]*gfs[IDX4S(8,2,3,i2d)]-gfs[IDX4S(3,2,3,i2d)]*gfs[IDX4S(6,2,3,i2d)]; //-THETA -> PI-|THETA|.. Nxx_plus_2NGHOSTS1-NGHOSTS-1 OR Nxx_plus_2NGHOSTS1-NGHOSTS-2???
        gfs[IDX4S(8,1,Nxx_plus_2NGHOSTS1-1,i2d)] = gfs[IDX4S(7,2,3,i2d)]/gfs[IDX4S(2,2,3,i2d)]-gfs[IDX4S(3,2,3,i2d)]*gfs[IDX4S(6,2,3,i2d)]/gfs[IDX4S(2,2,3,i2d)];
    } //end for phi loop (i2d)

    // PHI < -PI, THETA IN THE GRID
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++) {
        gfs[IDX4S(7,1,i1d,0)] = gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-2)]; //-THETA -> PI-|THETA|.. Nxx_plus_2NGHOSTS1-NGHOSTS-1 OR Nxx_plus_2NGHOSTS1-NGHOSTS-2???
        gfs[IDX4S(8,1,i1d,0)] = gfs[IDX4S(7,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-2)]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-2)]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-2)]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2-2)];
    } //end for phi loop (i1d)

    // PHI > PI, THETA IN THE GRID
    for(int i1d=NGHOSTS; i1d<Nxx_plus_2NGHOSTS1-NGHOSTS; i1d++) {
        gfs[IDX4S(7,1,i1d,Nxx_plus_2NGHOSTS1-1)] = gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(8,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2+1)]; //-THETA -> PI-|THETA|.. Nxx_plus_2NGHOSTS1-NGHOSTS-1 OR Nxx_plus_2NGHOSTS1-NGHOSTS-2???
        gfs[IDX4S(8,1,i1d,Nxx_plus_2NGHOSTS1-1)] = gfs[IDX4S(7,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2+1)]-gfs[IDX4S(3,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2+1)]*gfs[IDX4S(6,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2+1)]/gfs[IDX4S(2,2,Nxx_plus_2NGHOSTS1-NGHOSTS-i1d,Nxx_plus_2NGHOSTS2/2+1)];
    } //end for phi loop (i1d)

    // MISSING {PHI<-PI,THETA<0}, {PHI<-PI,THETA>PI}, {PHI>PI,THETA<0}, {PHI>PI,THETA>0}
    */



} // END function
