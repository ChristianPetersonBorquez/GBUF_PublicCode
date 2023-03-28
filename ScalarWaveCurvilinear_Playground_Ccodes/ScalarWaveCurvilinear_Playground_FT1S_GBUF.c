
// Step P0: Define REAL and NGHOSTS; and declare CFL_FACTOR. This header is generated in NRPy+.
#include "ScalarWaveCurvilinear_Playground_REAL__NGHOSTS__CFL_FACTOR.h"

#include "rfm_files/rfm_struct__declare.h"

#include "declare_Cparameters_struct.h"

// All SIMD intrinsics used in SIMD-enabled C code loops are defined here:
#include "SIMD/SIMD_intrinsics.h"

// Step P1: Import needed header files
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "stdint.h" // Needed for Windows GCC 6.x compatibility
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.707106781186547524400844362104849039L
#endif

// Step P2: Declare the IDX4S(gf,i,j,k) macro, which enables us to store 4-dimensions of
//           data in a 1D array. In this case, consecutive values of "i"
//           (all other indices held to a fixed value) are consecutive in memory, where
//           consecutive values of "j" (fixing all other indices) are separated by
//           Nxx_plus_2NGHOSTS0 elements in memory. Similarly, consecutive values of
//           "k" are separated by Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1 in memory, etc.
#define IDX4S(g,i,j,k) \
( (i) + Nxx_plus_2NGHOSTS0 * ( (j) + Nxx_plus_2NGHOSTS1 * ( (k) + Nxx_plus_2NGHOSTS2 * (g) ) ) )
#define IDX3S(i,j,k) ( (i) + Nxx_plus_2NGHOSTS0 * ( (j) + Nxx_plus_2NGHOSTS1 * ( (k) ) ) )
#define LOOP_REGION(i0min,i0max, i1min,i1max, i2min,i2max) \
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++)
#define LOOP_ALL_GFS_GPS(ii) _Pragma("omp parallel for") \
  for(int (ii)=0;(ii)<Nxx_plus_2NGHOSTS_tot*NUM_EVOL_GFS;(ii)++)

// Step P3: Set UUGF and VVGF macros, as well as xx_to_Cart()
#include "boundary_conditions/gridfunction_defines.h"

// Step P4: Set xx_to_Cart(const paramstruct *restrict params,
//                     REAL *restrict xx[3],
//                     const int i0,const int i1,const int i2,
//                     REAL xCart[3]),
//           which maps xx->Cartesian via
//    {xx[0][i0],xx[1][i1],xx[2][i2]}->{xCart[0],xCart[1],xCart[2]}
#include "xx_to_Cart.h"

// Step P5: Defines set_Nxx_dxx_invdx_params__and__xx(const int EigenCoord, const int Nxx[3],
//                                       paramstruct *restrict params, REAL *restrict xx[3]),
//          which sets params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for
//          the chosen Eigen-CoordSystem if EigenCoord==1, or
//          CoordSystem if EigenCoord==0.
#include "set_Nxx_dxx_invdx_params__and__xx.h"

// Step P6: Include basic functions needed to impose curvilinear
//          parity and boundary conditions.
#include "boundary_conditions/CurviBC_include_Cfunctions.h"

// Step P7: Find the CFL-constrained timestep
#include "find_timestep.h"

// Part P8: Declare the function for the exact solution at a single point. time==0 corresponds to the initial data.
#include "exact_solution_single_point.h"

// Part P9: Declare the function for the exact solution at all points. time==0 corresponds to the initial data.
#include "exact_solution_all_points.h"

// Part P10: Declare the function to evaluate the scalar wave RHSs
#include "rhs_eval_GBUF_Evans_FINAL.h"

// Part P11: Declare the function to evaluate the constrants.... NOT SURE IF THIS GOES HERE
#include "constraint.h"

// main() function:
// Step 0: Read command-line input, set up grid structure, allocate memory for gridfunctions, set up coordinates
// Step 1: Set up scalar wave initial data
// Step 2: Output relative error between numerical and exact solution.
// Step 3: Evolve scalar wave initial data forward in time using Method of Lines with chosen RK-like algorithm,
//         applying quadratic extrapolation outer boundary conditions.
// Step 4: Free all allocated memory
int main(int argc, const char *argv[]) {
    paramstruct params;
#include "set_Cparameters_default.h"

    // Step 0a: Read command-line input, error out if nonconformant
    if(argc != 5 || atoi(argv[1]) < NGHOSTS || atoi(argv[2]) < NGHOSTS || atoi(argv[3]) < NGHOSTS) {
        printf("Error: Expected one command-line argument: ./ScalarWaveCurvilinear_Playground Nx0 Nx1 Nx2 runtag,\n");
        printf("where Nx[0,1,2] is the number of grid points in the 0, 1, and 2 directions.\n");
        printf("Nx[] MUST BE larger than NGHOSTS (= %d) and runtag must be a string.\n",NGHOSTS);
        exit(1);
    }
    // Step 0b: Set up numerical grid structure, first in space...
    const int Nxx[3] = { atoi(argv[1]), atoi(argv[2]), atoi(argv[3]) };
    if(Nxx[0]%2 != 0 || Nxx[1]%2 != 0 || Nxx[2]%2 != 0) {
        printf("Error: Cannot guarantee a proper cell-centered grid if number of grid cells not set to even number.\n");
        printf("       For example, in case of angular directions, proper symmetry zones will not exist.\n");
        exit(1);
    }

    // Step 0c: Set free parameters, overwriting Cparameters defaults
    //          by hand or with command-line input, as desired.
#include "free_parameters.h"

   // Step 0d: Uniform coordinate grids are stored to *xx[3]
    REAL *xx[3];
    // Step 0d.i: Set bcstruct
    bc_struct bcstruct;
    {
        int EigenCoord = 1;
        // Step 0d.ii: Call set_Nxx_dxx_invdx_params__and__xx(), which sets
        //             params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
        //             chosen Eigen-CoordSystem.
        set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &params, xx);
        // Step 0d.iii: Set Nxx_plus_2NGHOSTS_tot
#include "set_Cparameters-nopointer.h"
        const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;
        // Step 0e: Find ghostzone mappings; set up bcstruct
#include "boundary_conditions/driver_bcstruct.h"
        // Step 0e.i: Free allocated space for xx[][] array
        for(int i=0;i<3;i++) free(xx[i]);
    }

    // Step 0f: Call set_Nxx_dxx_invdx_params__and__xx(), which sets
    //          params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
    //          chosen (non-Eigen) CoordSystem.
    int EigenCoord = 0;
    set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &params, xx);

    // Step 0g: Set all C parameters "blah" for params.blah, including
    //          Nxx_plus_2NGHOSTS0 = params.Nxx_plus_2NGHOSTS0, etc.
#include "set_Cparameters-nopointer.h"
    const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

    // Step 0h: Time coordinate parameters
    const REAL t_final = 4.0*domain_size; //2*domain_size;

    // Step 0i: Set timestep based on smallest proper distance between gridpoints and CFL factor
    //REAL dt = 0.05/Nxx[0]/Nxx[1]; //Works unless suuper high resolution is used
    REAL dt = 0.5*find_timestep(&params, xx); //# original
    printf("# Timestep set to = %e\n",(double)dt);
    int N_final = (int)(t_final / dt + 0.5); // The number of points in time.
                                             // Add 0.5 to account for C rounding down
                                             // typecasts to integers.
    int output_every_N = (int)((REAL)N_final/100.0);
    if(output_every_N == 0) output_every_N = 1;

    // Step 0j: Error out if the number of auxiliary gridfunctions outnumber evolved gridfunctions.
    //              This is a limitation of the RK method. You are always welcome to declare & allocate
    //              additional gridfunctions by hand.
    if(NUM_AUX_GFS > NUM_EVOL_GFS) {
        printf("Error: NUM_AUX_GFS > NUM_EVOL_GFS. Either reduce the number of auxiliary gridfunctions,\n");
        printf("       or allocate (malloc) by hand storage for *diagnostic_output_gfs. \n");
        exit(1);
    }

    // Step 0k: Allocate memory for gridfunctions
#include "MoLtimestepping/RK_Allocate_Memory.h"

    // Step 0l: Set up precomputed reference metric arrays
    // Step 0l.i: Allocate space for precomputed reference metric arrays.
#include "rfm_files/rfm_struct__malloc.h"

    // Step 0l.ii: Define precomputed reference metric arrays.
    {
#include "set_Cparameters-nopointer.h"
#include "rfm_files/rfm_struct__define.h"
    }

    // Step 1: Set up initial data to be exact solution at time=0:
    params.time = 0.0; exact_solution_all_points(&params, xx, y_n_gfs);
    apply_bcs_curvilinear(&params, &bcstruct, NUM_EVOL_GFS, evol_gf_parity, y_n_gfs);

    // Label tag of simulation:
    char tag[20];
    sprintf(tag,"%s",argv[4]);
    printf("Tag for output is: %s, final time is: %e.\n",argv[4],t_final);

    // Write grid to gridname, so it does not need to be included in the time evolution output.
    char gridname[100];
    sprintf(gridname,"grid%s%d.txt",tag,Nxx[0]);
    FILE *grid3D = fopen(gridname, "w");
    LOOP_REGION(NGHOSTS,Nxx_plus_2NGHOSTS0-NGHOSTS,
                NGHOSTS,Nxx_plus_2NGHOSTS1-NGHOSTS,
                NGHOSTS,Nxx_plus_2NGHOSTS2-NGHOSTS) {
        const int idx = IDX3S(i0,i1,i2);
        REAL xx0 = xx[0][i0];
        REAL xx1 = xx[1][i1];
        REAL xx2 = xx[2][i2];
        REAL xCart[3];
        xx_to_Cart(&params,xx,i0,i1,i2,xCart);
        fprintf(grid3D,"%e %e %e\n",
            //        xCart[0],xCart[1],xCart[2], // outputs Cartesian gridpoints
                xx0,xx1,xx2);
        }
        fclose(grid3D);

    for(int n=0;n<=N_final;n++)
       { // Main loop to progress forward in time.

        // Step 1a: Set current time to correct value & compute exact solution
        params.time = ((REAL)n)*dt;

        // Step 1b: Calculate the constraints ..... NOW I THINK THIS DOESNT GO HERE
        //#include "constraint.h"
        //constraint(&rfmstruct, &params, y_n_gfs, diagnostic_output_gfs);

        // Step 2: Code validation: Compute log of L2 norm of difference
        //         between numerical and exact solutions:
        //   log_L2_Norm = log10( sqrt[Integral( [numerical - exact]^2 * dV)] ),
        //         where integral is within 30% of the grid outer boundary (domain_size)
        if(n%output_every_N == 0) {
            printf("%e\n",(double)params.time);
            REAL integral = 0.0;
            REAL numpts   = 0.0;
            constraint(&rfmstruct, &params, y_n_gfs, diagnostic_output_gfs);
#pragma omp parallel for reduction(+:integral,numpts)
            LOOP_REGION(NGHOSTS,Nxx_plus_2NGHOSTS0-NGHOSTS,
                        NGHOSTS,Nxx_plus_2NGHOSTS1-NGHOSTS,
                        NGHOSTS,Nxx_plus_2NGHOSTS2-NGHOSTS) {
                REAL xCart[3]; xx_to_Cart(&params,xx,i0,i1,i2, xCart);
                if(sqrt(xCart[0]*xCart[0] + xCart[1]*xCart[1] + xCart[2]*xCart[2]) < domain_size*0.3) {
                    REAL gg_exact,gp_exact,
                    bb_exact,bp_exact,uu_exact,up_exact,fgauge_exact,fgaugep_exact,
                    chi_exact,chip_exact,chipp_exact,rcoord_exact,Radius_exact,Rprime_exact,sinth_exact,costh_exact,
                    psigD0_exact,psigD1_exact,psigD2_exact,psibD0_exact,psibD1_exact,psibD2_exact,
                    psiuD0_exact,psiuD1_exact,psiuD2_exact,psifD0_exact,psifD1_exact,psifD2_exact;
                    exact_solution_single_point(xx[0][i0],xx[1][i1],xx[2][i2],&params,
                    &gg_exact,&gp_exact,&bb_exact,&bp_exact, &uu_exact,&up_exact, &fgauge_exact,&fgaugep_exact,
                    &chi_exact,&chip_exact,&chipp_exact,&rcoord_exact,&Radius_exact,&Rprime_exact,&sinth_exact,&costh_exact,
                    &psigD0_exact,&psigD1_exact,&psigD2_exact,&psibD0_exact,&psibD1_exact,&psibD2_exact,
                    &psiuD0_exact,&psiuD1_exact,&psiuD2_exact,&psifD0_exact,&psifD1_exact,&psifD2_exact);
                    double num   = (double)y_n_gfs[IDX4S(GGGF,i0,i1,i2)];
                    double exact = (double)gg_exact;
                    integral += (num - exact)*(num - exact);
                    numpts   += 1.0;
                }
            }
            // Compute and output the log of the L2 norm.
//            REAL log_L2_Norm = log10(sqrt(integral/numpts));
//            printf("%e %e\n",(double)params.time,log_L2_Norm);

            // Write time evolution output of variables to files.
            char filename[100];
            sprintf(filename,"out%s%d-%08d.txt",tag,Nxx[0],n);
            FILE *out2D = fopen(filename, "w");
            fprintf(out2D,"%e\n",params.time);
            LOOP_REGION(NGHOSTS,Nxx_plus_2NGHOSTS0-NGHOSTS,
                        NGHOSTS,Nxx_plus_2NGHOSTS1-NGHOSTS,
                        NGHOSTS,Nxx_plus_2NGHOSTS2-NGHOSTS) {
                const int idx = IDX3S(i0,i1,i2);
                //REAL xx0 = xx[0][i0];
                //REAL xx1 = xx[1][i1];
                //REAL xx2 = xx[2][i2];
                //REAL xCart[3];
                //xx_to_Cart(&params,xx,i0,i1,i2,xCart);
                fprintf(out2D,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
                //        xCart[0],xCart[1],xCart[2], // outputs Cartesian gridpoints
                        //xx0,xx1,xx2, // outputs gridpoints in original coordinates
                        y_n_gfs[IDX4S(GGGF,i0,i1,i2)],
                        y_n_gfs[IDX4S(GPGF,i0,i1,i2)],
                        y_n_gfs[IDX4S(PSIGD0GF,i0,i1,i2)],
                        y_n_gfs[IDX4S(PSIGD1GF,i0,i1,i2)],
                        y_n_gfs[IDX4S(PSIGD2GF,i0,i1,i2)],
                        y_n_gfs[IDX4S(BBGF,i0,i1,i2)],
                        y_n_gfs[IDX4S(BPGF,i0,i1,i2)],
                        y_n_gfs[IDX4S(PSIBD0GF,i0,i1,i2)],
                        y_n_gfs[IDX4S(PSIBD1GF,i0,i1,i2)],
                        y_n_gfs[IDX4S(PSIBD2GF,i0,i1,i2)],
                        y_n_gfs[IDX4S(UUGF,i0,i1,i2)],
                        y_n_gfs[IDX4S(UPGF,i0,i1,i2)],
                        y_n_gfs[IDX4S(PSIUD0GF,i0,i1,i2)],
                        y_n_gfs[IDX4S(PSIUD1GF,i0,i1,i2)],
                        y_n_gfs[IDX4S(PSIUD2GF,i0,i1,i2)],
                        y_n_gfs[IDX4S(FGAUGEGF,i0,i1,i2)],
                        y_n_gfs[IDX4S(FGAUGEPGF,i0,i1,i2)],
                        y_n_gfs[IDX4S(PSIFD0GF,i0,i1,i2)],
                        y_n_gfs[IDX4S(PSIFD1GF,i0,i1,i2)],
                        y_n_gfs[IDX4S(PSIFD2GF,i0,i1,i2)],
                        diagnostic_output_gfs[IDX4S(CDG0GF,i0,i1,i2)],
                        diagnostic_output_gfs[IDX4S(CDG1GF,i0,i1,i2)],
                        diagnostic_output_gfs[IDX4S(CDG2GF,i0,i1,i2)],
                        diagnostic_output_gfs[IDX4S(CDB0GF,i0,i1,i2)],
                        diagnostic_output_gfs[IDX4S(CDB1GF,i0,i1,i2)],
                        diagnostic_output_gfs[IDX4S(CDB2GF,i0,i1,i2)],
                        diagnostic_output_gfs[IDX4S(CDU0GF,i0,i1,i2)],
                        diagnostic_output_gfs[IDX4S(CDU1GF,i0,i1,i2)],
                        diagnostic_output_gfs[IDX4S(CDU2GF,i0,i1,i2)],
                        diagnostic_output_gfs[IDX4S(CDF0GF,i0,i1,i2)],
                        diagnostic_output_gfs[IDX4S(CDF1GF,i0,i1,i2)],
                        diagnostic_output_gfs[IDX4S(CDF2GF,i0,i1,i2)]);
            }
            fclose(out2D);

        }

        // Step 3: Step forward one timestep (t -> t+dt) in time using
        //           chosen RK-like MoL timestepping algorithm
#include "MoLtimestepping/RK_MoL.h"

    } // End main loop to progress forward in time.
    printf("%d\n",Nxx_plus_2NGHOSTS1);
    printf("%d\n",Nxx_plus_2NGHOSTS2);
    printf("%d\n",Nxx[1]+NGHOSTS);
    //printf("%d",bcstruct->inner[1][0].inner_bc_dest_pt.i0);
    printf("%d\n",NGHOSTS);
    // Step 4: Free all allocated memory
#include "rfm_files/rfm_struct__freemem.h"
#include "boundary_conditions/bcstruct_freemem.h"
#include "MoLtimestepping/RK_Free_Memory.h"
    for(int i=0;i<3;i++) free(xx[i]);
    return 0;
}
