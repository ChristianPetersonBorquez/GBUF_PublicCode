	GBUF 3D-code generated from NRPy+

To compile and run the C codes with the current settings: 

- In the main directory (Public_GBUF_3Dcode/) compile in the terminal with:
gcc -std=gnu99 -Ofast -fopenmp -march=native -funroll-loops ScalarWaveCurvilinear_Playground_Ccodes/ScalarWaveCurvilinear_Playground_FT1S_GBUF.c -o output3Dhypwave/ScalarWaveCurvilinear_Playground_FT1S_GBUF -lm

- This will create an executable in the folder output3Dhypwave. To run, move to that directory and type:
./ScalarWaveCurvilinear_Playground_FT1S_GBUF N_r N_theta N_phi output_name
where N_r, N_theta and N_phi are the desired gridpoints in the corresponding directions and output_name is a given name for the output files.


In order to re-build the C codes you just have to run the notebook Tutorial-Start_to_Finish-ScalarWaveCurvilinear_FT1S_GBUF. The Initial Data is changed in the file ScalarWave/InitialData_FT1S_GBUF.py
Once the C files are built, you have to copy the file Public_GBUF_3Dcode/rhs_eval_GBUF_Evans_FINAL.h to the directory ScalarWaveCurvilinear_Playground_Ccodes/ and add it to the main C code ScalarWaveCurvilinear_Playground_FT1S_GBUF.c instead of rhs_eval.h (this file was modified by hand). Then compile and run this C code as previously explained.
