import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par    # NRPy+: Parameter interface
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import grid as gri                # NRPy+: Functions having to do with numerical grids
#import ScalarWave.ScalarWaveCurvilinear_RHSs_MetricComponents_FT1S_CharacteristicFields  as swrhs
import ScalarWave.ScalarWaveCurvilinear_RHSs_MetricComponents_FT1S_Hyperboloidal  as swrhs

thismodule = __name__

def Constraint():

    DIM = 3

# Alex: this should register the quantity - it may have to be modified / completed
    global CD
    for i in range(len(gri.glb_gridfcs_list)):
            if "CD0" in gri.glb_gridfcs_list[i].name:
                return
    CD = ixp.register_gridfunctions_for_single_rank1("AUX", "CD")

    CD = ixp.zerorank1()


# Alex: here we give the expression for the constraint

    CD[0] = -swrhs.gp/swrhs.chi +( sp.sympify(1)/(sp.sympify(2)*swrhs.Rprime-sp.sympify(1)) )*(swrhs.psiD[0] -swrhs.chip*swrhs.gg/swrhs.chi +sp.sympify(2)*swrhs.gg_dD[0] )
    CD[1] = swrhs.psiD[1]-swrhs.gg_dD[1]
    CD[2] = swrhs.sinth*swrhs.psiD[2]-swrhs.gg_dD[2]
