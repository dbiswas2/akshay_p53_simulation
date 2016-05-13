
import steps.model as smodel
import steps.geom as swm
import steps.rng as srng
import steps.solver as ssolver
import numpy
import pylab

def gen_model():

    mdl = smodel.Model()

    X = smodel.Spec('X', mdl)
    Y = smodel.Spec('Y', mdl)

    vsysA = smodel.Volsys('vsysA', mdl)
    vsysB = smodel.Volsys('vsysB', mdl)

    diff_X_A = smodel.Diff('diff_X_A', vsysA, dcst = 0.1e-9)
    diff_X_B = smodel.Diff('diff_X_B', vsysB, dcst = 0.1e-9)
    diff_Y_A = smodel.Diff('diff_Y_A', vsysA, dcst = 0.1e-9)
    diff_Y_B = smodel.Diff('diff_Y_B', vsysB, dcst = 0.1e-9)

    return mdl


