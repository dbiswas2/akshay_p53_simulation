
import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as solvmod
import steps.utilities.meshio as  meshio
import numpy as np
import matplotlib.pyplot as plt


def gen_model():


    mdl = smodel.Model()

    X = smodel.Spec('X', mdl)
    Y = smodel.Spec('Y', mdl)

    vsysA = smodel.Volsys('vsysA', mdl)
    vsysB = smodel.Volsys('vsysB', mdl)

    dcst_U_A = 0.1e-9
    dcst_U_B = 0.05e-9
    dcst_V_A = 0.05e-9
    dcst_V_B = 0.1e-9

    diff_U_A = smodel.Diff('diff_U_A', vsysA, X, dcst_U_A)
    diff_U_B = smodel.Diff('diff_U_B', vsysB, X, dcst_U_B)
    diff_V_A = smodel.Diff('diff_V_A', vsysA, Y, dcst_V_A)
    diff_V_B = smodel.Diff('diff_V_B', vsysB, Y, dcst_V_B)


    return mdl

def gen_geom():

    #n = 5
    #print n

    mesh = meshio.loadMesh('testMeshCylinder1')[0]

    ntets = mesh.countTets()

    tets_compA = []
    tets_compB = []
    tris_compA = set()
    tris_compB = set()

    z_max = mesh.getBoundMax()[2]
    z_min = mesh.getBoundMin()[2]
    print z_max
    print z_min
    z_mid = z_min + (z_max - z_min) / 2.0

    for t in range(ntets):

        barycz = mesh.getTetBarycenter(t)[2]

        tris = mesh.getTetTriNeighb(t)

        if (barycz < z_mid):
            tets_compA.append(t)
            tris_compA.add(tris[0])
            tris_compA.add(tris[1])
            tris_compA.add(tris[2])
            tris_compA.add(tris[3])
        else:
            tets_compB.append(t)
            tris_compB.add(tris[0])
            tris_compB.add(tris[1])
            tris_compB.add(tris[2])
            tris_compB.add(tris[3])

    compA = sgeom.TmComp('compA', mesh, tets_compA)
    compB = sgeom.TmComp('compB', mesh, tets_compB)

    compA.addVolsys('vsysA')
    compB.addVolsys('vsysB')

    tris_DB = tris_compA.intersection(tris_compB)
    tris_DB = list(tris_DB)

    diffusion_boundary = sgeom.DiffBoundary('diffusion_boundary', mesh, tris_DB)

    return mesh, tets_compA, tets_compB

print gen_geom()

mdl = gen_model()

mesh, tets_compA, tets_compB = gen_geom()

rng = srng.create('mt19937', 256)
rng.initialize(654)

sim = solvmod.Tetexact(mdl, mesh, rng)
sim.reset()

tpnts = np.arange(0.0,0.101,0.001)
ntpnts = tpnts.shape[0]

ntets = mesh.countTets()
resU = np.zeros((ntpnts, ntets))
resV = np.zeros((ntpnts, ntets))

tetU = mesh.findTetByPoint([0, 0, -4.99e-6])
tetV = mesh.findTetByPoint([0, 0, 4.99e-6])

sim.setTetCount(tetU, 'U', np.array([1000]).astype(np.uint16)[0])

# sim.setTetCount(tetY, 'V', 500)

# sim.setDiffBoundaryDiffusionActive('diffusion_boundary', 'U',True)
# sim.setDiffBoundaryDiffusionActive('diffusion_boundary', 'V',True)










