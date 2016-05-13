from steps.utilities.meshio import importAbaqus, stetmesh, saveMesh

MyMesh = importAbaqus('Cylinder1', 1e-7)[0]
print type(MyMesh)
saveMesh('testMeshCylinder1', MyMesh)
