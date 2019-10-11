# Tutorial for 3D case

Here flow inside a pipe with 45 degree cut as shown in figure is simulated, where first three rows of patches of curved surface 
from inlet side are moved tangentially with x-axis as axis of rotation.

![alt text](https://github.com/prstukumar/tangentVelocityAtPatches/blob/master/examples/3Dcase/figure1.png)

Boundary condition applied on walls is shown as follows (U file):

```
walls
{
    type            tangentVelocityAtPatches;
    axis            (1 0 0);
    rotatePoint     (0 0 0);
    velMagRotation  (1 1 1);
    is2D            "no";
    positions       ((0 1 2 3 4)(5 6 7 8 9)(10 11 12 13 14));
}
```
where, 

```positions``` keyword here is index of face centroids of each patches. To find the centroid indices of face centroids, user
can use following code:

```
label patchID = mesh.boundaryMesh().findPatchID("walls");
OFstream os("centFile.dat");
forAll(cPatch, i)
{
    os << mesh.Cf().boundaryField()[patchID][i] <<endl;
}
```
Above code will write the centroid values of "walls" boundary with index order in centFile.dat file. To run the case file, type
first ```blockMesh``` and then ```simpleFoam```. Visualize the movement of patches in paraview afterwards.

