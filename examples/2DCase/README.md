# Tutorial of library (2D case)

Here's an example to help users to get familiar with library (2D case). We would like to move the upper half portion of cylinder (radius = R) in clockwise sense and lower half portion in anticlockwise sense as shown in figure. 

![alt text](https://github.com/prstukumar/tangentVelocityAtPatches/blob/master/examples/2DCase/figure2.png)

This can be achieved by using this library on the cylinder wall in velocity file (U). Here's a case file attached which will simulate flow over cylinder (R = 0.5 m) with prescribed velocity boundary condition on cylinder wall as shown in figure. Boundary condition applied on cylinder in U file is shown here:

```
CYLINDER
{
    type              tangentVelocityAtPatches;
    axis              (0 0 1);
    rotatePoint       (0 0 0);
    velMagRotation    (-1 1);
    is2D              "yes";
    positions         ((0.5 0 0)(-0.5 0 0)(-0.5 0 0)(0.5 0 0));
 }
```
where, 
1. ```axis``` - axis of rotation
2. ```rotatePoint``` - This point must be inside body and if line drawn from this point, line must intersect body not more than two points.
3. ```velMagRotation``` - (-1 1) means upper patch is moving with velocity 1 unit in clockwise and lower patch with 1 unit velocity in anticlockwise direction
4. ```is2D``` - If body is 2D "yes", else "no"
5. ```positions``` - ((0.5 0 0)(-0.5 0 0)(-0.5 0 0)(0.5 0 0)) is a set of start and end coordinates of each patch, for upper patch start point is (0.5 0 0) while end point is (-0.5 0 0). Please note that start and end coordinates of each patch must be given in anticlockwise manner.

Also note for 3D case, user will have to provide face centroid indices to ```positions``` keyword as ```((first patch indices)(second patch indices) ... (nth patch indices))```. 

Use ```pisoFoam``` to run the case and visualize the simulation to see moving patches.
