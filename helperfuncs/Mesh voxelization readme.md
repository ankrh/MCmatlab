# mesh_voxelization
Forked from original code posted at  https://www.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation

Changes:
 - Primary functionality contained in Raster class.
   - To use, create a Grid object from three grid vectors (formerly gridX, gridY, gridZ).
   - Create an FV-struct, possibly from read_stl().
   - Choose ray directions, default is 'xyz'.
   - Create the Raster object as Raster(grid, fv, rays).
   - The previous gridOUTPUT is contained in the interior property of the Raster object.
   - There is also a faces property, which is a table of faces crossed by rays, and the voxel the crossing occurred in as a linear index.
   - There is also a normals property, which is a table like faces, but contains the x, y, and z coordinates of unit normal vectors for each face, as well as the ray direction the crossing occurred for.
   - The example has been improved to show off these new features.
   - Primary functionality is now about 2-3 times faster, but allocates 2-5 times as much memory.
 - read_stl() was improved.
   - Binary files are now fully buffered, and a vectorized approach is used using new fread() typecasting functionality.
   - ASCII files have improved vectorization, and use string arrays for improved code readability.
 - compute_normals() was vectorized, and the vertex winding fix was refactored into its own function, fix_vertex_ordering().
 - CONVERT_mesh_format() was renamed convert_triangle_geometry_format() and the interface was simplified.
   - Passing an FV-struct returns an Nx3x3 (like the old meshXYZ).
   - Passing an array of faces and array of vertices returns an Nx3x3.
   - Passing an Nx3x3 returns an FV-struct.

Original readme below.

- VOXELISE.m modified by William Warriner 19 July 2019 to fix edge intersection issue.

Polygon mesh voxelisation
=========================

Adam H. Aitkenhead
adam.aitkenhead@christie.nhs.uk
The Christie NHS Foundation Trust

Voxelise a closed (ie. watertight) triangular-polygon mesh.  The mesh can be in one of several formats:  in an STL file;  in a structure containing the faces and vertices data;  in three 3xN arrays containing the x,y,z coordinates;  or in a single Nx3x3 array defining the vertex coordinates for each of the N facets.
 
 
USAGE:
======

[gridOUTPUT,gridCOx,gridCOy,gridCOz] = VOXELISE(gridX,gridY,gridZ,STLin,raydirection)
  ...or...
[gridOUTPUT,gridCOx,gridCOy,gridCOz] = VOXELISE(gridX,gridY,gridZ,meshFV,raydirection)
  ...or...
[gridOUTPUT,gridCOx,gridCOy,gridCOz] = VOXELISE(gridX,gridY,gridZ,meshX,meshY,meshZ,raydirection)
  ...or...
[gridOUTPUT,gridCOx,gridCOy,gridCOz] = VOXELISE(gridX,gridY,gridZ,meshXYZ,raydirection)
 
 
INPUT PARAMETERS
================

    gridX   - Mandatory - 1xP array     - List of the grid X coordinates. 
                          OR an integer - Number of voxels in the grid in the X direction.

    gridY   - Mandatory - 1xQ array     - List of the grid Y coordinates.
                          OR an integer - Number of voxels in the grid in the Y direction.

    gridZ   - Mandatory - 1xR array     - List of the grid Z coordinates.
                          OR an integer - Number of voxels in the grid in the Z direction.

    STLin   - Optional - string      - Filename of the STL file.

    meshFV  - Optional - structure   - Structure containing the faces and vertices of the mesh, in the same format as that produced by the isosurface command.

    meshX   - Optional - 3xN array   - List of the mesh X coordinates for the 3 vertices of each of the N facets
    meshY   - Optional - 3xN array   - List of the mesh Y coordinates for the 3 vertices of each of the N facets
    meshZ   - Optional - 3xN array   - List of the mesh Z coordinates for the 3 vertices of each of the N facets

    meshXYZ - Optional - Nx3x3 array - The vertex positions for each facet, with:
                                       1 row for each facet
                                       3 columns for the x,y,z coordinates
                                       3 pages for the three vertices

    raydirection - Optional - String - Defines in which directions the
                                       ray-tracing is performed.  The
                                       default is 'xyz', which traces in
                                       the x,y,z directions and combines
                                       the results.
 
 
OUTPUT PARAMETERS
=================

    gridOUTPUT - Mandatory - PxQxR logical array - Voxelised data (1=>Inside the mesh, 0=>Outside the mesh)

    gridCOx    - Optional - 1xP array - List of the grid X coordinates.
    gridCOy    - Optional - 1xQ array - List of the grid Y coordinates.
    gridCOz    - Optional - 1xR array - List of the grid Z coordinates.
 
 
EXAMPLES
========

    To voxelise an STL file:
    >>  [gridOUTPUT] = VOXELISE(gridX,gridY,gridZ,STLin)

    To voxelise a mesh defined by a structure containing the faces and vertices:
    >>  [gridOUTPUT] = VOXELISE(gridX,gridY,gridZ,meshFV)

    To voxelise a mesh where the x,y,z coordinates are defined by three 3xN arrays:
    >>  [gridOUTPUT] = VOXELISE(gridX,gridY,gridZ,meshX,meshY,meshZ)

    To voxelise a mesh defined by a single Nx3x3 array:
    >>  [gridOUTPUT] = VOXELISE(gridX,gridY,gridZ,meshXYZ)

    To also output the lists of X,Y,Z coordinates:
    >>  [gridOUTPUT,gridCOx,gridCOy,gridCOz] = VOXELISE(gridX,gridY,gridZ,STLin)

    To use ray-tracing in only the z-direction:
    >>  [gridOUTPUT] = VOXELISE(gridX,gridY,gridZ,STLin,'z')
 
 
NOTES
=====

  - Defining raydirection='xyz' means that the mesh is ray-traced in each of the x,y,z directions, with the overall result being a combination of the result from each direction.  This gives the most reliable result at the expense of computation time.
  - Tracing in only one direction (eg. raydirection='z') is faster, but can potentially lead to artefacts if rays exactly cross a facet edge.
 
 
REFERENCES
==========

  - This code uses a ray intersection method similar to that described by:
    Patil S and Ravi B.  Voxel-based representation, display and thickness analysis of intricate shapes. Ninth International Conference on Computer Aided Design and Computer Graphics (CAD/CG 2005)
 
 