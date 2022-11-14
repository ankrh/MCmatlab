function normals = compute_normals( mesh )
% Modified by William Warriner 22 Aug 2019
% from COMPUTE_mesh_normals by Adam H. Aitkenhead

% Inputs:
%  - @mesh is either an FV-struct, or a real, finite, 3D double array with
% dimensions (:,3,3) representing face-vertex data with only triangular 
% faces. The first dimension represents faces, the second represents 
% dimensions (x,y,z), and the third represents vertices making up the 
% faces.
%
% Outputs:
%  - @normals is a real, finite, double array with dimensions (:,3)
%  representing the normal vectors of the faces in @mesh.

if isstruct( mesh )
    mesh = convert_triangle_geometry_format( mesh );
end

e1 = mesh( :, :, 2 ) - mesh( :, :, 1 );
e2 = mesh( :, :, 3 ) - mesh( :, :, 1 );
normals = cross( e1, e2 );
normals = normals ./ vecnorm( normals, 2, 2 );

end


