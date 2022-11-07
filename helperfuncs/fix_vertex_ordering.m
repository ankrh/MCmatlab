function mesh = fix_vertex_ordering( mesh )
% Written by William Warriner 22 Aug 2019
% from COMPUTE_mesh_normals by Adam H. Aitkenhead

% Inputs:
%  - @mesh is either an FV-struct, or a real, finite, 3D double array with
% dimensions (:,3,3) representing face-vertex data with only triangular 
% faces. The first dimension represents faces, the second represents 
% dimensions (x,y,z), and the third represents vertices making up the 
% faces.
%
% Outputs:
%  - @mesh is the same as the input, except the vertex ordering is such
%  that the mesh has a consistent winding.

fv_out = false;
if isstruct( mesh )
    fv_out = true;
    mesh = convert_triangle_geometry_format( mesh );
end

face_count = size( mesh, 1 );
warning_done = false;
start = 1;
point_a = 1;
checked = false( face_count, 1 );
waiting = false( face_count, 1 );
while min( checked ) == 0
    checked( start ) = 1;
    point_b = point_a + 1;
    if point_b == 4
        point_b = 1;
    end
    
    % Find points which match point_a
    same_x = mesh( :, 1, : ) == mesh( start, 1, point_a );
    same_y = mesh( :, 2, : ) == mesh( start, 2, point_a );
    same_z = mesh( :, 3, : ) == mesh( start, 3, point_a );
    [ lhs, rhs ] = find( same_x & same_y & same_z );
    match_a = [ lhs, rhs ];
    match_a = match_a( match_a( :, 1 ) ~= start, : );
    
    % Find points which match edge_b
    same_x = mesh( :, 1, : ) == mesh( start, 1, point_b );
    same_y = mesh( :, 2, : ) == mesh( start, 2, point_b );
    same_z = mesh( :, 3, : ) == mesh( start, 3, point_b );
    [lhs,rhs] = find( same_x & same_y & same_z );
    match_b = [ lhs, rhs ];
    match_b = match_b( match_b( :, 1 ) ~= start, : );
    
    % Find edges which match both edge_a and edge_b -> adjacent edge
    [ member_a, member_b ] = ismember( match_a( :, 1 ), match_b( :, 1 ) );
    match_face = match_a( member_a, 1 );
    
    if numel( match_face ) ~= 1
        if ~warning_done
            warning( 'Mesh is not manifold.' )
            warning_done = true;
        end
    else
        match_a = match_a( member_a, 2 );
        match_b = match_b( member_b( member_a ), 2 );
        if ~checked( match_face ) && ~waiting( match_face )
            % Ensure adjacent edge traveled in opposite direction
            if match_b - match_a == 1 || match_b - match_a == -2
                [ ...
                    mesh( match_face, :, match_a ), ...
                    mesh( match_face, :, match_b ) ...
                    ] = deal( ...
                    mesh( match_face, :, match_b ), ...
                    mesh( match_face, :, match_a ) ...
                    );
            end
        end
    end
    waiting( match_face ) = true;
    
    if point_a < 3
        point_a = point_a + 1;
    elseif point_a == 3
        point_a = 1;
        checked( start ) = true;
        start = find( waiting & ~checked, 1, 'first' );
    end
end

if fv_out
    mesh = convert_triangle_geometry_format( mesh );
end

end

