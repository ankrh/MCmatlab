function varargout = convert_triangle_geometry_format( varargin )
% Modified by William Warriner 22 Aug 2019
% from CONVERT_meshformat by Adam H. Aitkenhead.

% Inputs:
% - @mesh, a real, finite, 3D double array with dimensions (:,3,3)
% representing face-vertex data with only triangular faces. The first
% dimension represents faces, the second represents dimensions (x,y,z), and
% the third represents vertices making up the faces.
% - OR -
% - @vertices, a real, finite, 2D double array with dimensions (:,3) whose
% rows represent vertices, and whose columns represent dimensions (x,y,z).
% - @faces, a real, finite, positive, 2D double array with dimensions (:,3)
% whose rows represent facets, and whose columns represent indices into the
% rows of @vertices.
% - OR -
% - @fv, a struct with fields faces and vertices, which have the same
% requirements as the inputs @vertices and @faces above.

% Outputs:
% - @mesh, as in the inputs. Output if the input is either the pair @faces
% and @vertices, or the singular @fv.
% - OR -
% - @fv, as in the inputs. Output if the input is @mesh.

assert( ismember( nargin, [ 1 2 ] ) );

first = varargin{ 1 };
if nargin == 1 && isstruct( first )
    faces = first.faces;
    vertices = first.vertices;
    mesh = [];
elseif nargin == 2
    faces = varargin{ 2 };
    vertices = first;
    mesh = [];
elseif nargin == 1
    faces = [];
    vertices = [];
    mesh = first;
else
    assert( false )
end

if isempty( mesh )
    mesh = zeros( size( faces, 1 ), 3, 3 );
    for i = 1 : size( faces, 1 )
        mesh( i, :, 1 ) = vertices( faces( i, 1 ), : );
        mesh( i, :, 2 ) = vertices( faces( i, 2 ), : );
        mesh( i, :, 3 ) = vertices( faces( i, 3 ), : );
    end
    varargout{ 1 } = mesh;
elseif isempty( faces ) && isempty( vertices )
    vertices = [ mesh( :, :, 1 ); mesh( :, :, 2 ); mesh( :, :, 3 ) ];
    vertices = unique( vertices, 'rows' );
    faces = zeros( size( mesh, 1 ), 3 );
    for face = 1 : size( mesh, 1 )
        for vertex = 1 : 3
            index = find( vertices( :, 1 ) == mesh( face, 1, vertex ) );
            index = index( vertices( index, 2 ) == mesh( face, 2, vertex ) );
            index = index( vertices( index, 3 ) == mesh( face, 3, vertex ) );
            faces( face, vertex ) = index;
        end
    end
    fv.faces = faces;
    fv.vertices = vertices;
    varargout{ 1 } = fv;
else
    assert( false )
end

end

