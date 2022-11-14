classdef Raster < handle
    % Raster is a class that implements the ray intersection method of
    % Patil S. and Ravi B [see Patil S and Ravi B.  Voxel-based
    % representation, display and thickness analysis of intricate shapes.
    % Ninth International Conference on Computer Aided Design and Computer
    % Graphics (CAD/CG 2005)].
    %
    % It is adapted and heavily modified from the original work
    % mesh_voxelisation by Adam H. Aitkenhead. [see 
    % fileexchange/27390-mesh-voxelisation].
    %
    % Inputs:
    %  - @grid, a Grid object.
    %  - @fv, a typical MATLAB-style FV-struct. Must be manifold and
    %  watertight to ensure proper operation. Some corrections exist, but
    %  should not be relied on.
    %  - @rays, a character array containing any subset of 'xyz'. The empty
    %  subset is replaced by 'xyz'.
    %
    % Properties:
    %  - @interior, a 3D logical array with dimensions M x N x P, where M is
    %  the number of elements along the grid x-axis, N along y, and P along
    %  Z. True values represent points in the interior of @fv.
    %  - @faces, a table with three columns. The column "indices"
    %  represents the linear indices into @interior where ray intersections
    %  occurred. The column "faces" represents the faces involved in the
    %  intersections. The column "dimension" represents the ray direction
    %  the intersection occurred along.
    %  - @normals, a table with four columns. The column "indices"
    %  represents the linear indices into @interior where the normal is
    %  given for faces in @faces. The columns "x", "y", and "z" represent
    %  the components of the normal vector along each axis.
    
    properties ( SetAccess = private )
        interior(:,:,:) logical
        faces table
    end
    
    properties ( SetAccess = private, Dependent )
        normals table
    end
    
    methods
      function obj = Raster( grid, fv, filename, A, v ) %, rays )
            assert( isa( grid, 'Grid' ) );
            assert( grid.dimension_count == 3 );
            
            assert( isstruct( fv ) );
            assert( isfield( fv, 'faces' ) );
            assert( isfield( fv, 'vertices' ) );
            
%             if isstring( rays )
%                 rays = char( rays );
%             end
%             if isempty( rays )
%                 rays = 'xyz';
%             end
%             assert( ischar( rays ) );
%             assert( all( ismember( numel( rays ), 1 : 3 ) ) );
%             assert( all( ismember( rays, 'xyz' ) ) );
            
            mesh = convert_triangle_geometry_format( fv );
%             dimensions = find( ismember( 'xyz', rays ) );

            dimensions = find( 1 < cellfun( @numel, grid.points ) );
            
            obj.grid = grid;
            obj.mesh = mesh;
            obj.dimensions = dimensions;
            
            obj.generate();
        end
        
        function value = get.normals( obj )
            if isempty( obj.faces )
                value = table();
                return;
            end
            n = compute_normals( obj.mesh );
            f = obj.faces;
            value = f( :, : );
            value.faces = [];
            value.dimension = [];
            value.x = n( f.faces, 1 );
            value.y = n( f.faces, 2 );
            value.z = n( f.faces, 3 );
            value = varfun( ...
                @sum, ...
                value, ...
                'groupingvariables', { 'indices' } ...
                );
            value.GroupCount = [];
            value.Properties.VariableNames{'sum_x'} = 'x';
            value.Properties.VariableNames{'sum_y'} = 'y';
            value.Properties.VariableNames{'sum_z'} = 'z';
            t = value{ :, { 'x' 'y' 'z' } };
            v = vecnorm( t, 2, 2 );
            t = t ./ v;
            value{ :, { 'x' 'y' 'z' } } = t;
        end
    end
    
    properties ( Access = private )
        grid Grid
        mesh(:,3,3) double {mustBeReal,mustBeFinite}
        dimensions(1,:) double
    end
    
    methods ( Access = private )
        function generate( obj )
            v = false( obj.grid.shape );
            f = table();
            for i = obj.dimensions
                c = obj.cast( i );
                p_inv = obj.get_inverse_permutation( i );
                v = v + permute( c.interior_array, p_inv );
                face_list = c.get_face_list( p_inv );
                if ~isempty( face_list )
                    face_list.dimension(:) = i;
                    f = [ f; face_list ]; %#ok<AGROW>
                end
            end
            v = v >= numel( obj.dimensions ) ./ 2;
            
            obj.interior = v;
            obj.faces = f;
        end
        
        function rr = cast( obj, dimension )
            p = obj.get_permutation( dimension );
            rr = RayRaster( ...
                obj.grid.points{ p }, ...
                obj.mesh( :, p, : ) ...
                );
        end
    end
    
    properties ( Access = private, Constant )
        PERMUTATIONS = [ ...
            2 3 1; ...
            3 1 2; ...
            1 2 3 ...
            ];
        INVERSE_PERMUTATIONS = [ ...
            3 1 2; ...
            2 3 1; ...
            1 2 3 ...
            ];
    end
    
    methods ( Access = private, Static )
        function p = get_permutation( dimension )
            p = Raster.PERMUTATIONS( dimension, : );
        end
        
        function p = get_inverse_permutation( dimension )
            p = Raster.INVERSE_PERMUTATIONS( dimension, : );
        end
    end
    
end

