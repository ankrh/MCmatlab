classdef RayRaster < handle
    
    properties ( SetAccess = private, Dependent )
        interior_array(:,:,:) logical
    end
    
    methods
        function obj = RayRaster( X, Y, Z, mesh )
            shape = [ numel( X ) numel( Y ) numel( Z ) ];
            x_index_range = obj.compute_index_range( X, mesh( :, 1, : ) );
            y_index_range = obj.compute_index_range( Y, mesh( :, 2, : ) );
            z_value_range = obj.compute_value_range( mesh( :, 3, : ) );
            mesh_min = min( mesh, [], 3 );
            mesh_max = max( mesh, [], 3 );
            
            obj.shape = uint32( shape );
            obj.X = X;
            obj.Y = Y;
            obj.Z = Z;
            obj.mesh = mesh;
            obj.x_index_range = x_index_range;
            obj.y_index_range = y_index_range;
            obj.z_value_range = z_value_range;
            obj.mesh_min = mesh_min;
            obj.mesh_max = mesh_max;
            
            obj.prepare();
        end
        
        function value = get_face_list( obj, permutation )
            if isempty( obj.crossings )
                value = [];
                return;
            end
            
            t = [ ...
                obj.rays.xi( obj.crossings.rays ) ...
                obj.rays.yi( obj.crossings.rays ) ...
                obj.crossings.zi ...
                ];
            t = t( :, permutation );
            s = obj.shape( permutation );
            inds = sub2ind( s, t(:,1), t(:,2), t(:,3) );
            value = obj.crossings( :, 'faces' );
            value.indices = inds;
            value = sortrows( value, 'indices' );
            lhs = find( ~diff( value.indices ) );
            rhs = lhs + 1;
            value( unique( [ lhs; rhs ] ), : ) = [];
        end
        
        function value = get.interior_array( obj )
            if isempty( obj.rays ) || isempty( obj.crossings )
                value = false( obj.shape );
            else
                value = obj.construct_voxels();
                value = obj.correct_voxels( value );
            end
        end
    end
    
    properties ( Access = private )
        shape(1,3) uint32
        X(:,1) double
        Y(:,1) double
        Z(:,1) double
        mesh(:,3,3) double
        
        x_index_range(2,1) double
        y_index_range(2,1) double
        z_value_range(2,1) double
        mesh_min(:,3) double
        mesh_max(:,3) double
        rays table
        crossings table
    end
    
    methods ( Access = private )
        function prepare( obj )
            faces = obj.initialize_rays();
            if isempty( faces )
                return;
            end
            obj.initialize_crossings( faces );
            v_cross = obj.identify_vertex_crossings();
            f_cross = obj.identify_face_crossings();
            for i = 1 : height( obj.rays )
                f_cross{ i } = [ v_cross{ i }; f_cross{ i } ];
            end
            obj.initialize_crossings( f_cross );
            obj.compute_z_values();
        end
        
        function faces = initialize_rays( obj )
            ray_count = obj.compute_initial_ray_count();
            faces = cell( ray_count, 1 );
            xp = zeros( ray_count, 1 );
            xi = zeros( ray_count, 1 );
            yp = zeros( ray_count, 1 );
            yi = zeros( ray_count, 1 );
            i = 1;
            for y_index = obj.y_index_range( 1 ) : obj.y_index_range( 2 )
                y = obj.Y( y_index );
                cross_y = find( y <= obj.mesh_max( :, 2 ) & obj.mesh_min( :, 2 ) <= y );
                for x_index = obj.x_index_range( 1 ) : obj.x_index_range( 2 )
                    x = obj.X( x_index );
                    check = x <= obj.mesh_max( cross_y, 1 ) & obj.mesh_min( cross_y, 1 ) <= x;
                    if any( check )
                        faces{ i } = cross_y( check );
                        xp( i ) = x;
                        xi( i ) = x_index;
                        yp( i ) = y;
                        yi( i ) = y_index;
                        i = i + 1;
                    end
                end
            end
            faces( i : end ) = [];
            xp( i : end ) = [];
            xi( i : end ) = [];
            yp( i : end ) = [];
            yi( i : end ) = [];
            correction = false( i - 1, 1 );
            if isempty( faces )
                obj.rays = table();
            else
                obj.rays = table( ...
                    xp, xi, yp, yi, correction, ...
                    'variablenames', { 'xp', 'xi', 'yp', 'yi', 'correction' } ...
                    );
            end
        end
        
        function initialize_crossings( obj, faces )
            face_counts = arrayfun( @(x)numel(x{1}), faces );
            face_indices = cell2mat( faces );
            ray_indices = zeros( numel( face_indices ), 1 );
            start = 1;
            for i = 1 : height( obj.rays )
                finish = face_counts( i ) + start - 1;
                ray_indices( start : finish ) = i;
                start = finish + 1;
            end
            z = nan( numel( face_indices ), 1 );
            obj.crossings = table( ...
                ray_indices, face_indices, z, ...
                'variablenames', { 'rays', 'faces', 'z' } ...
                );
        end
        
        function faces = identify_vertex_crossings( obj )
            faces = cell( height( obj.rays ), 1 );
            xpe = obj.rays.xp( obj.crossings.rays );
            ype = obj.rays.yp( obj.crossings.rays );
            check = ( obj.mesh(obj.crossings.faces,1,1)==xpe & obj.mesh(obj.crossings.faces,2,1)==ype ) ...
                | ( obj.mesh(obj.crossings.faces,1,2)==xpe & obj.mesh(obj.crossings.faces,2,2)==ype ) ...
                | ( obj.mesh(obj.crossings.faces,1,3)==xpe & obj.mesh(obj.crossings.faces,2,3)==ype );
            check = find( check );
            f = obj.crossings.faces;
            r = obj.crossings.rays;
            c = obj.rays.correction;
            for i = 1 : numel( check )
                check_index = check( i );
                [ correction, faces_crossed ] = ...
                    obj.find_vertex_crossings( f( check_index ) );
                ri = r( check_index );
                c( ri ) = correction;
                faces{ ri } = [ faces{ ri }; faces_crossed ];
            end
            obj.rays.correction = c;
            obj.crossings( obj.rays.correction, : ) = [];
        end
        
        function faces = identify_face_crossings( obj )
            inds = ( 1 : height( obj.crossings ) ).';
            rem = obj.crossings.faces;
            xx = obj.rays.xp( obj.crossings.rays );
            yy = obj.rays.yp( obj.crossings.rays );
            fy = obj.mesh(rem,2,2) - ((obj.mesh(rem,2,2)-obj.mesh(rem,2,3)) .* (obj.mesh(rem,1,2)-obj.mesh(rem,1,1))./(obj.mesh(rem,1,2)-obj.mesh(rem,1,3)));
            ry_slope = (obj.mesh(rem,1,2)-xx)./(obj.mesh(rem,1,2)-obj.mesh(rem,1,3));
            ry = obj.mesh(rem,2,2) - ((obj.mesh(rem,2,2)-obj.mesh(rem,2,3)) .* ry_slope);
            vert_check = isnan( ry_slope ) & ( ( xx >= obj.mesh(rem,1,2) & fy >= obj.mesh(rem,2,1) ) | ( xx <= obj.mesh(rem,1,2) & fy <= obj.mesh(rem,2,1) ) );
            check = (fy >= obj.mesh(rem,2,1) & ry >= yy) | (fy <= obj.mesh(rem,2,1) & ry <= yy) | vert_check;
            
            inds = inds( check );
            rem = rem( check );
            xx = xx( check );
            yy = yy( check );
            fy = obj.mesh(rem,2,3) - ((obj.mesh(rem,2,3)-obj.mesh(rem,2,1)) .* (obj.mesh(rem,1,3)-obj.mesh(rem,1,2))./(obj.mesh(rem,1,3)-obj.mesh(rem,1,1)));
            ry_slope = (obj.mesh(rem,1,3)-xx)./(obj.mesh(rem,1,3)-obj.mesh(rem,1,1));
            ry = obj.mesh(rem,2,3) - ((obj.mesh(rem,2,3)-obj.mesh(rem,2,1)) .* ry_slope);
            vert_check = isnan( ry_slope ) & ( ( xx >= obj.mesh(rem,1,3) & fy >= obj.mesh(rem,2,2) ) | ( xx <= obj.mesh(rem,1,3) & fy <= obj.mesh(rem,2,2) ) );
            check = (fy >= obj.mesh(rem,2,2) & ry >= yy) | (fy <= obj.mesh(rem,2,2) & ry <= yy) | vert_check;
            
            inds = inds( check );
            rem = rem( check );
            xx = xx( check );
            yy = yy( check );
            fy = obj.mesh(rem,2,1) - ((obj.mesh(rem,2,1)-obj.mesh(rem,2,2)) .* (obj.mesh(rem,1,1)-obj.mesh(rem,1,3))./(obj.mesh(rem,1,1)-obj.mesh(rem,1,2)));
            ry_slope = (obj.mesh(rem,1,1)-xx)./(obj.mesh(rem,1,1)-obj.mesh(rem,1,2));
            ry = obj.mesh(rem,2,1) - ((obj.mesh(rem,2,1)-obj.mesh(rem,2,2)) .* ry_slope);
            vert_check = isnan( ry_slope ) & ( ( xx >= obj.mesh(rem,1,1) & fy >= obj.mesh(rem,2,3) ) | ( xx <= obj.mesh(rem,1,1) & fy <= obj.mesh(rem,2,3) ) );
            check = (fy >= obj.mesh(rem,2,3) & ry >= yy) | (fy <= obj.mesh(rem,2,3) & ry <= yy) | vert_check;
            
            inds = inds( check );
            rem = rem( check );
            faces = cell( height( obj.rays ), 1 );
            r = obj.crossings.rays;
            for i = 1 : numel( inds )
                index = inds( i );
                ri = r( index );
                faces{ ri } = [ faces{ ri }; rem( i ) ];
            end
        end
        
        function compute_z_values( obj )
            faces = obj.crossings.faces;
            A = obj.mesh(faces,2,1) .* ( obj.mesh(faces,3,2) - obj.mesh(faces,3,3) ) ...
                + obj.mesh(faces,2,2) .* ( obj.mesh(faces,3,3) - obj.mesh(faces,3,1) ) ...
                + obj.mesh(faces,2,3) .* ( obj.mesh(faces,3,1) - obj.mesh(faces,3,2) );
            B = obj.mesh(faces,3,1) .* ( obj.mesh(faces,1,2) - obj.mesh(faces,1,3) ) ...
                + obj.mesh(faces,3,2) .* ( obj.mesh(faces,1,3) - obj.mesh(faces,1,1) ) ...
                + obj.mesh(faces,3,3) .* ( obj.mesh(faces,1,1) - obj.mesh(faces,1,2) );
            C = obj.mesh(faces,1,1) .* ( obj.mesh(faces,2,2) - obj.mesh(faces,2,3) ) ...
                + obj.mesh(faces,1,2) .* ( obj.mesh(faces,2,3) - obj.mesh(faces,2,1) ) ...
                + obj.mesh(faces,1,3) .* ( obj.mesh(faces,2,1) - obj.mesh(faces,2,2) );
            D = - obj.mesh(faces,1,1) .* ( obj.mesh(faces,2,2).*obj.mesh(faces,3,3) - obj.mesh(faces,2,3).*obj.mesh(faces,3,2) ) ...
                - obj.mesh(faces,1,2) .* ( obj.mesh(faces,2,3).*obj.mesh(faces,3,1) - obj.mesh(faces,2,1).*obj.mesh(faces,3,3) ) ...
                - obj.mesh(faces,1,3) .* ( obj.mesh(faces,2,1).*obj.mesh(faces,3,2) - obj.mesh(faces,2,2).*obj.mesh(faces,3,1) );
            
            C( abs( C ) < obj.TOL ) = 0;
            xpe = obj.rays.xp( obj.crossings.rays );
            ype = obj.rays.yp( obj.crossings.rays );
            z = -( D + A.*xpe + B.*ype ) ./ C;
            ROUND_TOL = 1e12;
            obj.crossings.z = round( z * ROUND_TOL ) / ROUND_TOL;
            check = obj.z_value_range( 1 ) <= obj.crossings.z & obj.crossings.z <= obj.z_value_range( 2 );
            obj.crossings = obj.crossings( check, : );
            if ~isempty( obj.crossings )
                [ ~, ia ] = unique( obj.crossings( :, { 'rays' 'z' } ), 'stable', 'rows' );
                obj.crossings = obj.crossings( ia, : );
            end
            %obj.crossings = unique( obj.crossings, 'stable', 'rows' );
        end
        
        function v = construct_voxels( obj )
            obj.crossings.xie = obj.rays.xi( obj.crossings.rays );
            obj.crossings.yie = obj.rays.yi( obj.crossings.rays );
            obj.crossings = sortrows( obj.crossings, { 'rays', 'xie', 'yie', 'z' } );
            v = false( obj.shape );
            fill_count = height( obj.crossings );
            i = 1;
            r = obj.crossings.rays;
            z = obj.crossings.z;
            xie = obj.crossings.xie;
            yie = obj.crossings.yie;
            c = obj.rays.correction;
            zi = zeros( height( obj.crossings ), 1 );
            while i <= fill_count
                if i == fill_count
                    [ ~, m ] = min( abs( obj.Z - z( i ) ) );
                    zi( i ) = m;
                    c( r( i ) ) = true;
                    break;
                elseif r( i ) == r( i + 1 )
                    inside = z( i ) < obj.Z & obj.Z < z( i + 1 );
                    if any( inside )
                        v( xie( i ), yie( i ), inside ) = true;
                        zi( i ) = find( inside, 1, 'first' );
                        zi( i + 1 ) = find( inside, 1, 'last' );
                    else
                        [ ~, m ] = min( abs( obj.Z - z( i ) ) );
                        zi( i ) = m;
                        [ ~, m ] = min( abs( obj.Z - z( i + 1 ) ) );
                        zi( i + 1 ) = m;
                        c( r( i ) ) = true;
                        c( r( i + 1 ) ) = true;
                    end
                    i = i + 2;
                else
                    [ ~, m ] = min( abs( obj.Z - z( i ) ) );
                    zi( i ) = m;
                    c( r( i ) ) = true;
                    i = i + 1;
                end
            end
            obj.crossings.rays = r;
            obj.crossings.z = z;
            obj.crossings.zi = zi;
            obj.crossings.xie = [];
            obj.crossings.yie = [];
            obj.rays.correction = c;
        end
        
        function v = correct_voxels( obj, v )
            correction_count = sum( obj.rays.correction );
            if correction_count <= 0
                return;
            end
            
            xi = obj.rays.xi;
            yi = obj.rays.yi;
            
            padSize = size(v) + [2 2 0];
            vpad = zeros(padSize);
            vpad(2:end-1,2:end-1,:) = v;
            v = vpad;

            xi = xi + 1;
            yi = yi + 1;
            corrections = obj.rays.correction;
            for i = 1 : correction_count
                c = corrections( i );
                xc = xi( c );
                yc = yi( c );
                vc = squeeze( sum( [ ...
                    v( xc - 1, yc - 1, : ),...
                    v( xc - 1, yc, : ),...
                    v( xc - 1, yc + 1, : ),...
                    v( xc, yc - 1, : ),...
                    v( xc, yc + 1, : ),...
                    v( xc + 1, yc - 1, : ),...
                    v( xc + 1, yc, : ),...
                    v( xc + 1, yc + 1, : ),...
                    ] ) );
                v( xc, yc, vc >= 4 ) = true;
            end
            v = v( 2 : end-1, 2 : end-1, : );
        end
        
        function [ correction, faces ] = find_vertex_crossings( obj, cross )
            correction = false;
            faces = [];
            if isempty( cross )
                return;
            end
            
            check = zeros( 1, numel( cross ) );
            while min( check ) == 0
                index = find( check == 0, 1, 'first' );
                check( index ) = 1;
                
                fv = convert_triangle_geometry_format( obj.mesh( cross, :, : ) );
                adjacent = ismember( fv.faces, fv.faces( index, : ) );
                adjacent = max( adjacent, [], 2 );
                check( adjacent ) = 1;
                
                normals = compute_normals( obj.mesh( cross( adjacent ), :, : ) );
                if max( normals( :, 3 ) ) < 0 || min( normals( :, 3 ) ) > 0
                    faces = [ faces cross( index ) ]; %#ok<AGROW>
                else
                    correction = true;
                    break;
                end
            end
        end
        
        function value = compute_initial_ray_count( obj )
            value = obj.compute_initial_x_count() ...
                * obj.compute_initial_y_count();
        end
        
        function value = compute_initial_x_count( obj )
            value = diff( obj.x_index_range ) + 1;
        end
        
        function value = compute_initial_y_count( obj )
            value = diff( obj.y_index_range ) + 1;
        end
    end
    
    properties ( Access = private, Constant )
        TOL = 1e-12;
    end
    
    methods ( Access = private, Static )
        function range = compute_index_range( x, mesh_x )
            r = [ min( mesh_x(:) ) max( mesh_x(:) ) ];
            r_min = abs( x - r( 1 ) );
            p_min = find( r_min == min( r_min ), 1, 'first' );
            r_max = abs( x - r( 2 ) );
            p_max = find( r_max == min( r_max ), 1, 'last' );
            range = sort( [ p_min p_max ] );
        end
        
        
        function range = compute_value_range( mesh_x )
            range = [ ...
                min( mesh_x(:) ) - RayRaster.TOL ...
                max( mesh_x(:) ) + RayRaster.TOL...
                ];
        end
    end
    
end
