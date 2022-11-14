function [ fv, normals, name ] = read_stl( file, A, v )
% Modified by William Warriner 22 Aug 2019
% from READ_stl by Adam H. Aitkenhead

% Inputs:
%  - @file, a string containing the path to the desired STL file.
%  - @format, optional, default "auto", a string containing one of "ascii",
% "auto", or "binary". The string "auto" is preferred and allows selection
% of either "ascii" or "binary" from the input @file. If "ascii" or
% "binary" are supplied, that method of reading the STL file is forced,
% which may fail if the supplied value is incorrect.

% Outputs:
%  - @fv, a MATLAB FV-struct representing the data held in @file.
%  - @normals, the normals of the faces represented in @fv as a real,
% finite, 2D double array. Each row represents one normal (and one face of
% @fv, and each column one axis. Not recommended, use compute_normals()
% instead, and fix_vertex_ordering() as needed.
%  - @name, the name stored in the ASCII representation of @file, if
% applicable. Not recommended, though the file name itself may be more
% useful and relevant.

%{
% STL ASCII FILE FORMAT
%
% ASCII STL files have the following structure.  Technically each facet
% could be any 2D shape, but in practice only triangular facets tend to be
% used.  The present code ONLY works for meshes composed of triangular
% facets.
%
% solid object_name
% facet normal x y z
%   outer loop
%     vertex x y z
%     vertex x y z
%     vertex x y z
%   endloop
% endfacet
%
% <Repeat for all facets...>
%
% endsolid object_name
%}

%{
% STL BINARY FILE FORMAT
%
% Binary STL files have an 84 byte header followed by 50-byte records, each
% describing a single facet of the mesh.  Technically each facet could be
% any 2D shape, but that would screw up the 50-byte-per-facet structure, so
% in practice only triangular facets are used.  The present code ONLY works
% for meshes composed of triangular facets.
%
% HEADER:
% 80 bytes:  Header text
% 4 bytes:   (int) The number of facets in the STL mesh
%
% DATA:
% 4 bytes:  (float) normal x
% 4 bytes:  (float) normal y
% 4 bytes:  (float) normal z
% 4 bytes:  (float) vertex1 x
% 4 bytes:  (float) vertex1 y
% 4 bytes:  (float) vertex1 z
% 4 bytes:  (float) vertex2 x
% 4 bytes:  (float) vertex2 y
% 4 bytes:  (float) vertex2 z
% 4 bytes:  (float) vertex3 x
% 4 bytes:  (float) vertex3 y
% 4 bytes:  (float) vertex3 z
% 2 bytes:  Padding to make the data for each facet 50-bytes in length
%   ...and repeat for next facet...
%}

if ischar( file )
    file = string( file );
end
assert( isstring( file ) );
assert( isscalar( file ) );
assert( isfile( file ) );

format = identify_format( file );

if strcmpi( format, "ascii" )
    [ mesh, normals, name ] = read_ascii( file );
elseif strcmpi( format, "binary" )
    [ mesh, normals ] = read_binary( file );
    name = "unknown";
else
    assert( false );
end

%% This section written by Anders K. Hansen for MCmatlab integration
%======================================================
% SCALE, ROTATE AND/OR MIRROR WITH THE PROVIDED MULTIPLICATION MATRIX AND
% TRANSLATE WITH THE PROVIDED VECTOR
%======================================================

for iFacet = 1:size(mesh,1)
  for iVertex = 1:size(mesh,3)
    mesh(iFacet,:,iVertex) = A*mesh(iFacet,:,iVertex).' + v;
  end
end

h_f = figure(31);
if any(strcmp(h_f.UserData,file))
  clf reset;
end
if iscell(h_f.UserData)
  h_f.UserData = [h_f.UserData, file];
else
  h_f.UserData = file;
end
patch(squeeze(mesh(:,1,:)).',...
      squeeze(mesh(:,2,:)).',...
      squeeze(mesh(:,3,:)).',...
      [0 0 0],'FaceColor','none');
axis equal tight;
grid on; grid minor;
view(3);
set(gca,'ZDir','reverse');
xlabel('x [cm]');
ylabel('y [cm]');
zlabel('z [cm]');
set(h_f,'WindowStyle','Docked');
h_f.Name = 'STL wireframes';
title('STL wireframes');
if ~verLessThan('matlab','9.0')
  setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera');
end
rotate3d on

%%
fv = convert_triangle_geometry_format( mesh );

end


function format = identify_format( file )
% Test whether an stl file is ascii or binary
fid = fopen( file, "r" );
cleaner = onCleanup( @()fclose( fid ) );

% Binary files have size 84 + ( 50 * n )
fseek( fid, 0, "eof" );
file_size = ftell( fid );
if mod( file_size - 84, 50 ) ~= 0
    format = "ascii";
    return;
end

% Ascii files start with "solid", end with "endsolid <name>"
fseek( fid, 0, "bof" );
head = fread( fid, 80, "uchar" );
head = char( head.' );
head = strtrim( head );

fseek( fid, -80, "eof" );
tail = fread( fid, 80, "uchar" );
tail = char( tail.' );

HEAD_PATTERN = "solid";
TAIL_PATTERN = "endsolid";
if startsWith( head, HEAD_PATTERN ) && contains( tail, TAIL_PATTERN )
    format = "ascii";
    return;
else
    format = "binary";
    return;
end

end


function [ mesh, normals, name ] = read_ascii( file )

lines = read_contents( file );

% Read the STL name
% first line starts with "solid "
PREFIX = "solid";
name = lines( startsWith( lines, PREFIX ) );
if isempty( name )
    name = "unknown";
else
    name = strtrim( replace( name, PREFIX, "" ) );
end

% Read the vector normals
PREFIX = "facet normal";
normals = lines( startsWith( lines, PREFIX ) );
normals = replace( normals, PREFIX, "" );
normals = str2num( char( normals ) ); %#ok<ST2NM>

% Read the vertex coordinates
PREFIX = "vertex";
mesh = lines( startsWith( lines, PREFIX ) );
mesh = replace( mesh, PREFIX, "" );
mesh = str2num( char( mesh ) ); %#ok<ST2NM>
PREFIX = "endfacet";
face_count = numel( lines( startsWith( lines, PREFIX ) ) );
mesh = reshape( mesh, 3, face_count, 3 );
mesh = shiftdim( mesh, 1 );

end


function lines = read_contents( file )

fid = fopen( file, "r" );
cleaner = onCleanup( @()fclose( fid ) );
lines = textscan( fid, "%s", "delimiter", "\n" );
lines = string( lines{ 1 } );
lines = strtrim( lines );
lines = lines( lines ~= "" );

end


function [ mesh, normals ] = read_binary( file )

fid = fopen( file );
cleaner = onCleanup( @()fclose( fid ) );

fseek( fid, 80, 'bof' );
face_count = fread( fid, 1, "int32" );

values_per_face = 12;
value_count = face_count * values_per_face;
precision = sprintf( "%i*single=>single", values_per_face );
pad_byte_count = 2;
values = fread( fid, value_count, precision, pad_byte_count );
values = reshape( values, [ values_per_face face_count ] );
values = values.';

NORMAL = 1 : 3;
normals = values( :, NORMAL );
if any( 1 < abs( normals ), 'all' )
    warning( 'File may be corrupted: unexpected normal values.' );
end

VERTEX_1 = 4 : 6;
VERTEX_2 = 7 : 9;
VERTEX_3 = 10 : 12;
mesh = cat( ...
    3, ...
    values( :, VERTEX_1 ), ...
    values( :, VERTEX_2 ), ...
    values( :, VERTEX_3 ) ...
    );

end

