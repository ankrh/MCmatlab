function inside = findInsideVoxels(X,Y,Z,filename,A,v)
if ~isequal(size(A),[3 3])
  error('The multiplication matrix A (argument 5) must have dimensions 3x3. If you don''t want any scaling, rotation or mirroring, you can input the identity matrix, eye(3).');
end
if numel(v) ~= 3
  error('The translation vector v (argument 6) must be either a 1x3 or 3x1 array. If you don''t want to translate the points, you can input [0 0 0].');
end

inside = VOXELISE(X(1:end,1,1),Y(1,1:end,1),squeeze(Z(1,1,1:end)),filename,A,v,'xyz');
end