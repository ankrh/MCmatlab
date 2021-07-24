function inside = findInsideVoxels(X,Y,Z,filename,xyzRotAngles,translation)
inside = VOXELISE(X(1:end,1,1),Y(1,1:end,1),squeeze(Z(1,1,1:end)),filename,xyzRotAngles,translation,'xyz');
end