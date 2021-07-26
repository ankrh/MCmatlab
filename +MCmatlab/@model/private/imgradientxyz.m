function [Gx, Gy, Gz] = imgradientxyz(I, method)
%IMGRADIENTXYZ Find the directional gradients of a volume.
%   [Gx, Gy, Gz] = IMGRADIENTXYZ(V) takes a grayscale or binary volume V
%   as input and returns the gradient along the X axis, Gx, the Y axis, Gy
%   and the Z axis, Gz. X axis points in the direction of increasing column
%   subscripts, Y axis points in the direction of increasing row subscripts
%   and Z axis points in the direction of increasing third dimension
%   subscripts. Gx, Gy and Gz are the same size as the input volume V.
%
%   [Gx, Gy, Gz] = IMGRADIENTXYZ(V, METHOD) calculates the directional
%   gradients of the volume V using the specified METHOD. Supported METHODs
%   are:
%
%       'sobel'                 : Sobel gradient operator (default)
%
%       'prewitt'               : Prewitt gradient operator
%
%       'central'               : Central difference gradient dV/dx = (V(x+1)- V(x-1))/ 2
%
%       'intermediate'          : Intermediate difference gradient dV/dx = V(x+1) - V(x)
%
%   Class Support
%   -------------
%   The input volume V can be numeric or logical three-dimensional matrix, and
%   it must be nonsparse. Gx, Gy and Gz are of class double, unless the
%   input volume V is of class single, in which case Gx, Gy and Gz will be of
%   class single.
%
%   Notes
%   -----
%   1. When applying the gradient operator at the boundaries of the volume,
%      values outside the bounds of the volume are assumed to equal the
%      nearest volume border value. This is similar to the 'replicate'
%      boundary option in IMFILTER.
%
%
%   Example
%   -------
%   % Computes 3-D directional image gradients using Sobel's gradient
%   % operator on an MRI volume.
%       volData = load('mri');
%       sz = volData.siz;
%       vol = squeeze(volData.D);
%
%       [Gx, Gy, Gz] = imgradientxyz(vol);
%
%       % Visualize directional gradients as a montage.
%       figure, montage(reshape(Gx,sz(1),sz(2),1,sz(3)),'DisplayRange',[])
%       title('Gradient magnitude along X')
%
%       figure, montage(reshape(Gy,sz(1),sz(2),1,sz(3)),'DisplayRange',[])
%       title('Gradient magnitude along Y')
%
%       figure, montage(reshape(Gz,sz(1),sz(2),1,sz(3)),'DisplayRange',[])
%       title('Gradient magnitude along Z')
%
%   See also imgradient3, imgradient, imgradientxy.

% Copyright 2015-2018 The MathWorks, Inc.

narginchk(1,2);

validateattributes(I,{'numeric','logical'},{'3d','nonsparse','real'}, ...
                   mfilename,'I',1);

if nargin==1
    method = 'sobel'; % Default method
end

if (nargin > 1)
    validateattributes(method,{'char','string'}, ...
        {'scalartext'},mfilename,'METHOD',2);
    methodstrings = {'sobel', 'prewitt', 'central', 'intermediate'};
    method = validatestring(method, methodstrings, ...
        mfilename, 'METHOD', 2);
end


if ~isfloat(I)
    I = double(I);
end

switch method
    case 'sobel'

        % 3-D Kernel for sobel along X, Y and Z direction
        hx(:,:,1) = [-1 0 1; -3 0 3; -1 0 1];
        hx(:,:,2) = [-3 0 3; -6 0 6; -3 0 3];
        hx(:,:,3) = [-1 0 1; -3 0 3; -1 0 1];

        hy(:,:,1) = [-1 -3 -1; 0 0 0; 1 3 1];
        hy(:,:,2) = [-3 -6 -3; 0 0 0; 3 6 3];
        hy(:,:,3) = [-1 -3 -1; 0 0 0; 1 3 1];

        hz(:,:,1) = [-1 -3 -1; -3 -6 -3; -1 -3 -1];
        hz(:,:,2) = [0 0 0; 0 0 0; 0 0 0];
        hz(:,:,3) = [1 3 1; 3 6 3; 1 3 1];

%         Gx = imfilter(I,hx,'replicate');
%         if nargout > 1
%             Gy = imfilter(I,hy,'replicate');
%         end
%         if nargout > 2
%             Gz = imfilter(I,hz,'replicate');
%         end
        
        % ================
        % Addition by ankrh, add replicated padding:
        I_pad = NaN(size(I)+[2 2 2]);
        
        I_pad([1 end],[1 end],[1 end]) = I([1 end],[1 end],[1 end]);

        I_pad(2:end-1,[1 end],[1 end]) = I(:,[1 end],[1 end]);
        I_pad([1 end],2:end-1,[1 end]) = I([1 end],:,[1 end]);
        I_pad([1 end],[1 end],2:end-1) = I([1 end],[1 end],:);

        I_pad([1 end],2:end-1,2:end-1) = I([1 end],:,:);
        I_pad(2:end-1,[1 end],2:end-1) = I(:,[1 end],:);
        I_pad(2:end-1,2:end-1,[1 end]) = I(:,:,[1 end]);
        
        I_pad(2:end-1,2:end-1,2:end-1) = I;
        
        % Run the filter assuming zero-boundaries:
        Gx = imfilter(I_pad,hx);
        Gy = imfilter(I_pad,hy);
        Gz = imfilter(I_pad,hz);
        
        % Remove padding again:
        Gx = Gx(2:end-1,2:end-1,2:end-1);
        Gy = Gy(2:end-1,2:end-1,2:end-1);
        Gz = Gz(2:end-1,2:end-1,2:end-1);
        % ================
        
    case 'prewitt'

        % 3-D Kernel for prewitt along X, Y and Z direction
        hx(:,:,1) = [-1 0 1; -1 0 1; -1 0 1];
        hx(:,:,2) = [-1 0 1; -1 0 1; -1 0 1];
        hx(:,:,3) = [-1 0 1; -1 0 1; -1 0 1];

        hy(:,:,1) = [-1 -1 -1; 0 0 0; 1 1 1];
        hy(:,:,2) = [-1 -1 -1; 0 0 0; 1 1 1];
        hy(:,:,3) = [-1 -1 -1; 0 0 0; 1 1 1];

        hz(:,:,1) = [-1 -1 -1; -1 -1 -1; -1 -1 -1];
        hz(:,:,2) = [0 0 0; 0 0 0; 0 0 0];
        hz(:,:,3) = [1 1 1; 1 1 1; 1 1 1];


        Gx = imfilter(I,hx,'replicate');
        if nargout == 2
            Gy = imfilter(I,hy,'replicate');
        elseif nargout == 3
            Gy = imfilter(I,hy,'replicate');
            Gz = imfilter(I,hz,'replicate');
        end

    case 'central'

        if isrow(I)
            Gx = gradient(I);
            if nargout > 1
                Gy = zeros(size(I),'like', I);
            end
            if nargout > 2
                Gz = zeros(size(I),'like', I);
            end

        elseif iscolumn(I)
            Gx = zeros(size(I),'like', I);
            if nargout > 1
                Gy = gradient(I);
            end
            if nargout > 2
                Gz = zeros(size(I),'like', I);
            end

        elseif ismatrix(I)
            [Gx, Gy] = gradient(I);
            Gz = zeros(size(I),'like', I);

        else
            [Gx, Gy, Gz] = gradient(I);
        end


    case 'intermediate'
        Gx = zeros(size(I), 'like', I);
        if (size(I,2) > 1)
            Gx(:,1:end-1, :) = diff(I,1,2);
        end

        if nargout >= 2
            Gy = zeros(size(I),'like', I);
            if (size(I,1) > 1)
                Gy(1:end-1,:, :) = diff(I,1,1);
            end
        end

        if nargout == 3
            Gz = zeros(size(I),'like', I);
            if (size(I,3) > 1)
                Gz(:, :, 1:end-1) = diff(I,1,3);
            end
        end

end

end
