function I2=apply_H(I, H, size_mode)

% Apply the transformation H to the image I
%
% INPUTS:
%   I: image to be transformed (color or gray scale image)
%   H: 3x3 matrix that specifies the desired transformation
%   size_mode: variable that specifies the size of the output image.
%              There are two possible options:
%              'original': the size of the input image is maintained
%              'fit': the size of the output image is changed so that
%                     the whole transformed image is contained in the
%                     output
%
% OUTOUT:
%   I2: transformed image
%

% input image size
[nrows, ncols, nchan] = size(I);

% check input value
switch size_mode %la funcio switch mira tots els possibles casos d'un input de la funcio
    case 'fit'
        % the new image size should contain the whole transformed imag
        
        % image corners in homogeneous coordinates
        c1 = [1 1 1]';
        c2 = [ncols 1 1]';
        c3 = [1 nrows 1]';
        c4 = [ncols nrows 1]';
        
        % transform corners according to H
        % *** TO COMPLETE ***
        %     Compute the transformed homogeneous image corners according to H.
        %     Normalize the homogeneous coordinates so as the third
        %     coordinate is always 1.
        %     Call the transformed corners: Hc1, Hc2, Hc3, Hc4
        %
        Hc1=H*c1; Hc1=Hc1/(Hc1(3)); % transformed corner c1
        Hc2=H*c2; Hc2=Hc2/(Hc2(3)); % transformed corner c2
        Hc3=H*c3; Hc3=Hc3/(Hc3(3)); % transformed corner c3
        Hc4=H*c4; Hc4=Hc4/(Hc4(3)); % transformed corner c4
        %
        % *** ***
        
        % compute extremal transformed corner coordinates
        xmin = round(min([Hc1(1) Hc2(1) Hc3(1) Hc4(1)]));
        xmax = round(max([Hc1(1) Hc2(1) Hc3(1) Hc4(1)]));
        ymin = round(min([Hc1(2) Hc2(2) Hc3(2) Hc4(2)]));
        ymax = round(max([Hc1(2) Hc2(2) Hc3(2) Hc4(2)]));
    
    case 'original'
        % the image size is the same as the original one
        xmin = 1;
        xmax = ncols;
        ymin = 1;
        ymax = nrows;
    
    otherwise
        error('size_mode should be fit/original ');
 end

% create matrices of homogeneous coordinates
[X,Y] = meshgrid(xmin:xmax, ymin:ymax);
Hncols = xmax - xmin + 1;
Hnrows = ymax - ymin + 1;
Z = ones(Hnrows,Hncols);


% create a 3x(Hnrows*Hncols) matrix in order to transform all the 
% coordinate points with matrix operations (thus avoiding the use of a for loop)
XYZs = [X(:) Y(:) Z(:)]';

% transform image
Hi = inv(H); %tenim la inversa de la homografia que passem
HiXYZs = Hi * XYZs; %
HX = reshape(HiXYZs(1,:), Hnrows, Hncols);
HY = reshape(HiXYZs(2,:), Hnrows, Hncols);
HZ = reshape(HiXYZs(3,:), Hnrows, Hncols);
HX = HX ./ HZ;
HY = HY ./ HZ;
% reshape(X,M,N) returns the M-by-N matrix whose elements
% are taken columnwise from X.  An error results if X does
% not have M*N elements.

I2 = zeros(Hnrows, Hncols, nchan);
for c=1:nchan,
    I2(:,:,c) = interp2(double(I(:,:,c)), HX, HY, 'linear', 0);
end

