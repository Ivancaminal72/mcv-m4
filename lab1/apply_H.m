function It=apply_H(I, H)
% I: image to be transformed
% H: homography [3x3 matrix]
% It: transformed image

% image size
[nr, nc, nchan] = size(I);

%image corners
c1 = [1 1 1]';
c2 = [nc 1 1]';
c3 = [1 nr 1]';
c4 = [nc nr 1]';
        
% transform corners using H
H_c1=H*c1; H_c1=H_c1/(H_c1(3));
H_c2=H*c2; H_c2=H_c2/(H_c2(3)); 
H_c3=H*c3; H_c3=H_c3/(H_c3(3)); 
H_c4=H*c4; H_c4=H_c4/(H_c4(3)); 
        
% compute extremal transformed corner coordinates
x_min = round(min([H_c1(1) H_c2(1) H_c3(1) H_c4(1)]));
x_max = round(max([H_c1(1) H_c2(1) H_c3(1) H_c4(1)]));
y_min = round(min([H_c1(2) H_c2(2) H_c3(2) H_c4(2)]));
y_max = round(max([H_c1(2) H_c2(2) H_c3(2) H_c4(2)]));


% create matrices of homogeneous coordinates
[X,Y] = meshgrid(x_min:x_max, y_min:y_max);
Hnc = x_max - x_min + 1;
Hnr = y_max - y_min + 1;
Z = ones(Hnr,Hnc);


% Matrix with all image points, to be transformed, in projective space
XYZ = [X(:) Y(:) Z(:)]';

% transform image
Hi = inv(H); 
Hi_XYZ = Hi * XYZ; %
HX = reshape(Hi_XYZ(1,:), Hnr, Hnc);
HY = reshape(Hi_XYZ(2,:), Hnr, Hnc);
HZ = reshape(Hi_XYZ(3,:), Hnr, Hnc);
HX = HX ./ HZ;
HY = HY ./ HZ;

%Transformed image
It = zeros(Hnr, Hnc, nchan);
for l=1:nchan
    It(:,:,l) = interp2(double(I(:,:,l)), HX, HY, 'linear', 0);
end
