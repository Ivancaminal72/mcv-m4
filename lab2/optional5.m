%% 5. OPTIONAL: Calibration with a planar pattern
clear all;

%% Read template and images.
T     = imread('Data/calib/template.jpg');
I{1}  = imread('Data/calib/graffiti1.tif');
I{2}  = imread('Data/calib/graffiti2.tif');
I{3}  = imread('Data/calib/graffiti3.tif');
%I{4}  = imread('Data/calib/graffiti4.tif');
%I{5}  = imread('Data/calib/graffiti5.tif');
Tg = sum(double(T), 3) / 3 / 255;
Ig{1} = sum(double(I{1}), 3) / 3 / 255;
Ig{2} = sum(double(I{2}), 3) / 3 / 255;
Ig{3} = sum(double(I{3}), 3) / 3 / 255;

N = length(I);

%% Compute keypoints.
addpath('sift')
fprintf('Computing sift points in template... ');
[pointsT, descrT] = sift(Tg, 'Threshold', 0.05);
fprintf(' done\n');

points = cell(N,1);
descr = cell(N,1);
for i = 1:N
    fprintf('Computing sift points in image %d... ', i);
    [points{i}, descr{i}] = sift(Ig{i}, 'Threshold', 0.05);
    fprintf(' done\n');
end

%% Match and compute homographies.
H = cell(N,1);
for i = 1:N
    % Match against template descriptors.
    fprintf('Matching image %d... ', i);
    matches = siftmatch(descrT, descr{i});
    fprintf('done\n');

    % Fit homography and remove outliers
    x1 = [pointsT(1:2, matches(1, :)); ones(1, length(matches))];
    x2 = [points{i}(1:2, matches(2, :)); ones(1, length(matches))];
    H{i} = 0;
    [H{i}, inliers] =  ransac_homography_adaptive_loop((x1), (x2), 3, 1000);

    % Plot inliers.
    figure;
    plotmatches(Tg, Ig{i}, pointsT(1:2,:), points{i}(1:2,:), matches(:, inliers));

    % Play with the homography
    %vgg_gui_H(T, I{i}, H{i});
end

%% Compute the Image of the Absolute Conic
V = zeros(2 * N, 6);
for i = 1:N
    h = H{i};
    V(2*i - 1,:) = [h(1,1) * h(1,2), ...
                    h(1,1) * h(2,2) + h(2,1) * h(1,2), ...
                    h(1,1) * h(3,2) + h(3,1) * h(1,2), ...
                    h(2,1) * h(2,2), ...
                    h(2,1) * h(3,2) + h(3,1) * h(2,2), ...
                    h(3,1) * h(3,2)];
    V(2*i,:) = [h(1,1) * h(1,1), ...
                h(1,1) * h(2,1) + h(2,1) * h(1,1), ...
                h(1,1) * h(3,1) + h(3,1) * h(1,1) , ...
                h(2,1) * h(2,1), ...
                h(2,1) * h(3,1) + h(3,1) * h(2,1), ...
                h(3,1) * h(3,1)] ...
                -  [h(1,2) * h(1,2), ...
                h(1,2) * h(2,2) + h(2,2) * h(1,2), ...
                h(1,2) * h(3,2) + h(3,2) * h(1,2), ...
                h(2,2) * h(2,2), ...
                h(2,2) * h(3,2) + h(3,2) * h(2,2), ...
                h(3,2) * h(3,2)];
end
[U,S,UT]=svd(V);
w = [UT(:,6)];
W= [w(1),w(2),w(3); w(2),w(4),w(5); w(3),w(5),w(6)]./w(6);
 
 
%% Recover the camera calibration.

K_INV= chol(W);
K= inv(K_INV);
K=K/K(3,3);
    
% ToDo: in the report make some comments related to the obtained internal
%       camera parameters and also comment their relation to the image size

%% Compute camera position and orientation.
R = cell(N,1);
t = cell(N,1);
P = cell(N,1);
figure;hold;
for i = 1:N
    % ToDo: compute r1, r2, and t{i}
    h = H{i};
    r1 = K\h(:,1)
    r2 = K\h(:,2)
    t{i} = K\h(:,3)
    
    % Solve the scale ambiguity by forcing r1 and r2 to be unit vectors.
    s = sqrt(norm(r1) * norm(r2)) * sign(t{i}(3));
    r1 = r1 / s;
    r2 = r2 / s;
    t{i} = t{i} / s;
    R{i} = [r1, r2, cross(r1,r2)];
    
    % Ensure R is a rotation matrix
    [U S V] = svd(R{i});
    R{i} = U * eye(3) * V';
   
    P{i} = K * [R{i} t{i}];
    plot_camera(P{i}, 800, 600, 200);
end

% ToDo: in the report explain how the optical center is computed in the
%       provided code

[ny,nx] = size(T);
p1 = [0 0 0]';
p2 = [nx 0 0]';
p3 = [nx ny 0]';
p4 = [0 ny 0]';
% Draw planar pattern
vgg_scatter_plot([p1 p2 p3 p4 p1], 'g');
% Paint image texture
surface('XData',[0 nx; 0 nx],'YData',[0 0; 0 0],'ZData',[0 0; -ny -ny],'CData',T,'FaceColor','texturemap');
colormap(gray);
axis equal;

%% Plot a static camera with moving calibration pattern.
figure; hold;
plot_camera(K * eye(3,4), 800, 600, 200);
% ToDo: complete the call to the following function with the proper
%       coordinates of the image corners in the new reference system
for i = 1:N
    vgg_scatter_plot( [[R{i}*p1+t{i} R{i}*p2+t{i} R{i}*p3+t{i} R{i}*p4+t{i} R{i}*p1+t{i}]], 'r');
end

%% Augmented reality: Plot some 3D points on every camera.
[Th, Tw] = size(Tg);
cube = [0 0 0; 1 0 0; 1 0 0; 1 1 0; 1 1 0; 0 1 0; 0 1 0; 0 0 0; 0 0 1; 1 0 1; 1 0 1; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 0 1; 0 0 0; 1 0 0; 1 0 0; 1 0 1; 1 0 1; 0 0 1; 0 0 1; 0 0 0; 0 1 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 1 0; 0 0 0; 0 1 0; 0 1 0; 0 1 1; 0 1 1; 0 0 1; 0 0 1; 0 0 0; 1 0 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 1 0 1; 1 0 1; 1 0 0 ]';

X = (cube - .5) * Tw / 4 + repmat([Tw / 2; Th / 2; -Tw / 8], 1, length(cube));

for i = 1:N
    figure; colormap(gray);
    imagesc(Ig{i});
    hold on;
    x = euclid(P{i} * homog(X));
    vgg_scatter_plot(x, 'g');
end

