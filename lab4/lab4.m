%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 4: Reconstruction from two views (knowing internal camera parameters) 

addpath('../lab2/sift'); % ToDo: change 'sift' to the correct path where you have the sift functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Triangulation

% ToDo: create the function triangulate.m that performs a triangulation
%       with the homogeneous algebraic method (DLT)
%
%       The entries are (x1, x2, P1, P2, imsize), where:
%           - x1, and x2 are the Euclidean coordinates of two matching 
%             points in two different images.
%           - P1 and P2 are the two camera matrices
%           - imsize is a two-dimensional vector with the image size

%% Test the triangulate function
% Use this code to validate that the function triangulate works properly

P1 = eye(3,4);
c = cosd(15); s = sind(15);
R = [c -s 0; s c 0; 0 0 1];
t = [.3 0.1 0.2]';
P2 = [R t];
n = 8;
X_test = [rand(3,n); ones(1,n)] + [zeros(2,n); 3 * ones(1,n); zeros(1,n)];
x1_test = euclid(P1 * X_test);
x2_test = euclid(P2 * X_test);

N_test = size(x1_test,2);
X_trian = zeros(4,N_test);
for i = 1:N_test
    X_trian(:,i) = triangulate(x1_test(:,i), x2_test(:,i), P1, P2, [2 2]);
end

% error
euclid(X_test) - euclid(X_trian)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Reconstruction from two views

%% Read images
Irgb{1} = imread('Data/0001_s.png');
Irgb{2} = imread('Data/0002_s.png');
I{1} = sum(double(Irgb{1}), 3) / 3 / 255;
I{2} = sum(double(Irgb{2}), 3) / 3 / 255;
[h,w] = size(I{1});


%% Compute keypoints and matches.
points = cell(2,1);
descr = cell(2,1);
for i = 1:2
    [points{i}, descr{i}] = sift(I{i}, 'Threshold', 0.01);
    points{i} = points{i}(1:2,:);
end

matches = siftmatch(descr{1}, descr{2});

% Plot matches.
figure();
plotmatches(I{1}, I{2}, points{1}, points{2}, matches, 'Stacking', 'v');


%% Fit Fundamental matrix and remove outliers.
x1 = points{1}(:, matches(1, :));
x2 = points{2}(:, matches(2, :));
[F, inliers] = ransac_fundamental_matrix(homog(x1), homog(x2), 2.0);

% Plot inliers.
inlier_matches = matches(:, inliers);
figure;
plotmatches(I{1}, I{2}, points{1}, points{2}, inlier_matches, 'Stacking', 'v');

x1 = points{1}(:, inlier_matches(1, :));
x2 = points{2}(:, inlier_matches(2, :));

vgg_gui_F(Irgb{1}, Irgb{2}, F');


%% Compute candidate camera matrices.

% Camera calibration matrix
K = [2362.12 0 1520.69; 0 2366.12 1006.81; 0 0 1];
scale = 0.3;
H = [scale 0 0; 0 scale 0; 0 0 1];
K = H * K;


% ToDo: Compute the Essential matrix from the Fundamental matrix
E = K'*F*K;

% ToDo: write the camera projection matrix for the first camera
P1 = K*eye(3,4);
%P1(1:3,1:3) = K; %????

% ToDo: write the four possible matrices for the second camera
[U,D,V] = svd(E);

W = [0, -1, 0;
     1,  0, 0;
     0,  0, 1];
% HINT: You may get improper rotations; in that case you need to change
%       their sign.
% Let Rx be rotation matrices, you may check:
R1= U*W*V';
if det(R1) < 0
    R1 = -R1;
end
R2= U*W'*V';
if det(R2) < 0
    R2 = -R2;
end
% Let Tx be translation vectors
T1=U(:,3);
T2=-U(:,3);

Pc2 = {};
Pc2{1} = K*[R1,T1];
Pc2{2} = K*[R1,T2];
Pc2{3} = K*[R2,T1];
Pc2{4} = K*[R2,T2];

% plot the first camera and the four possible solutions for the second
figure;
plot_camera(P1,w,h);
plot_camera(Pc2{1},w,h);
plot_camera(Pc2{2},w,h);
plot_camera(Pc2{3},w,h);
plot_camera(Pc2{4},w,h);


%% Reconstruct structure
% ToDo: Choose a second camera candidate by triangulating a match.

for c=1:4
    P2 = Pc2{c};
    % Triangulate all matches.
    N = size(x1,2);
    X = zeros(4,N);
    for i = 1:N
        X(:,i) = triangulate(x1(:,i), x2(:,i), P1, P2, [w h]);
        Xc1 = P1*X;
        Xc2 = P2*X;
    end
    if(all(Xc1(3,:) > 0) && all(Xc2(3,:) > 0))
        break;
    end
end

%% Plot with colors
r = interp2(double(Irgb{1}(:,:,1)), x1(1,:), x1(2,:));
g = interp2(double(Irgb{1}(:,:,2)), x1(1,:), x1(2,:));
b = interp2(double(Irgb{1}(:,:,3)), x1(1,:), x1(2,:));
Xe = euclid(X);
figure; hold on;
plot_camera(P1,w,h);
plot_camera(P2,w,h);
for i = 1:length(Xe)
    scatter3(Xe(1,i), Xe(3,i), -Xe(2,i), 5^2, [r(i) g(i) b(i)]/255, 'filled');
end
axis equal;


%% Compute reprojection error.

% ToDo: compute the reprojection errors
%   Reprojection of the points (referenciated in the world coordinate system)
%   into the two camera planes.
    X1p = euclid(P1*X);
    X2p = euclid(P2*X);
    dsts1 = zeros(1,length(X1p));
    dsts2 = zeros(1,length(X1p));
    for i = 1:length(X1p)
        dsts1(i) = pdist([X1p(:,i)';x1(:,i)'], 'euclidean');
        dsts2(i) = pdist([X2p(:,i)';x2(:,i)'], 'euclidean');
    end
    
%   Plot the histogram of reprojection errors, and
%   plot the mean reprojection error
    mu1 = mean(dsts1);
    mu2 = mean(dsts2);
    hold on;
    histogram(dsts1)
    histogram(dsts2)
    line([mu1 mu1],ylim,'LineWidth', 2, 'Color','b');
    line([mu2 mu2],ylim,'LineWidth', 2, 'Color','r');
    hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Depth map computation with local methods (SSD)

% Data images: 'scene1.row3.col3.ppm','scene1.row3.col4.ppm'
% Disparity ground truth: 'truedisp.row3.col3.pgm'

% Write a function called 'stereo_computation' that computes the disparity
% between a pair of rectified images using a local method based on a matching cost 
% between two local windows.
% 
% The input parameters are 5:
% - left image
% - right image
% - minimum disparity
% - maximum disparity
% - window size (e.g. a value of 3 indicates a 3x3 window)
% - matching cost (the user may able to choose between SSD and NCC costs)
%
% In this part we ask to implement only the SSD cost
%
% Evaluate the results changing the window size (e.g. 3x3, 9x9, 20x20,
% 30x30) and the matching cost. Comment the results.
%
% Note 1: Use grayscale images
% Note 2: For this first set of images use 0 as minimum disparity 
% and 16 as the the maximum one.

% Read images
imL_rgb = imread('Data/scene1.row3.col3.ppm');
imR_rgb = imread('Data/scene1.row3.col4.ppm');
imL = sum(double(imL_rgb), 3) / 3 / 255;
imR = sum(double(imR_rgb), 3) / 3 / 255;


window_size = 9;
%Compute disparity between imL and imR
disparity = stereo_computation(imL, imR, 0, 16, window_size, 'SSD');

figure;
imshow(disparity);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Depth map computation with local methods (NCC)

% Complete the previous function by adding the implementation of the NCC
% cost.
%
% Evaluate the results changing the window size (e.g. 3x3, 9x9, 20x20,
% 30x30) and the matching cost. Comment the results.

window_size = 9;
%Compute disparity between imL and imR
disparity = stereo_computation(imL, imR, 0, 16, window_size, 'NCC');

figure;
imshow(disparity);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Depth map computation with local methods

% Data images: '0001_rectified_s.png','0002_rectified_s.png'
imL_rgb = imread('Data/0001_rectified_s.png');
imR_rgb = imread('Data/0002_rectified_s.png');
imL = sum(double(imL_rgb), 3) / 3 / 255;
imR = sum(double(imR_rgb), 3) / 3 / 255;

% Test the functions implemented in the previous section with the facade
% images. Try different matching costs and window sizes and comment the
% results.
% Notice that in this new data the minimum and maximum disparities may
% change.

window_size = 9;
%Compute disparity between imL and imR
disparity = stereo_computation(imL, imR, 0, 16, window_size, 'NCC');

figure;
imshow(disparity);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Bilateral weights

% Modify the 'stereo_computation' so that you can use bilateral weights (or
% adaptive support weights) in the matching cost of two windows.
% Reference paper: Yoon and Kweon, "Adaptive Support-Weight Approach for Correspondence Search", IEEE PAMI 2006
%
% Comment the results and compare them to the previous results (no weights).
%
% Note: Use grayscale images (the paper uses color images)

% Read images
imL_rgb = imread('Data/scene1.row3.col3.ppm');
imR_rgb = imread('Data/scene1.row3.col4.ppm');
imL = sum(double(imL_rgb), 3) / 3 / 255;
imR = sum(double(imR_rgb), 3) / 3 / 255;

window_size = 30;
%Compute disparity between imL and imR
disparity = stereo_computation(imL, imR, 0, 16, window_size, "BW");

figure;
imshow(disparity);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  Stereo computation with Belief Propagation

% Use the UGM library used in module 2 and implement a  
% stereo computation method that minimizes a simple stereo energy with 
% belief propagation. 
% For example, use an L2 or L1 pixel-based data term (SSD or SAD) and 
% the same regularization term you used in module 2. 
% Or pick a stereo paper (based on belief propagation) from the literature 
% and implement it. Pick a simple method or just simplify the method they propose.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  Depth computation with Plane Sweeping

% Implement the plane sweeping method explained in class.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  Depth map fusion 

% In this task you are asked to implement the depth map fusion method
% presented in the following paper:
% B. Curless and M. Levoy. A Volumetric Method for Building Complex
% Models from Range Images. In Proc. SIGGRAPH, 1996.
%
% 1. Use the set of facade images 00xx_s.png to compute depth maps 
% corresponding to different views (and optionally from different pairs of 
% images for the same view).
% 2. Then convert each depth map to a signed distance function defined in 
% a disretized volume (using voxels).
% 3. Average the different signed distance functions, the resulting 
% signed distance is called D.
% 4. Set as occupied voxels (those representing the surface) those 
% where D is very close to zero. The rest of voxels will be considered as 
% empty.
%
% For that you need to compute a depth map from a pair of views in general
% position (non rectified). Thus, you may either use the plane sweep
% algorithm (if you did it) or the local method for estimating depth
% (mandatory task) together with the following rectification method which 
% has an online demo available: 
% http://demo.ipol.im/demo/m_quasi_euclidean_epipolar_rectification/


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  New view synthesis

% In this task you are asked to implement part of the new view synthesis method
% presented in the following paper:
% S. Seitz, and C. Dyer, View morphing, Proc. ACM SIGGRAPH 1996.

% You will use a pair of rectified stereo images (no need for prewarping
% and postwarping stages) and their corresponding ground truth disparities
% (folder "new_view").
% Remember to take into account occlusions as explained in the lab session.
% Once done you can apply the code to the another pair of rectified images 
% provided in the material and use the estimated disparities with previous methods.
