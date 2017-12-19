%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 1: Image rectification


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying image transformations

% ToDo: create the function  "apply_H" that gets as input a homography and
% an image and returns the image transformed by the homography.
% The size of the transformed image has to be automatically set so as to 
% contain the whole transformed image.
% At some point you will need to interpolate the image values at some points,
% you may use the Matlab function "interp2" for that.


%% 1.1. Similarities
I=imread('Data/0005_s.png'); % we have to be in the proper folder

%---------------------------------------------------------------------
% ToDo: generate a matrix H which produces a similarity transformation
theta=1;
s=1;
t=[30 30];
H = [s*cosd(theta) -s*sind(theta) t(1);
    s*sind(theta) s*cosd(theta) t(2);
    0 0 1];
%---------------------------------------------------------------------

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));


%% 1.2. Affinities
% ToDo: generate a matrix H which produces an affine transformation
theta=90;
fi=50;
s=[2, 1];
t=[30, 30];

R1 = [cosd(theta), -sind(theta), 0;
      sind(theta),  cosd(theta), 0;
      0          ,  0          , 1];
  
R2a = [cosd(fi), -sind(fi), 0;
       sind(fi),  cosd(fi), 0;
       0       ,  0       , 1];
  
R2b = [cosd(-fi), -sind(-fi), 0;
       sind(-fi),  cosd(-fi), 0;
       0        ,  0        , 1];
  
S = [s(1), 0   , 0;
     0   , s(2), 0;
     0   , 0   , 1];

T = [0, 0, t(1);
     0, 0, t(2);
     0, 0,  0  ];

H = (R1*R2a*S*R2b)+T;
I2 = apply_H(I, H);
hold on;
subplot(1, 3, 1); imshow(I); 
subplot(1, 3, 2); imshow(uint8(I2));

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation
[U,D,V] = svd(H(1:2,1:2));
R1_d = U*V';
R2_d = V';
S_d = D;
T_d = H(1:2,3);

% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above
A = R1_d*V*S_d*R2_d;
H2 = zeros(3,3);
H2(3,3) = 1;
H2(1:2,3) = T_d;
H2(1:2,1:2) = A;
if(~any(any(round(H-H2))))
    disp('same matrices')
end

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before
I2 = apply_H(I, H2);
subplot(1, 3, 3); imshow(uint8(I2));
hold off;

%% 1.3 Projective transformations (homographies)

% ToDo: generate a matrix H which produces a projective transformation

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification


% choose the image points
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
i = 424;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 240;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 712;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 565;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';
% imshow(I);
% hold on;
% plot(p8(1), p8(2), 'r*', 'LineWidth', 2, 'MarkerSize', 15);


% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points

l1=cross(p1,p2);
l2=cross(p3,p4);
l3=cross(p5,p6);
l4=cross(p7,p8);

% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% ToDo: compute the homography that affinely rectifies the image
% Compute the vanishing points 
v1 = cross(l1,l2);
v2 = cross(l3,l4);

% Compute the line at infinity
linf = cross(v1,v2);
linf = linf / linf(3);

H =[1 0 0;
    0 1 0;
    linf(1) linf(2) linf(3)];

I2 = apply_H(I, H);
figure; imshow(uint8(I2));

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4

lr1=inv(H')*l1;
lr2=inv(H')*l2;
lr3=inv(H')*l3;
lr4=inv(H')*l4;

% show the transformed lines in the transformed image
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'r');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'r');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.

% 2 pairs of ortogonal lines
L1=lr1; M1=lr3;
L2=lr2; M2=lr4;

% system of equations
Eq=[L1(1)*M1(1), L1(1)*M1(2) + L1(2)*M1(1), L1(2)*M1(2); L2(1)*M2(1), L2(1)*M2(2) + L2(2)*M2(1), L2(2)*M2(2)];
%computing the null vector and S matrix
s=null(Eq);
S=[s(1),s(2);s(2),s(3)];
%Cholesky decomposition to find K matrix
K=chol(S);
%Compute the inverse of K to fing de homography Hs<-a
K=inv(K);
H=[K(1,1) K(1,2) 0; K(2,1) K(2,2) 0; 0 0 1];

%Computing the homography to the lines
L1=inv(H')*L1;
M1=inv(H')*M1;
L2=inv(H')*L2;
M2=inv(H')*M2;

I3 = apply_H(I2, H);
figure; imshow(uint8(I3));
hold on;
plot(t, -(L1(1)*t + L1(3)) / L1(2), 'y');
plot(t, -(M1(1)*t + M1(3)) / M1(2), 'r');
plot(t, -(L2(1)*t + L2(3)) / L2(2), 'y');
plot(t, -(M2(1)*t + M2(3)) / M2(2), 'r');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that use in every step.

I1=imread('Data/0001_s.png');
I1=I1(1:end,1:474,:);
A = load('Data/0001_s_info_lines.txt');

% indices of lines
i = 614;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 541;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 159;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 645;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';
%lines
l1=cross(p1,p2);
l2=cross(p3,p4);
l3=cross(p5,p6);
l4=cross(p7,p8);

% show the chosen lines in the image
figure;imshow(I1);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)



