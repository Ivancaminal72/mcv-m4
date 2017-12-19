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

I2 = apply_H(I, H,'fit');
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
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

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
L1=lr1;
M1=lr2;

L2=lr3;
M2=lr4;

Eq=[L1(1)*M1(1), L1(1)*M1(2) + L1(2)*M1(1), L1(2)*M1(2); L2(1)*M2(1), L2(1)*M2(2) + L2(2)*M2(1), L2(2)*M2(2)];
s=null(Eq);
S=[s(1),s(2);s(2),s(3)];
K=chol(S);
K=inv(K);