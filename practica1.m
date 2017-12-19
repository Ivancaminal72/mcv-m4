%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Practica 1: Transformacions d'imatges


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying transformations

%% Question 1
c1 = [1 1 1]';
H = [1 0 2; 
     0 1 1;
     0 0 1];
Hc1 = H*c1; Hc1=Hc1/Hc1(3)
plot(c1(1), c1(2), 'r+', Hc1(1), Hc1(2), 'b+');


%% 2. Similarities
I=imread('mondrian.jpg');
figure()
imshow(I)


%% Question 2
alpha2=30;
H = [cos(alpha2) -sin(alpha2) 0;
    sin(alpha2) cos(alpha2) 0;
    0 0 1];
I2 = apply_H(I, H, 'fit');
figure; imagesc(I); figure; imagesc(uint8(I2));
%imagesc(...) is the same as IMAGE(...) except the data 
%is scaled to use the full colormap.

%% Question 3
alpha3=-70;
H = [cos(alpha3) -sin(alpha3) 10;
    sin(alpha3) cos(alpha3) 350;
    0 0 1];
I2 = apply_H(I, H, 'fit');
figure; imagesc(I); figure; imagesc(uint8(I2));


%% Question 4
s=1.5;
alpha4=10;
H = [s*cos(alpha4) -s*sin(alpha4) 0;
    s*sin(alpha4) s*cos(alpha4) 0;
    0 0 1];
I2 = apply_H(I, H, 'fit');
figure; imagesc(I); figure; imagesc(uint8(I2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Affinities

%% Question 5
H=[2 -1 1;
    1 1 3;
    0 0 1];
I25 = apply_H(I, H, 'fit');
figure; imagesc(I); figure; imagesc(uint8(I25));


%% Question 6
I=imread('mondrian.jpg');
A=[2 -1;
   1 1];
[U D V]=svd(A);

R1=U*V'*inv(V');    %first rotation
R2=V';      %second rotation
S=D;        %scalinng
T=[1,3];    %translation

HR1=[R1(1,1) R1(1,2) 0;
    R1(2,1) R1(2,2) 0;
    0 0 1]; % first rotation matrix 3x3

HR2=[R2(1,1) R2(1,2) 0;
    R2(2,1) R2(2,2) 0;
    0 0 1]; %second rotatin matrix 3x3
HS=[S(1,1) S(1,2) 0;
    S(2,1) S(2,2) 0;
    0 0 1]; % scalinng matrix
HT=[1 0 T(1);
    0 1 T(2);
    0 0 1]; % translation matrix

% applying affinities separately
I2 = apply_H(I, HR1, 'fit');
I3 = apply_H(I2, HS, 'fit');
I4 = apply_H(I3, HR2, 'fit');
I5 = apply_H(I4, HT, 'fit');

%Hreconstruct
TT=[0 0 T(1);
    0 0 T(2);
    0 0 1];
Hrect=(HR1*HR2)*HS+TT
Hrect=Hrect-[0 0 0;0 0 0; 0 0 1];
Irecns=apply_H(I, Hrect,'fit');
figure()
imshow(uint8(Irecns))
title('image with H reconstructed')

figure()
subplot(3,2,1)
imshow(I)
title('original image')
subplot(3,2,2)
imshow(uint8(I25))
title('affinity of question 5')
subplot(3,2,3)
imshow(uint8(I2))
title('first rotation')
subplot(3,2,4)
imshow(uint8(I3))
title('scalinng')
subplot(3,2,5)
imshow(uint8(I4))
title('second rotation')
subplot(3,2,6)
imshow(uint8(I5))
title('translation matrix')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Homographies

%% Question 7, 8, 9
H = [1 0 0; 0 1 0; 0.001 0.001 1];
% I2 = apply_H(I, H, 'fit');
% 8- Show that an affinity maps points at infinity to points at infinity
l=[H(3,1) H(3,2) H(3,3)];
Ha=[1 0 0; 0 1 0; 0 0 1];

Hafi=Ha*H;

linf=(inv(Hafi))'*l'

% 9- show that a projectivity maps points at infinity to finite points
% in this part we have to show that the homography applyed to de infinity
% line points at a finite points

Linf=[0 0 1]';
L=(Hafi)'*Linf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Image rectification
I = imread('Data/0000_s.png');
imshow(I);


%% Question 10
A = load('Data/0000_s_info_lines.txt');
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
figure; imshow(I); hold;
plot(p1(1), p1(2), 'y+');
plot(p2(1), p2(2), 'y+');
plot(p3(1), p3(2), 'y+');
plot(p4(1), p4(2), 'y+');


%% Question 11
% q1 = ...
% q2 = ...
% q3 = ...
% q4 = ...


%% Question 12
% X = ...
% Y = ...
H = homography2d(X,Y);
I2 = apply_H(I, H, 'fit');
I3 = apply_H(I, H, 'original');

%% Show results
figure;imshow(uint8(I))
hold on
for i=1:4,
  plot(X(1,i),X(2,i),'y*');
end
figure;imshow(uint8(I3))
hold on
for i=1:4,
  plot(Y(1,i),Y(2,i),'y*');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Affine Rectification

%% Question 13
% l1 = cross([74 265 1],[128 278 1]);
% l2 = cross([67 342 1],[125 347 1]);
% l3 = cross([74 265 1],[67 342 1]);
% l4 = cross([128 278 1],[125 347 1]);
l1=cross(p1,p2);
l2=cross(p3,p4);
l3=cross(p5,p6);
l4=cross(p7,p8);
I = imread('Data/0000_s.png');
figure()
imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');


%% Question 14
v1 = cross(l1,l2);
v2 = cross(l3,l4);


%% Question 15
linf = cross(v1,v2);
linf = linf / linf(3);  % Prevent large values in linf


%% Question 16, 17
H =[1 0 0;
    0 1 0;
    linf(1) linf(2) linf(3)];
I2 = apply_H(I, H);

figure()
imshow(I2)





