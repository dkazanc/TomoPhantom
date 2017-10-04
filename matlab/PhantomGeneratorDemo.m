% Script to generate analytical phantoms and sinograms (parallel beam)
% that can be used to test reconstruction algorithms without Inverse Crime
% If one needs to modify/add phantoms just edit PhantomLibrary.dat

close all;clc;clear all;
% adding paths
addpath('../models/'); addpath('supp/'); 

ModelNo = 39; % Select a model (0 - 43  )
% Define phantom dimensions
N = 512; % x-y size (squared image)

% generate 2D phantom:
[G] = buildPhantom(ModelNo,N);
figure(1); imshow(G, []);

%%
fprintf('%s \n', 'Generating sinogram analytically and numerically...');
% generate angles
max_anlges = round(sqrt(2)*N);
angles = linspace(0,179.9,max_anlges); % projection angles

% lets use Matlab's radon function
tic;
[F_d,xp] = radon(G,angles); % discrete sinogram
toc;
P = size(F_d,1); %detectors dimension

% generate the 2D analytical parallel beam sinogram
tic;
[F_a] = buildSino(ModelNo,G,P,angles);
toc;
F_a = F_a.*(N/2); % scaling

figure(2); 
subplot(1,2,1); imshow(F_a, []); title('Analytical Sinogram');
subplot(1,2,2); imshow(F_d, []); title('Numerical Sinogram');


% reconstructing with FBP (iradon)
FBP_F_a = iradon(F_a,angles,N);
FBP_F_d = iradon(F_d,angles,N);
figure(3); 
subplot(1,2,1); imshow(FBP_F_a, []); title('Analytical Sinogram Reconstruction');
subplot(1,2,2); imshow(FBP_F_d, []); title('Numerical Sinogram Reconstruction');
%%
% run this once to compile
cd supp
mex DeformObject_C.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
cd ..
%%
% perform deformation according to the tranform proposed in the paper [1] (see readme file)
FocalP = 2.5; % focal point distance
RFP = 1./FocalP;
AngleTransform = 15; % angle in degrees

% deform forward
DeformType = 0; % 0 - forward deformation
G_deformed = DeformObject_C(G, RFP, AngleTransform, DeformType);
% deform back
DeformType = 1; % 1 - inverse 
G_inv = DeformObject_C(G_deformed, RFP, AngleTransform, DeformType);
figure; 
subplot(1,2,1); imshow(G_deformed, []); title('Deformed Phantom');
subplot(1,2,2); imshow(G_inv, []); title('Inversely Deformed Phantom');
MSE = norm(G_inv(:) - G(:))./norm(G(:));
fprintf('%s %f \n', 'Error (MSE)', MSE);
%%