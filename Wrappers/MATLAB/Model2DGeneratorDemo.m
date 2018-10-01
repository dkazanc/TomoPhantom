% GPLv3 license (ASTRA toolbox)
% Note that the TomoPhantom package is released under Apache License, Version 2.0

% Script to generate 2D analytical phantoms and sinograms (parallel beam)
% that can be used to test reconstruction algorithms without the "Inverse
% Crime"
% If one needs to modify/add phantoms, please edit Phantom2DLibrary.dat
% >>>> Prerequisites: ASTRA toolbox, if one needs to do reconstruction <<<<<

close all;clc;clear;
% adding paths
fsep = '/';
pathtoModels = sprintf(['..' fsep 'functions' fsep 'models' fsep], 1i);
addpath(pathtoModels);
addpath('compiled'); addpath('supplem'); 

ModelNo = 4; % Select a model from Phantom2DLibrary.dat
% Define phantom dimensions
N = 512; % x-y size (squared image)

% Generate 2D phantom:
curDir   = pwd;
mainDir  = fileparts(curDir);
pathtoLibrary = sprintf([fsep 'functions' fsep 'models' fsep 'Phantom2DLibrary.dat'], 1i);
pathTP = strcat(mainDir, pathtoLibrary); % path to TomoPhantom parameters file
[G] = TomoP2DModel(ModelNo,N,pathTP); 
figure; imagesc(G, [0 1]); daspect([1 1 1]); colormap hot;
%%
fprintf('%s \n', 'Generating sinogram analytically and numerically with Matlab (radon)...');
% generate angles
angles = linspace(0,180,N); % projection angles

% lets use Matlab's radon function
[F_d,xp] = radon(G,angles); % discrete sinogram
P = size(F_d,1); %detectors dimension
F_d = F_d';

% generate the 2D analytical parallel beam sinogram
[F_a] = TomoP2DModelSino(ModelNo, N, P, single(angles), pathTP, 'radon'); 

figure; 
subplot(1,2,1); imshow(F_a, []); title('Analytical Sinogram');
subplot(1,2,2); imshow(F_d, []); title('Numerical Sinogram');

% calculate residiual norm (the error is expected since projection models not the same)
err_diff = norm(F_a(:) - F_d(:))./norm(F_d(:));
fprintf('%s %.4f\n', 'NMSE for sino residuals:', err_diff);

% reconstructing with FBP (iradon)
FBP_F_a = iradon(F_a',angles,N);
FBP_F_d = iradon(F_d',angles,N);
figure; 
subplot(1,2,1); imagesc(FBP_F_a, [0 1]); title('Analytical Sinogram Reconstruction'); daspect([1 1 1]); colormap hot;
subplot(1,2,2); imagesc(FBP_F_d, [0 1]); title('Numerical Sinogram (radon) Reconstruction'); daspect([1 1 1]); colormap hot;
%%
fprintf('%s \n', 'Use the ASTRA-toolbox to generate numerical sinogram...');
% generate 2D analytical parallel beam sinogram (note the 'astra' opton)
[F_a] = TomoP2DModelSino(ModelNo, N, P, single(angles), pathTP, 'astra'); 
[F_num_astra] = sino2Dastra(G, (angles*pi/180), P, N, 'cpu');

% calculate residiual norm (the error is expected since projection models not the same)
err_diff = norm(F_a(:) - F_num_astra(:))./norm(F_num_astra(:));
fprintf('%s %.4f\n', 'NMSE for sino residuals:', err_diff);

figure; 
subplot(1,2,1); imshow(F_a, []); title('Analytical Sinogram');
subplot(1,2,2); imshow(F_num_astra, []); title('Numerical Sinogram (ASTRA)');
%%
fprintf('%s \n', 'Reconstruction using the ASTRA-toolbox (FBP)...');

rec_an = rec2Dastra(F_a, (angles*pi/180), P, N, 'cpu');
rec_num = rec2Dastra(F_num_astra, (angles*pi/180), P, N, 'cpu');

figure; 
subplot(1,2,1); imagesc(rec_an, [0 1]); daspect([1 1 1]); colormap hot; title('Analytical Sinogram Reconstruction [ASTRA]');
subplot(1,2,2); imagesc(rec_num, [0 1]); daspect([1 1 1]); colormap hot; title('Numerical Sinogram Reconstruction [ASTRA]');
%%
