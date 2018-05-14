% GPLv3 license (ASTRA toolbox)
% Note that the TomoPhantom package is released under Apache License, Version 2.0

% Script to generate 2D analytical phantoms and their sinograms with added noise and artifacts
% Sinograms then reconstructed using ASTRA TOOLBOX 

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
fprintf('%s \n', 'Generating sinogram analytically and numerically using ASTRA-toolbox...');
angles = linspace(0,180,N); % projection angles
P = round(sqrt(2)*N); %detectors dimension
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
fprintf('%s \n', 'Adding noise and artifacts to analytical sinogram...');
dose =  1e4; % photon flux (controls noise level)
[sino_noise] = add_noise(F_a, dose, 'Poisson'); % adding Poisson noise
[sino_noise_zingers] = add_zingers(sino_noise, 0.5, 10); % adding zingers
[sino_noise_zingers_stripes] = add_stripes(sino_noise_zingers, 1, 1); % adding stripes

figure; imshow(sino_noise_zingers_stripes, []); title('Analytical Sinogram degraded with noise and artifacts');
%%
fprintf('%s \n', 'Reconstruction using ASTRA-toolbox (FBP)...');

rec_an = rec2Dastra(sino_noise_zingers_stripes, (angles*pi/180), P, N, 'cpu');
rec_num = rec2Dastra(F_num_astra, (angles*pi/180), P, N, 'cpu');

figure; 
subplot(1,2,1); imagesc(rec_an, [0 1]); daspect([1 1 1]); colormap hot; title('Degraded analytical sinogram reconstruction');
subplot(1,2,2); imagesc(rec_num, [0 1]); daspect([1 1 1]); colormap hot; title('Numerical sinogram reconstruction');
%%