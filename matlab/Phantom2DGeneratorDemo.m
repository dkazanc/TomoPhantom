% Script to generate 2D analytical phantoms and sinograms (parallel beam)
% that can be used to test reconstruction algorithms without Inverse Crime
% If one needs to modify/add phantoms just edit Phantom2DLibrary.dat

close all;clc;clear;
% adding paths
addpath('../functions/models/'); addpath('supp/'); 

ModelNo = 06; % Select a model from Phantom2DLibrary.dat
% Define phantom dimensions
N = 512; % x-y size (squared image)

% generate 2D phantom:
pathTP = '/home/algol/Documents/MATLAB/TomoPhantom/functions/models/Phantom2DLibrary.dat'; % path to TomoPhantom parameters file
[G] = buildPhantom2D(ModelNo,N,pathTP);
figure; imagesc(G, [0 1]); daspect([1 1 1]); colormap hot;
%%
fprintf('%s \n', 'Generating sinogram analytically and numerically...');
% generate angles
angles = linspace(0,180,N); % projection angles

% lets use Matlab's radon function
[F_d,xp] = radon(G,angles); % discrete sinogram
P = size(F_d,1); %detectors dimension

% generate the 2D analytical parallel beam sinogram
[F_a] = buildSino2D(ModelNo, N, P, single(angles), pathTP, 'radon'); 

figure; 
subplot(1,2,1); imshow(F_a, []); title('Analytical Sinogram');
subplot(1,2,2); imshow(F_d, []); title('Numerical Sinogram');

% calculate residiual norm (the error is expected since projection models not the same)
err_diff = norm(F_a(:) - F_d(:))./norm(F_d(:));
fprintf('%s %.4f\n', 'NMSE for sino residuals:', err_diff);

% reconstructing with FBP (iradon)
FBP_F_a = iradon(F_a,angles,N);
FBP_F_d = iradon(F_d,angles,N);
figure; 
subplot(1,2,1); imagesc(FBP_F_a, [0 1]); title('Analytical Sinogram Reconstruction'); daspect([1 1 1]); colormap hot;
subplot(1,2,2); imagesc(FBP_F_d, [0 1]); title('Numerical Sinogram Reconstruction'); daspect([1 1 1]); colormap hot;
%%
fprintf('%s \n', 'Use ASTRA-toolbox to generate numerical sinogram...');
% >>>> Requirements: ASTRA toolbox, if one needs to do reconstruction <<<<<
% generate the 2D analytical parallel beam sinogram (note 'astra' option)
[F_a] = buildSino2D(ModelNo, N, P, single(angles), pathTP, 'astra'); 

proj_geom = astra_create_proj_geom('parallel', 1, P, (angles*pi/180));
vol_geom = astra_create_vol_geom(N,N);

[sinogram_id, sino_astra] = astra_create_sino_cuda(G, proj_geom, vol_geom);
astra_mex_data2d('delete', sinogram_id);

sinT = sino_astra';
% calculate residiual norm (the error is expected since projection models not the same)
err_diff = norm(F_a(:) - sinT(:))./norm(sinT(:));
fprintf('%s %.4f\n', 'NMSE for sino residuals:', err_diff);

figure; 
subplot(1,2,1); imshow(F_a, []); title('Analytical Sinogram');
subplot(1,2,2); imshow(sino_astra', []); title('Numerical Sinogram');
%%