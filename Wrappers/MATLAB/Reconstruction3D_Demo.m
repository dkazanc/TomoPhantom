% GPLv3 license (ASTRA toolbox)
% Note that the TomoPhantom package is released under Apache License, Version 2.0

% Script to generate 3D analytical phantoms and their sinograms with added noise and artifacts
% Sinograms then reconstructed using ASTRA TOOLBOX 

% If one needs to modify/add phantoms, please edit Phantom2DLibrary.dat
% >>>> Prerequisites: ASTRA toolbox, if one needs to do reconstruction <<<<<
% adding paths
clear 
close all
fsep = '/';
addpath('compiled'); addpath('supplem'); 

% generate 3D phantom (modify your PATH bellow):
curDir   = pwd;
mainDir  = fileparts(curDir);
pathtoLibrary = sprintf([fsep '..' fsep 'PhantomLibrary' fsep 'models' fsep 'Phantom3DLibrary.dat'], 1i);
pathTP = strcat(mainDir, pathtoLibrary); % path to TomoPhantom parameters file

% generate a 3D phantom 
N = 256;
ModelNo = 13;
[G] = TomoP3DModel(ModelNo,N,pathTP);
figure; 
slice = round(0.5*N);
subplot(1,3,1); imagesc(G(:,:,slice), [0 1]); daspect([1 1 1]); colormap hot; title('Axial Slice');
subplot(1,3,2); imagesc(squeeze(G(:,slice,:)), [0 1]); daspect([1 1 1]); colormap hot; title('Y-Slice');
subplot(1,3,3); imagesc(squeeze(G(slice,:,:)), [0 1]); daspect([1 1 1]); colormap hot; title('X-Slice');

angles = linspace(0,179.99,360); % projection angles
DetectorWidth = 300; % detector column count (horizontal)
DetectorHeight = N; % detector row count (vertical)
%%
% using astra-toolbox to generate data
proj_geom = astra_create_proj_geom('parallel3d', 1, 1, DetectorHeight, DetectorWidth, angles*pi/180);
vol_geom = astra_create_vol_geom(N,N,N);

tic; [sinogram_id,proj3D_astra] = astra_create_sino3d_cuda(G, proj_geom, vol_geom); toc;
astra_mex_data3d('delete', sinogram_id);
%%

% figure; imshow(squeeze(sino3D_astra(:,1,:)),[]);
tic; proj3D_tomophant = TomoP3DModelSino(ModelNo, DetectorHeight, DetectorWidth, N, single(angles), pathTP); toc;
% figure; imshow(sino_tomophan3D, []);

% load astraprojs.mat
% comparing 2D analytical projections with ASTRA numerical projections
%compar_im = proj250_astra;
% sino_tomophan3D = reshape(sino_tomophan3D, [DetectorWidth,length(angles),DetectorHeight]);
slice2 = 160;
compar_im = squeeze(proj3D_astra(:,slice2,:));
sel_im = proj3D_tomophant(:,:,slice2);
disp(norm(sel_im(:) - compar_im(:))/norm(compar_im(:)))

% figure;imshow(squeeze(sino_tomophan3D(:,150,:)), []);
max_val = 100;
figure; 
subplot(1,3,1); imagesc(sel_im, [0 max_val]); title('Analytical projection');
subplot(1,3,2); imagesc(compar_im, [0 max_val]); title('Numerical projection');
subplot(1,3,3); imagesc(abs(sel_im - compar_im), [0 max_val]); title('image error');

figure; 
subplot(1,3,1); imagesc(squeeze(proj3D_astra(:,slice2,:)), [0 max_val]); title('Astra projection');
subplot(1,3,2); imagesc(squeeze(proj3D_astra(slice2,:,:)), [0 max_val]); title('Tangentogram');
subplot(1,3,3); imagesc(squeeze(proj3D_astra(:,:,slice2)), [0 max_val]); title('Sinogram');

figure; 
subplot(1,3,1); imagesc(squeeze(proj3D_tomophant(:,:,slice2)), [0 max_val]); title('Analytical projection');
subplot(1,3,2); imagesc(squeeze(proj3D_tomophant(slice2,:,:))', [0 max_val]); title('Tangentogram');
subplot(1,3,3); imagesc(squeeze(proj3D_tomophant(:,slice2,:)), [0 max_val]); title('Sinogram');


%%