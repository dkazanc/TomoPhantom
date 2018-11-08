% GPLv3 license (ASTRA toolbox)
% Note that the TomoPhantom package is released under Apache License, Version 2.0

% Script to generate 3D analytical phantoms and projection data then reconstructed using ASTRA TOOLBOX 

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

disp('Using TomoPhantom to generate 3D phantom');
N = 256;
ModelNo = 13;
[G] = TomoP3DModel(ModelNo,N,pathTP);
figure; 
slice = round(0.5*N);
subplot(1,3,1); imagesc(G(:,:,slice), [0 1]); daspect([1 1 1]); colormap hot; title('Axial Slice');
subplot(1,3,2); imagesc(squeeze(G(:,slice,:)), [0 1]); daspect([1 1 1]); colormap hot; title('Y-Slice');
subplot(1,3,3); imagesc(squeeze(G(slice,:,:)), [0 1]); daspect([1 1 1]); colormap hot; title('X-Slice');

angles_num = round(0.5*pi*N); % angles number
angles = linspace(0,179.99,angles_num); % projection angles
Horiz_det = round(sqrt(2)*N); % detector column count (horizontal)
Vert_det = N; % detector row count (vertical) (no reason for it to be > N, so fixed)
%%
disp('Using astra-toolbox (GPU) to generate 3D projection data');
proj3D_astra = sino3Dastra(G, angles, Vert_det, Horiz_det);
%%
disp('Using TomoPhantom to generate 3D projection data');
proj3D_tomophant = TomoP3DModelSino(ModelNo, Vert_det, Horiz_det, N, single(angles), pathTP);
%%

% comparing 2D analytical projections with ASTRA numerical projections
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
disp('Reconstructing data using ASTRA toolbox (CGLS method)');
reconstructon = rec3Dastra(proj3D_tomophant, angles, Vert_det, Horiz_det);
figure; 
slice = round(0.5*N);
subplot(1,3,1); imagesc(reconstructon(:,:,slice), [0 1]); daspect([1 1 1]); colormap hot; title('Axial Slice (reconstruction)');
subplot(1,3,2); imagesc(squeeze(reconstructon(:,slice,:)), [0 1]); daspect([1 1 1]); colormap hot; title('Y-Slice');
subplot(1,3,3); imagesc(squeeze(reconstructon(slice,:,:)), [0 1]); daspect([1 1 1]); colormap hot; title('X-Slice');

%%