% Script to generate 3D analytical phantoms
% If one needs to modify/add phantoms just edit Phantom3DLibrary.dat

close all;clc;clear all;
% adding paths
addpath('../functions/'); addpath('../functions/models/');


ModelNo = 04; % Select a model
% Define phantom dimensions
N = 256; % x-y-z size (cubic image)

% generate 3D phantom:
cd ../functions/
[G] = buildPhantom3D(ModelNo,N);
cd ../matlab/

% check the cenral slice
figure(1); imshow(G(:,:,round(0.5*N)), []);

% see the whole 3D Phantom
% figure(2);
% filename = strcat('ModelNo',num2str(ModelNo));
% counter = 1;
% for i = 1:N
%     imshow(G(:,:,i), [0 1]);
%     pause(0.01);
    
    % write tiff images
%     IM = im2uint16(G(:,:,i));
%     setStrNo = num2str(counter);
%     if (counter < 10)
%         filename_save = strcat(filename,'_','000',setStrNo, '.tiff');
%     elseif ((counter >= 10) &&  (counter < 100))
%         filename_save = strcat(filename,'_','00',setStrNo,'.tiff');
%     elseif ((counter >= 100) &&  (counter < 1000))
%         filename_save = strcat(filename,'_','0',setStrNo, '.tiff');
%     else
%         filename_save = strcat(filename,'_',setStrNo, '.tiff');
%     end
%     imwrite(IM,filename_save,'tiff');
%     counter = counter + 1;
% end
% close (figure(2));
%%
% obtaining 3D parallel-beam exact sinogram using TomoPhantom
angles = 1:0.5:180;
det = round(sqrt(2)*N);
tic; 
cd ../functions/
sino_tomophan3D = buildSino3D(ModelNo, N, det, single(angles)); 
cd ../matlab/
toc;
%%
% obtaining 3D parallel-beam sinogram of a phantom using ASTRA-toolbox
proj_geom = astra_create_proj_geom('parallel', 1, det, (angles*pi/180));
vol_geom = astra_create_vol_geom(N,N);
sino_astra3D = zeros(length(angles),det,N,'single');

tic;
for i = 1:N
[sinogram_id, sino_astra] = astra_create_sino_cuda(G(:,:,i), proj_geom, vol_geom);
sino_astra3D(:,:,i) = single(sino_astra);
astra_mex_data2d('delete', sinogram_id);
end
toc;

% calculate residiual norm (the error is expected since projection models not the same)
err_diff = norm(sino_tomophan3D(:) - sino_astra3D(:))./norm(sino_astra3D(:));
fprintf('%s %.4f\n', 'NRMSE for sino residuals:', err_diff);
figure; 
subplot(1,2,1); imshow(sino_tomophan3D(:,:,128), []); colormap hot; colorbar; title('Exact sinogram');
subplot(1,2,2); imshow(sino_astra3D(:,:,128), []); colormap hot;colorbar; title('Discrete sinogram');
%%
% Reconstruction using ASTRA-toolbox (FBP)
% Create a data object for the reconstruction
cfg = astra_struct('FBP_CUDA');
cfg.FilterType = 'Ram-Lak';

FBP3D = zeros(N,N,N,'single');

% choose which sinogram to substitute
for i = 1:N
rec_id = astra_mex_data2d('create', '-vol', vol_geom);
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sino_tomophan3D(:,:,i));
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;

alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', alg_id);
 
% Get the result
rec = astra_mex_data2d('get', rec_id);

% % Clean up. Note that GPU memory is tied up in the algorithm object,
% % and main RAM in the data objects.
astra_mex_algorithm('delete', alg_id);
astra_mex_data2d('delete', sinogram_id);
astra_mex_data2d('delete', rec_id);
FBP3D(:,:,i) = single(rec);
end
figure(4); imshow(FBP3D(:,:,round(0.5*N)), []);

%%
% angles = 1:0.5:180;
% sino = radon(G(:,:,128), angles);
% det = 367;
% recon_size = 256;
% 
% cd ../functions/
% F = buildSino3D(01, recon_size, det, single(angles));
% cd ../matlab/
% 
% figure(3); imshow(sino, []);
% figure(4); imshow(F(:,:,128)', []);
% 
% 
% FBP_analyt = iradon(F(:,:,128)',angles,recon_size);
% FBP_discr = iradon(sino,angles,recon_size);
% %%