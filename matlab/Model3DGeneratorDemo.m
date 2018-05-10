% GPLv3 license (ASTRA toolbox)
% Note that the TomoPhantom package is released under Apache License, Version 2.0

% Script to generate 3D analytical phantoms
% If one needs to modify/add phantoms, please edit Phantom3DLibrary.dat
% >>>> Prerequisites: ASTRA toolbox, if one needs to do reconstruction <<<<<

close all;clc;clear;
% adding paths
fsep = '/';
pathtoModels = sprintf(['..' fsep 'functions' fsep 'models' fsep], 1i);
addpath(pathtoModels);
addpath('compiled'); addpath('supplem'); 

ModelNo = 10; % Select a model
% Define phantom dimensions
N = 256; % x-y-z size (cubic image)

% generate 3D phantom (modify your PATH bellow):
curDir   = pwd;
mainDir  = fileparts(curDir);
pathtoLibrary = sprintf([fsep 'functions' fsep 'models' fsep 'Phantom3DLibrary.dat'], 1i);
pathTP = strcat(mainDir, pathtoLibrary); % path to TomoPhantom parameters file
tic; [G] = TomoP3DModel(ModelNo,N,pathTP); toc;

% check 3 projections
figure; 
slice = round(0.5*N);
subplot(1,3,1); imagesc(G(:,:,slice), [0 1]); daspect([1 1 1]); colormap hot; title('Axial Slice');
subplot(1,3,2); imagesc(squeeze(G(:,slice,:)), [0 1]); daspect([1 1 1]); colormap hot; title('Y-Slice');
subplot(1,3,3); imagesc(squeeze(G(slice,:,:)), [0 1]); daspect([1 1 1]); colormap hot; title('X-Slice');

% visualise/save the whole 3D Phantom
% % figure(2);
% filename = strcat('ModelNo',num2str(ModelNo));
% counter = 1;
% for i = 1:N
%     imshow(G(:,:,i), [0 1]);
%     pause(0.01);
% %     
% % % %     write png images
%     IM = im2uint16(G(:,:,i));    
%     setStrNo = num2str(counter);
%     if (counter < 10)
%         filename_save = strcat(filename,'_','000',setStrNo, '.png');
%     elseif ((counter >= 10) &&  (counter < 100))
%         filename_save = strcat(filename,'_','00',setStrNo,'.png');
%     elseif ((counter >= 100) &&  (counter < 1000))
%         filename_save = strcat(filename,'_','0',setStrNo, '.png');
%     else
%         filename_save = strcat(filename,'_',setStrNo, '.png');
%     end
%     imwrite(IM,filename_save,'png');
%     counter = counter + 1;
% end
% close (figure(2));
%%
fprintf('%s \n', 'Calculating 3D parallel-beam sinogram of a phantom using ASTRA-toolbox...');
angles = linspace(0,pi,N); % projection angles
det = round(sqrt(2)*N);
sino_astra3D = zeros(length(angles),det,N,'single');

tic;
for i = 1:N
sino_astra3D(:,:,i) = sino2Dastra(G(:,:,i), angles, det, N, 'gpu');
end
toc;

% calculate residiual norm (the error is expected since projection models not the same)
% err_diff = norm(sino_tomophan3D(:) - sino_astra3D(:))./norm(sino_astra3D(:));
% fprintf('%s %.4f\n', 'NRMSE for sino residuals:', err_diff);
% figure; 
% subplot(1,2,1); imagesc(sino_tomophan3D(:,:,slice)', [0 70]); colormap hot; colorbar; daspect([1 1 1]); title('Exact sinogram');
% subplot(1,2,2); imagesc(sino_astra3D(:,:,slice)', [0 70]); colormap hot; colorbar; daspect([1 1 1]); title('Discrete sinogram');
%%
fprintf('%s \n', 'Reconstruction using ASTRA-toolbox (FBP)...');
FBP3D = zeros(N,N,N,'single');

% choose which sinogram to substitute
for i = 1:N
FBP3D(:,:,i) = rec2Dastra(sino_astra3D(:,:,i), angles, det, N);
end
figure(3); 
subplot(1,3,1); imagesc(FBP3D(:,:,slice), [0 1]); daspect([1 1 1]); colormap hot; title('Axial Slice');
subplot(1,3,2); imagesc(squeeze(FBP3D(:,slice,:)), [0 1]); daspect([1 1 1]); colormap hot; title('Y-Slice');
subplot(1,3,3); imagesc(squeeze(FBP3D(slice,:,:)), [0 1]); daspect([1 1 1]); colormap hot; title('X-Slice');
%%
