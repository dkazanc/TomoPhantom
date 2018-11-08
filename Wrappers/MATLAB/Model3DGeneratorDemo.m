% GPLv3 license (ASTRA toolbox)
% Note that the TomoPhantom package is released under Apache License, Version 2.0

% Script to generate 3D analytical phantoms
% If one needs to modify/add phantoms, please edit Phantom3DLibrary.dat
% >>>> Prerequisites: ASTRA toolbox, if one needs to do reconstruction <<<<<

close all;clc;clear;
% adding paths
fsep = '/';
addpath('compiled'); addpath('supplem'); 

ModelNo = 13; % Select a model
% Define phantom dimensions
N = 256; % x-y-z size (cubic image)

% generate 3D phantom (modify your PATH bellow):
curDir   = pwd;
mainDir  = fileparts(curDir);
pathtoLibrary = sprintf([fsep '..' fsep 'PhantomLibrary' fsep 'models' fsep 'Phantom3DLibrary.dat'], 1i);
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
angles_num = round(0.5*pi*N); % angles number
angles = linspace(0,179.99,angles_num); % projection angles
Horiz_det = round(sqrt(2)*N); % detector column count (horizontal)
Vert_det = N; % detector row count (vertical) (no reason for it to be > N, so fixed)

disp('Using TomoPhantom to generate 3D projection data');
proj3D_tomophant = TomoP3DModelSino(ModelNo, Vert_det, Horiz_det, N, single(angles), pathTP);

figure; 
subplot(1,3,1); imagesc(squeeze(proj3D_tomophant(:,:,slice2)), [0 max_val]); title('Analytical projection');
subplot(1,3,2); imagesc(squeeze(proj3D_tomophant(slice2,:,:))', [0 max_val]); title('Tangentogram');
subplot(1,3,3); imagesc(squeeze(proj3D_tomophant(:,slice2,:)), [0 max_val]); title('Sinogram');

%%
