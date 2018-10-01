% GPLv3 license (ASTRA toolbox)
% Note that the TomoPhantom package is released under Apache License, Version 2.0

% Script to generate 2D analytical temporal (2D + time) phantoms and sinograms (parallel beam)
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

ModelNo = 100; % Select a model from Phantom2DLibrary.dat
% Define phantom dimensions
N = 512; % x-y size (squared image)

% Generate 2D+t phantom:
curDir   = pwd;
mainDir  = fileparts(curDir);
pathtoLibrary = sprintf([fsep 'functions' fsep 'models' fsep 'Phantom2DLibrary.dat'], 1i);
pathTP = strcat(mainDir, pathtoLibrary); % path to TomoPhantom parameters file
[G] = TomoP2DModel(ModelNo,N,pathTP);
figure(1); imagesc(G, [0 1]); daspect([1 1 1]); title('2D+t model, t=3 here'); colormap hot;
%%
% Lets look at more finely discretized temporal model and related sinograms
ModelNo = 101; % Select a model from Phantom2DLibrary.dat
% Define phantom dimensions
N = 512; % x-y size (squared image)

% Generate 2D phantom:
curDir   = pwd;
mainDir  = fileparts(curDir);
pathtoLibrary = sprintf([fsep 'functions' fsep 'models' fsep 'Phantom2DLibrary.dat'], 1i);
pathTP = strcat(mainDir, pathtoLibrary); % path to TomoPhantom parameters file
[G] = TomoP2DModel(ModelNo,N,pathTP);

angles = linspace(0,180,N); % projection angles
P = round(sqrt(2)*N);
F = TomoP2DModelSino(ModelNo, N, P, single(angles), pathTP, 'radon');

figure(2);
for i = 1:350
    subplot(1,2,1); imagesc(G(:,:,i), [0 1]); daspect([1 1 1]); title('2D+t phantom'); colormap hot;
    subplot(1,2,2); imagesc(F(:,:,i), [0 165]); daspect([1 1 1]); title('Corresponding sinogram'); colormap hot;
    pause(0.01);
end
%%
% another temporal phantom with stationary and dynamic features
ModelNo = 102; % Select a model from Phantom2DLibrary.dat
% Define phantom dimensions
N = 512; % x-y size (squared image)

% Generate 2D+t phantom:
timeFrames = 25; %must be the same as in model
curDir   = pwd;
mainDir  = fileparts(curDir);
pathtoLibrary = sprintf([fsep 'functions' fsep 'models' fsep 'Phantom2DLibrary.dat'], 1i);
pathTP = strcat(mainDir, pathtoLibrary); % path to TomoPhantom parameters file
[G] = TomoP2DModel(ModelNo,N,pathTP);

angles = linspace(0,180,N); % projection angles
P = round(sqrt(2)*N);
F = TomoP2DModelSino(ModelNo, N, P, single(angles), pathTP, 'radon');

% reconstruct
FBP_F_a = zeros(N,N,timeFrames);
for i = 1:timeFrames
    FBP_F_a(:,:,i) = iradon(F(:,:,i)',angles,N);
end

fig_num = 2;
for i = 1:timeFrames
    figure(fig_num);
    subplot(1,3,1); imagesc(G(:,:,i), [0 1]); daspect([1 1 1]); title('2D + time phantom'); colormap hot;
    subplot(1,3,2); imagesc(F(:,:,i), [0 175]); daspect([1 1 1]); title('Corresponding sinogram'); colormap hot;
    subplot(1,3,3); imagesc(FBP_F_a(:,:,i), [0 1]); daspect([1 1 1]); title('FBP reconstruction'); colormap hot;
    pause(0.15);
    
    % if one needs animated gif here: 
%     filename = 'animat.gif';
%     del = 0.1; % time between animation frames
%     drawnow
%     frame = getframe(fig_num);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if (i == 1)
%         imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',del);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
%     end    
end
%%
