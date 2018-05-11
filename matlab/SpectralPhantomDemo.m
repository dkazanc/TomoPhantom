% Note that the TomoPhantom package is released under Apache License, Version 2.0

% Script to generate material-specific phantom

% This phantom has been used in paper: 
% Joint image reconstruction method with correlative multi-channel prior
% for X-ray spectral computed tomography by D. Kazantsev et al.
% Inverse Problems 2018. 

close all;clc;clear;
% adding paths
fsep = '/';
pathtoModels = sprintf(['..' fsep 'functions' fsep 'models' fsep], 1i);
addpath(pathtoModels);
addpath('compiled'); addpath('supplem'); 


ModelNo = 11; 
% Define phantom dimension
N = 512; % x-y size (squared image)

% generate the 2D phantom:
curDir   = pwd;
mainDir  = fileparts(curDir);
pathtoLibrary = sprintf([fsep 'functions' fsep 'models' fsep 'Phantom2DLibrary.dat'], 1i);
pathTP = strcat(mainDir, pathtoLibrary); % path to TomoPhantom parameters file
[G] = TomoP2DModel(ModelNo,N,pathTP);

% create 4 phantoms with dedicated materials
G(G >= 1.29) = 0.8;
G1 = zeros(N,N);
G2 = zeros(N,N);
G3 = zeros(N,N);
G4 = zeros(N,N);

G1((G > 0.2) & (G < 0.4)) = 1;
G2((G > 0.7) & (G < 0.9)) = 2;
G3((G > 1.0) & (G < 1.15)) = 3;
G4(G > 1.17) = 4;
G1(G > 1) = 0;

figure;
subplot(2,2,1)       % add first plot in 2 x 2 grid
imagesc(G1, [0 1]); daspect([1 1 1]); colormap hot;
title('Material 1')

subplot(2,2,2)       % add second plot in 2 x 2 grid
imagesc(G2, [0 2]); daspect([1 1 1]); colormap hot;
title('Material 2')

subplot(2,2,3)       % add third plot in 2 x 2 grid
imagesc(G3, [0 3]); daspect([1 1 1]); colormap hot;
title('Material 3')

subplot(2,2,4)       % add fourth plot in 2 x 2 grid
imagesc(G4, [0 4]); daspect([1 1 1]); colormap hot;
title('Material 4')
