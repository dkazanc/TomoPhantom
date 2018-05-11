% GPLv3 license (ASTRA toolbox)
% Note that the TomoPhantom package is released under Apache License, Version 2.0

% Script to generate 3D analytical object without using Phantom3DLibrary.dat


close all;clc;clear;
% adding paths
fsep = '/';
pathtoModels = sprintf(['..' fsep 'functions' fsep 'models' fsep], 1i);
addpath(pathtoModels);
addpath('compiled'); addpath('supplem'); 


% Define phantom dimensions
N = 256; % x-y-z size (cubic image)

% define parameters
paramsObject.Ob = 'gaussian';
paramsObject.C0 = 1; 
paramsObject.x0 = -0.25;
paramsObject.y0 = 0.1;
paramsObject.z0 = 0.0;
paramsObject.a = 0.2;
paramsObject.b = 0.35;
paramsObject.c = 0.7;
paramsObject.phi1 = 30;
paramsObject.phi2 = 60;
paramsObject.phi3 = -25;

% generate 3D phantom [N x N x N]:
tic; [G] = TomoP3DObject(paramsObject.Ob,paramsObject.C0, paramsObject.x0, paramsObject.y0, paramsObject.z0, paramsObject.a, paramsObject.b, paramsObject.c, paramsObject.phi1, paramsObject.phi2, paramsObject.phi3, N); toc;

% check 3 projections
figure; 
slice = round(0.5*N);
subplot(1,3,1); imagesc(G(:,:,slice), [0 1]); daspect([1 1 1]); colormap hot; title('Axial Slice');
subplot(1,3,2); imagesc(squeeze(G(:,slice,:)), [0 1]); daspect([1 1 1]); colormap hot; title('Y-Slice');
subplot(1,3,3); imagesc(squeeze(G(slice,:,:)), [0 1]); daspect([1 1 1]); colormap hot; title('X-Slice');