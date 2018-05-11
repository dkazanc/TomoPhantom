% GPLv3 license (ASTRA toolbox)
% Note that the TomoPhantom package is released under Apache License, Version 2.0

% Script to generate 2D analytical object and the corresponding sinogram 
% without using Phantom2DLibrary.dat

close all;clc;clear;
% adding paths
fsep = '/';
pathtoModels = sprintf(['..' fsep 'functions' fsep 'models' fsep], 1i);
addpath(pathtoModels);
addpath('compiled'); addpath('supplem'); 

% Define object dimensions
N = 256; % x-y-z size (cubic image)

% define parameters
paramsObject.Ob = 'gaussian';
paramsObject.C0 = 1; 
paramsObject.x0 = -0.25;
paramsObject.y0 = 0.1;
paramsObject.a = 0.2;
paramsObject.b = 0.35;
paramsObject.phi1 = 30;

% generate 2D phantom [N x N ]:
[G1] = TomoP2DObject(paramsObject.Ob,paramsObject.C0, paramsObject.x0, paramsObject.y0, paramsObject.a, paramsObject.b, paramsObject.phi1, N);
figure; imagesc(G1, [0 1]); daspect([1 1 1]); colormap hot;

% generate corresponding 2D sinogram 
angles = single(linspace(0,180,N)); % projection angles
P = round(sqrt(2)*N);

[F1] = TomoP2DObjectSino(paramsObject.Ob,paramsObject.C0, paramsObject.x0, paramsObject.y0, paramsObject.a, paramsObject.b, paramsObject.phi1, N, P, angles);
figure; imagesc(F1, [0 50]); colormap hot;

% generate another object 
paramsObject.Ob = 'rectangle';
paramsObject.C0 = 0.75; 
paramsObject.x0 = 0.25;
paramsObject.y0 = -0.1;
paramsObject.a = 0.2;
paramsObject.b = 0.4;
paramsObject.phi1 = -45;

[G2] = TomoP2DObject(paramsObject.Ob,paramsObject.C0, paramsObject.x0, paramsObject.y0, paramsObject.a, paramsObject.b, paramsObject.phi1, N);
figure; imagesc(G2, [0 1]); daspect([1 1 1]); colormap hot;
% generate corresponding 2D sinogram 
[F2] = TomoP2DObjectSino(paramsObject.Ob,paramsObject.C0, paramsObject.x0, paramsObject.y0, paramsObject.a, paramsObject.b, paramsObject.phi1, N, P, angles);
figure; imagesc(F2, [0 50]); colormap hot;

% build a new model
G = G1 + G2; 
F = F1 + F2; 

figure;
subplot(1,2,1); imagesc(G, [0 1]); daspect([1 1 1]); colormap hot; title('New model');
subplot(1,2,2); imagesc(F, [0 70]); daspect([1 1 1]); colormap hot; title('Sinogram');

% reconstruct the model
REC = rec2Dastra(F, double(angles*pi/180), P, N, 'cpu');
figure; imagesc(REC, [0 1]); daspect([1 1 1]); colormap hot;