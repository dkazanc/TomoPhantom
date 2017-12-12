%Copyright 2017 Daniil Kazantsev
%Licensed under the Apache License, Version 2.0 (the "License");

% Perform 2D deformation according to the tranform proposed in the paper:
% D. Kazantsev, V. Pickalov "New iterative reconstruction methods for fan-beam tomography", IPSE, 2017

close all;clc;clear;
% adding paths
addpath('../functions/models/'); addpath('install/'); 

ModelNo = 7; % Select a model (0 - 43  )
% Define phantom dimensions
N = 512; % x-y size (squared image)

% generate 2D phantom:
curDir   = pwd;
mainDir  = fileparts(curDir);
pathTP = strcat(mainDir,'/functions/models/Phantom2DLibrary.dat'); % path to TomoPhantom parameters file
[G] = buildPhantom2D(ModelNo,N,pathTP);
figure; imagesc(G, [0 1]); daspect([1 1 1]); colormap hot;


FocalP = 2.5; % focal point distance
RFP = 1./FocalP;
AngleTransform = 15; % angle in degrees

% deform forward
DeformType = 0; % 0 - forward deformation
G_deformed = DeformObject_C(double(G), RFP, AngleTransform, DeformType);
% deform back
DeformType = 1; % 1 - inverse 
G_inv = DeformObject_C(double(G_deformed), RFP, AngleTransform, DeformType);
figure; 
subplot(1,2,1); imagesc(G_deformed, [0 1]); title('Deformed Phantom'); daspect([1 1 1]); colormap hot
subplot(1,2,2); imagesc(G_inv, [0 1]); title('Inversely Deformed Phantom'); daspect([1 1 1]); colormap hot
MSE = norm(G_inv(:) - G(:))./norm(G(:));
fprintf('%s %f \n', 'Error (NMSE)', MSE);
