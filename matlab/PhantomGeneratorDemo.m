% Script to generate analytical phantoms that can be used to test
% reconstruction algorithms.
% If one needs to modify/add phantoms just edit PhantomLibrary.dat
% ver 0.1, 1.07.17

close all;clc;clear all;
% adding paths
addpath('../models/'); addpath('supp/'); 

ModelNo = 43; % Select a model (0 - 42  )
% Define phantom dimension
N = 400; % x-y size (squared image)

% generate the 2D phantom:
[G] = buildPhantom(ModelNo,N);
figure(1); imshow(G, []);

%%
% perform deformation according to the tranform proposed in paper [1] (see readme)
FocalP = 4; % focal point distance
RFP = 1./FocalP;
DeformType = 0; % 0 - forward, 1 - inverse 
angle = 15;
angles_rad = (pi/180)*angle;

% deform forward
G_deformed = DeformObject(G, RFP, angles_rad, DeformType);
% deform back
DeformType = 1;
G_inv = DeformObject(G_deformed, RFP, angles_rad, DeformType);
figure; 
subplot(1,2,1); imshow(G_deformed, []); title('Deformed Phantom');
subplot(1,2,2); imshow(G_inv, []); title('Inversely Deformed Phantom');
%%