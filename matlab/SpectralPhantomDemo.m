% generating material-specific spectral phantom

close all;clc;clear;
% adding paths
addpath('../functions/models/'); addpath('supp/'); 

ModelNo = 11; 
% Define phantom dimension
N = 512; % x-y size (squared image)

% generate the 2D phantom:
pathTP = '/home/algol/Documents/MATLAB/TomoPhantom/functions/models/Phantom2DLibrary.dat'; % path to TomoPhantom parameters file
[G] = buildPhantom2D(ModelNo,N,pathTP);
figure; imagesc(G, [0 1]); daspect([1 1 1]); colormap hot;

% create 4 phantoms with dedicated materials
G(G >= 1.29) = 0.8;
G1 = zeros(N,N);
G2 = zeros(N,N);
G3 = zeros(N,N);
G4 = zeros(N,N);

G1((G > 0.2) & (G < 0.4)) = 1;
G2((G > 0.7) & (G < 0.9)) = 2;
G3((G > 1.0) & (G < 1.15)) = 3;
G3(G > 1.17) = 4;
G1(G > 1) = 0;
