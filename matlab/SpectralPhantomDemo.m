% spectral phantom

close all;clc;clear all;
% adding paths
addpath('../models/'); addpath('supp/'); 

ModelNo = 43; % Select a model (0 - 42  )
% Define phantom dimension
N = 512; % x-y size (squared image)

% generate the 2D phantom:
[G] = buildPhantom(ModelNo,N);
figure(1); imshow(G, [0 1]);

% create 4 phantoms with dedicated materials
G(G ==1.3) = 0.8;
G1 = zeros(N,N);
G2 = zeros(N,N);
G3 = zeros(N,N);
G4 = zeros(N,N);

G1(G == 0.3) = 1;
G2(G == 0.8) = 2;
G3(G == 1.1) = 3;
G4(G == 1.2) = 4;
G1(G > 1) = 0;
