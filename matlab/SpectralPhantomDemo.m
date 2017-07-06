%spectral phantom

close all;clc;clear all;
% adding paths
addpath('../models/'); addpath('supp/'); 

ModelNo = 43; % Select a model (0 - 42  )
% Define phantom dimension
N = 400; % x-y size (squared image)

% generate the 2D phantom:
[G] = buildPhantom(ModelNo,N);
G(G ==1.3) = 0.8;
G(G == 0.3) = 1;
G(G == 0.8) = 2;
G(G == 1.1) = 3;
G(G == 1.2) = 4;
figure(1); imshow(G, [0 4]);