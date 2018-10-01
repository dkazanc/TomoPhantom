% GPLv3 license (ASTRA toolbox)
% Note that the TomoPhantom package is released under Apache License, Version 2.0

% Script to generate 4D analytical phantoms (3D + time)
% If one needs to modify/add phantoms, please edit Phantom3DLibrary.dat
% >>>> Prerequisites: ASTRA toolbox, if one needs to do reconstruction <<<<<

close all;clc;clear;
% adding paths
fsep = '/';
pathtoModels = sprintf(['..' fsep 'functions' fsep 'models' fsep], 1i);
addpath(pathtoModels);
addpath('compiled'); addpath('supplem'); 

ModelNo = 100; % Select a model
% Define phantom dimensions
N = 256; % x-y-z size (cubic image)

% generate 4D phantom (modify your PATH bellow):
curDir   = pwd;
mainDir  = fileparts(curDir);
pathtoLibrary = sprintf([fsep 'functions' fsep 'models' fsep 'Phantom3DLibrary.dat'], 1i);
pathTP = strcat(mainDir, pathtoLibrary); % path to TomoPhantom parameters file
[G] = TomoP3DModel(ModelNo,N,pathTP);

sliceM = round(0.5*N);
figure(1);
for i = 1:5
    imagesc(G(:,:,sliceM, i), [0 1]); daspect([1 1 1]); title('3D+t phantom'); colormap hot;   
    pause(0.1); 
end
%%
% another 3D + time model
ModelNo = 101; % Select a model
% Define phantom dimensions
N = 256; % x-y-z size (cubic image)

% generate 4D phantom (modify your PATH bellow):
curDir   = pwd;
mainDir  = fileparts(curDir);
pathtoLibrary = sprintf([fsep 'functions' fsep 'models' fsep 'Phantom3DLibrary.dat'], 1i);
pathTP = strcat(mainDir, pathtoLibrary); % path to TomoPhantom parameters file
[G] = TomoP3DModel(ModelNo,N,pathTP);

%%
sliceM = round(0.5*N);
figure(2);
ll = 10;
counter  = 1;
for ll = 1:5
    imagesc(G(:,:,sliceM, ll), [0 1]); daspect([1 1 1]); title('3D+t phantom'); colormap hot;   
    pause(0.1); 
end
%%

