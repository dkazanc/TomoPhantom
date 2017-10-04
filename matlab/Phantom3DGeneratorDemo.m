% Script to generate 3D analytical phantoms
% If one needs to modify/add phantoms just edit Phantom3DLibrary.dat

close all;clc;clear all;
% adding paths
addpath('../models/'); addpath('supp/');

ModelNo = 01; % Select a model
% Define phantom dimensions
N = 512; % x-y-z size (cubic image)

% generate the 3D phantom:
[G] = buildPhantom3D(ModelNo,N);

% check the cenral slice
figure(1); imshow(G(:,:,round(0.5*N)), []);

% see the whole 3D Phantom
figure(2);
filename = strcat('ModelNo',num2str(ModelNo));
counter = 1;
for i = 1:N
    imshow(G(:,:,i), [0 1]);
    pause(0.01);
    
    % write tiff images
%     IM = im2uint16(G(:,:,i));
%     setStrNo = num2str(counter);
%     if (counter < 10)
%         filename_save = strcat(filename,'_','000',setStrNo, '.tiff');
%     elseif ((counter >= 10) &&  (counter < 100))
%         filename_save = strcat(filename,'_','00',setStrNo,'.tiff');
%     elseif ((counter >= 100) &&  (counter < 1000))
%         filename_save = strcat(filename,'_','0',setStrNo, '.tiff');
%     else
%         filename_save = strcat(filename,'_',setStrNo, '.tiff');
%     end
%     imwrite(IM,filename_save,'tiff');
%     counter = counter + 1;
end
%%