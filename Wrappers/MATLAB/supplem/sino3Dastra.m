function [proj3D_astra] = sino3Dastra(phantom, angles, DetectorHeight, DetectorWidth)

proj_geom = astra_create_proj_geom('parallel3d', 1, 1, DetectorHeight, DetectorWidth, angles*pi/180);
vol_geom = astra_create_vol_geom(DetectorHeight, DetectorHeight, DetectorHeight);

tic; [sinogram_id, proj3D_astra] = astra_create_sino3d_cuda(phantom, proj_geom, vol_geom); toc;
astra_mex_data3d('delete', sinogram_id);