function [sino] = sino2Dastra(phantom, angles, detsize, N)

proj_geom = astra_create_proj_geom('parallel', 1, detsize, angles);
vol_geom = astra_create_vol_geom(N,N);

[sinogram_id, sino] = astra_create_sino_cuda(phantom, proj_geom, vol_geom);
astra_mex_data2d('delete', sinogram_id);

