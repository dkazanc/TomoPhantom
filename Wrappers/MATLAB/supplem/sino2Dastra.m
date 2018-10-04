function [sino] = sino2Dastra(phantom, angles, detsize, N, device)

proj_geom = astra_create_proj_geom('parallel', 1, detsize, angles);
vol_geom = astra_create_vol_geom(N,N);

if (strcmp(device, 'cpu') == 1)
    proj_id = astra_create_projector('linear', proj_geom, vol_geom);
    [sinogram_id, sino] = astra_create_sino(phantom, proj_id);
    astra_mex_data2d('delete', proj_id);
elseif (strcmp(device, 'gpu') == 1)
    [sinogram_id, sino] = astra_create_sino_cuda(phantom, proj_geom, vol_geom);
else
    error('Device is not selected, please choose between cpu and gpu')
end
astra_mex_data2d('delete', sinogram_id);

