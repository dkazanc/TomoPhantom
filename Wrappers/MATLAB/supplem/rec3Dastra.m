function [reconstructon] = rec3Dastra(proj_data, angles, detU, detV, N, device)


proj_geom = astra_create_proj_geom('parallel3d', 1, 1, detU, detV, angles);
vol_geom = astra_create_vol_geom(N,N,N);

if (strcmp(device, 'cpu') == 1)
    cfg = astra_struct('FBP');
    proj_id = astra_create_projector('strip', proj_geom, vol_geom);
    cfg.ProjectorId = proj_id;
    fprintf('%s \n', 'CPU FBP...');
elseif (strcmp(device, 'gpu') == 1)    
    cfg = astra_struct('CGLS3D_CUDA');
    fprintf('%s \n', 'GPU FBP...');
else
    error('Device is not selected, please choose between cpu and gpu')    
end
    
rec_id = astra_mex_data3d('create', '-vol', vol_geom);
sinogram_id = astra_mex_data3d('create', '-sino', proj_geom, proj_data);

cfg.FilterType = 'Ram-Lak';
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;

alg_id = astra_mex_algorithm('create', cfg);

astra_mex_algorithm('iterate', alg_id, 10);
%astra_mex_algorithm('run', alg_id);
 
% Get the result
reconstructon = astra_mex_data3d('get', rec_id);

% % Clean up. Note that GPU memory is tied up in the algorithm object,
% % and main RAM in the data objects.
astra_mex_algorithm('delete', alg_id);
astra_mex_data3d('delete', sinogram_id);
astra_mex_data3d('delete', rec_id);

if (strcmp(device, 'cpu') == 1)
astra_mex_data3d('delete', proj_id);
end