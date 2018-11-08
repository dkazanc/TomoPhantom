function [reconstructon] = rec3Dastra(proj_data, angles, DetectorHeight, DetectorWidth)


proj_geom = astra_create_proj_geom('parallel3d', 1, 1, DetectorHeight, DetectorWidth, angles);
vol_geom = astra_create_vol_geom(DetectorHeight,DetectorHeight,DetectorHeight);


cfg = astra_struct('CGLS3D_CUDA');    
rec_id = astra_mex_data3d('create', '-vol', vol_geom);
sinogram_id = astra_mex_data3d('create', '-sino', proj_geom, proj_data);
alg_id = astra_mex_algorithm('create', cfg);

iterations = 10;
astra_mex_algorithm('iterate', alg_id, iterations);
%astra_mex_algorithm('run', alg_id);
 
% Get the result
reconstructon = astra_mex_data3d('get', rec_id);

% % Clean up. Note that GPU memory is tied up in the algorithm object,
% % and main RAM in the data objects.
astra_mex_algorithm('delete', alg_id);
astra_mex_data3d('delete', sinogram_id);
astra_mex_data3d('delete', rec_id);