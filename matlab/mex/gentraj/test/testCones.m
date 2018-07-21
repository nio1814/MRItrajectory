%% trajectory 1
addpath('../..')

params.fieldOfView = [28 28 14];
params.spatialResolution = [1 1 1]*1.2;
trajectory = gentraj('cones', params);
