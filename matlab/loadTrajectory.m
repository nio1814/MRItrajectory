function trajectory = loadTrajectory(filename, endian)

if(nargin<2)
  endian = 'l';
end
  
fileID = fopen(filename, 'r', endian);
if(fileID==-1)
   error("saveTrajectory: Error opening %s for read\n", filename);
end

fileVersion = fread(fileID, 1, 'int32');
dimensions = fread(fileID, 1, 'int32');
trajectory = [];
trajectory.imageDimensions = fread(fileID, dimensions, 'int32');
trajectory.spatialResolution = fread(fileID, dimensions, 'single');
trajectory.fieldOfView = fread(fileID, dimensions, 'single');
readouts = fread(fileID, 1, 'int32');
bases = fread(fileID, 1, 'int32');
trajectory.maxGradientAmplitude = fread(fileID, 1, 'single');
trajectory.maxReadoutGradientAmplitude = fread(fileID, 1, 'single');
trajectory.maxSlewRate = fread(fileID, 1, 'single');
waveformPoints = fread(fileID, 1, 'int32');
readoutPoints = fread(fileID, 1, 'int32');
trajectory.samplingInterval = fread(fileID, 1, 'single');
storage = fread(fileID, 1, 'int32');

points = fread(fileID, 1, 'ulong');
if(points)
  trajectory.gradientWaveforms = fread(fileID, points, 'single');
  waveforms = points/(waveformPoints*dimensions);
  trajectory.gradientWaveforms = reshape(trajectory.gradientWaveforms, [waveformPoints dimensions waveforms]);  
end

points = fread(fileID, 1, 'ulong');
if(points)
  trajectory.gradientWaveformsShort = fread(fileID, points, 'short');
  waveforms = points/(waveformPoints*dimensions);
  trajectory.gradientWaveformsShort = reshape(trajectory.gradientWaveformsShort, [waveformPoints dimensions waveforms]);  
end

points = fread(fileID, 1, 'ulong');
if(points)
  trajectory.kSpaceCoordinates = fread(fileID, points, 'single');
  waveforms = points/(readoutPoints*dimensions);
  trajectory.kSpaceCoordinates = reshape(trajectory.kSpaceCoordinates, [readoutPoints dimensions waveforms]);  
end

points = fread(fileID, 1, 'ulong');
if(points)
  waveforms = points/(readoutPoints);
  trajectory.densityCompensation = fread(fileID, points, 'single');
  trajectory.densityCompensation = reshape(trajectory.densityCompensation, [readoutPoints waveforms]);  
end

fclose(fileID);