function [trajectoryCoordinates, densityCompensation] = loadKSpaceFile(filePath, numDimensions, endian, numReadouts, numReadoutPoints)
  if(nargin < 3)
    endian = 'l';
  end
  
  fileID = fopen(filePath, 'r', endian);
  if(fileID == -1)
    error('Error opening file %s', filePath);
  end
  
  data = fread(fileID, inf, 'float32');
  fclose(fileID);
  
  if(nargin < 4)
    data = reshape(data, numDimensions + 1, []);
    trajectoryCoordinates = data(1:numDimensions, :);
    densityCompensation = data(end, :);
  else
    data = reshape(data, numDimensions + 1, numReadoutPoints, numReadouts);
    trajectoryCoordinates = data(1:numDimensions, :, :);
    densityCompensation = data(end, :, :);
  end
