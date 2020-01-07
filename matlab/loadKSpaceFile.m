function [trajectoryCoordinates, densityCompensation] = loadKSpaceFile(filePath, numDimensions, endian, numReadouts, numReadoutPoints)
%    LOADKSPACEFILE  Load a k-space tractory and density compensation weights from a file.
%    
%    Inputs:
%      filePath - The file path to the trajectory file.
%      numDimensions - The number of dimensions in the trajectory.
%      endian - The byte format of the saved file ('b' for big, 'l' for little).
%      numReadouts - (Optional) The number of readouts of the trajectory.
%      numReadoutPoints - (Required with `numReadouts`) The number of readout points per readout.
    
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
