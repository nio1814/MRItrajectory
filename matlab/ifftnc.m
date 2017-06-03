function a = ifftnc(a, N, doShift, doScale);
% Written by Corey Baron
% Calculates the multidimensional ifft of a matrix with DC at the
% center of the matrix. 
%
% Inputs:
%   a: input matrix
%   N: do fft's along first "N" dim
% 
% Outputs:
%   A: output matrix
% 

if(nargin<4) || isempty(doScale)
  doScale = 1;
end
if(nargin<3) || isempty(doShift)
  doShift = 1;
end
if(nargin<2) || isempty(N)
  N = length(size(a));
end

if (N > length(size(a)))
  N = length(size(a))
end

scale = 1;
if doShift
  for n=1:N
    scale = scale*sqrt(size(a,n));
    a = ifftshift(a,n);
    a = ifft(a,[],n);
    a = fftshift(a,n);
  end
else
  for n=1:N
    scale = scale*sqrt(size(a,n));
    a = ifft(a,[],n);
  end
end

if doScale
  a = a*scale;
end	
