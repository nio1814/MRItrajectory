classdef PhantomSimulator < handle
properties
FOV = [];
resolution = [];
kScale = [];
FT = [];
phantom = [];
dimensions = 0;
end

methods
% 		Constructor
function obj = PhantomSimulator(numDimensions, phantomSize, translation)
  if(nargin<2)
		assert(numDimensions>1, 'At least 2 dimensional');
		assert(numDimensions<=3, 'At most 3 dimensional');
		phantomSize = [20, 24, 18];
  end
  if(nargin<3) 
		translation = zeros(numDimensions, 1); 
	end
%					Get scaling factor
	phantomSizeScaleOriginal = [1.38 1.84 1.62]; % Phantom size in original units
	phantomSizeScale = phantomSize(1:numDimensions)'./phantomSizeScaleOriginal(1:numDimensions);	%scaling factor to cm
	
	switch(numDimensions)
		case 2
			disp('2D not yet implemented');
		case 3
			[~, E] = shepplogan3D;
			E(:,[8:10]) = E(:,[8:10])*pi/180;	% scale angles to radians
			% Shift the object off center
			E(:,5) = E(:,5)+translation(1); 
			E(:,6) = E(:,6)+translation(2); 
			E(:,7) = E(:,7)+translation(3); 

%            % Kill outer bright ring
%            E = E(2:end,:);
%            E(1,1) = 0.2;
%  
%            % Add wide flat disks
%            dw = 0.03;
%            Eadd =[0.3, dw, 0.8, 0.8, 0,    0, 0, 0, 0, 0;...
%                   0.3, dw, 0.8, 0.8, 0.4,  0, 0, 0, 0, 0;...
%                   0.3, dw, 0.8, 0.8, -0.4, 0, 0, 0, 0, 0;...
%                   0.3, 0.8, dw, 0.8, 0, 0, 0, 0, 0, 0;...
%                   0.3, 0.8, dw, 0.8, 0, 0.4, 0, 0, 0, 0;...
%                   0.3, 0.8, dw, 0.8, 0, -0.4, 0, 0, 0, 0;...
%                   0.3, 0.8, 0.8, dw, 0, 0, 0, 0, 0, 0;...
%                   0.3, 0.8, 0.8, dw, 0, 0, 0.4, 0, 0, 0;...
%                   0.3, 0.8, 0.8, dw, 0, 0, -0.4, 0, 0, 0;...
%                   ]; 
%            E = cat(1,E,Eadd);

% 						Scale to specified size
		for n=1:3
			idx = [1,4] + n;
			E(:,idx) = phantomSizeScale(n)*E(:,idx);
		end
	  obj.phantom = SheppLogan3D(E);
  end
end
end
end
