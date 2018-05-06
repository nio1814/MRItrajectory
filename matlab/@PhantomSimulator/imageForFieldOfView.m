function m = imageForFieldOfView(simulator, fov, imageSize);
% Description
%	Returns an image for a specified field of view and resolution

numDimensions = length(fov);

assert(length(fov)==length(imageSize), fprintf('Number of field of view dimensions %d in do not match image size dimensions %d\n', length(fov), length(imageSize)));

for d=1:numDimensions
	N2 = floor(imageSize(d)/2);
	if(mod(imageSize,2))
		n{d} = [-N2:N2]/imageSize(d)*fov(d);
	else
		n{d} = [-N2+1:N2]/imageSize(d)*fov(d);
	end
end

switch(numDimensions)
	case 1
		x = n{1};
		m = simulator.phantom.imageDomainSignal(x);
	case 2
		[x y] = meshgrid(n{1}, n{2});
		m = simulator.phantom.imageDomainSignal(x,y);
	case 3
		[x y z] = meshgrid(n{1}, n{2}, n{3});
		m = simulator.phantom.imageDomainSignal(x,y,z);
end

m = reshape(m, size(x));