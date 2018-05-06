function m = img(p, x, y, z);
% Description
%	Returns the pixel value at the specified image location in cm

if(p.dimensions==3)
	m = p.phantom.imageDomainSignal(x,y,z);
else
	m = p.phantom.imageDomainSignal(x,y);
end

m = reshape(m, size(x));