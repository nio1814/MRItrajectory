function signal = acquire(p, k, varargin)
% Inputs
% 	p:	Phantom simulator object
% 	kx: 2D or 3D k-space coordinates [readout,acquisition,axis]
% 	 or
% 	kx: 2D k-space coordinates [readout,acquisition] (real) kx, (imag) ky
% 	 or
% 	kx:	x k-space coordinate (1/cm)
% 	ky:	y k-space coordinate (1/cm)
% 	kz:	z k-space coordinate (1/cm)
% 
% Output
% 	signal:	Acquired signal

doReformat = 0;

ip = inputParser;

ip.addOptional('motion', []);
ip.addOptional('tr', []);
ip.addOptional('te', []);
ip.addOptional('flip', []);
ip.addOptional('noise', []);
ip.addOptional('time', []);
ip.addOptional('bits', 16);
ip.addOptional('maxsig', []);
ip.parse(varargin{:});

doMotion = ~isempty(ip.Results.motion);

if(doMotion)
    d = ip.Results.motion;
end

% 	Assume 2D trajectory
doReformat = 1;
N(1) = size(k,1);	% Readout length
N(2) = size(k,2);	% Number of acquisitions

if(isreal(k))
% 		Trajectory input as [readout,acquisition,axis]
	k_x = k(:,:,1);
	if(p.dimensions>1)
		k_y = k(:,:,2);
	end
	
	if(p.dimensions>2)
		k_z = k(:,:,3);
	end
else
% 		Trajectory input as single complex array
	k = k;
	k_y = imag(k);
	k_x = real(k);
end

% Simulate signal
if(doMotion)
	switch(p.dimensions)
		case 2
			signal = p.phantom.fourierDomainSignalMotion(d, k_x, k_y);
		case 3
			signal = p.phantom.fourierDomainSignalMotion(k_x, k_y, k_z, d);
	end
else
	switch(p.dimensions)
		case 1
			signal = p.phantom.fourierDomainSignal(k_x);
		case 2
			signal = p.phantom.fourierDomainSignal(k_x, k_y);
		case 3
			signal = p.phantom.fourierDomainSignal(k_x, k_y, k_z);
	end
end

% Reshape signal to trajectory dimensions
if(doReformat)
	signal = reshape(signal, N);
end

% If specified, add noise
if(~isempty(ip.Results.noise))
	noiseSignal = sqrt(ip.Results.noise)*(randn(size(signal)) + i*randn(size(signal)));
else
	noiseSignal = 0*signal;
end

% Quantize to n bit number to simulate scanner
B = ip.Results.bits-1;
if(isempty(ip.Results.maxsig))
	maxsig = (2^B);
else
	maxsig = ip.Results.maxsig;
end

signal = round(maxsig*signal/max(abs([real(signal(:)); imag(signal(:))])) + noiseSignal);

end