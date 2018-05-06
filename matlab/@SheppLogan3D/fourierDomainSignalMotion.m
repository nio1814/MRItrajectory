function signal = fourierDomainSignalMotion(sl, kx, ky, kz, dmotion)
%     sl:	Shepp Logan phantom simulator object
Nread = size(kx,1);	% # readout points per acquisition
Nacq = size(kx,2);	% # total # of readouts
npts = length(kx(:));	% total # of points

%         double[] k = {kx,ky,kz};
k = [kx(:) ky(:) kz(:)]';

signal = zeros(1,npts);
for e=1:sl.NE   
abc = sl.a(e)*sl.b(e)*sl.c(e);
K = sqrt(sum((diag([sl.a(e) sl.b(e) sl.c(e)])*sl.RT(:,:,e)*k).^2, 1));
%             
%              arg = 2.0 * Math.PI * K;
argK = 2*pi*K;
cargK = cos(argK);
sargK = sin(argK);
%         

K2 = K.^2;
K3 = K.^3;
K4 = K.^4;

for na=1:Nacq
% 	if(mod(na,1000)==0)
% 	fprintf(1, '%d ', na);
% 	end
	n = [1:Nread] + Nread*(na-1);
	dm = sl.d(e,:) + dmotion(na,:);
	if(norm(dm)==0.0)
% 	k=0
		 signalAdd0 = (4/3)*pi* sl.rho(e)*abc;
		 
		 temp = 4.1887902047863905 - 16.5366808961599*K2(n) + 23.315785507450016*K4(n);
		 signalAdd1 = sl.rho(e)*abc*temp;
		 
		temp = (sargK(n)-argK(n).*cargK(n))./(2*(pi.^2)*K3(n));
		signalAdd2 = sl.rho(e)*abc*temp;
	else
		kd = dm*k(:,n);
		arg = 2.0 * pi * kd;
		carg = cos(arg);
		sarg = sin(arg);
	
%                      double temp = (4./3.)*Math.PI* rho[i]*abc[i][0]*abc[i][1]*abc[i][2];
	temp = (4/3)*pi* sl.rho(e)*abc;
%                      signal[0] += temp * Math.cos(2.0 * Math.PI * kd);
%                      signal[1] -= temp * Math.sin(2.0 * Math.PI * kd);
	signalAdd0  = temp * (carg - 1i*sarg);
	
% 	k<=.002
%                      double kd = MAT.dot(k,d[i]);
%                      double temp1 = 4.1887902047863905 - 16.5366808961599*Math.pow(K,2) + 23.315785507450016*Math.pow(K,4);
	temp1 = 4.1887902047863905 - 16.5366808961599*K2(n) + 23.315785507450016*K4(n);
%                      double temp2 = rho[i]*abc[i][0]*abc[i][1]*abc[i][2]*temp1;
temp2 = sl.rho(e)*abc*temp1;
%                      
%                      signal[0] += temp2 * Math.cos(2.0 * Math.PI * kd);
%                      signal[1] -= temp2 * Math.sin(2.0 * Math.PI * kd);
signalAdd1 = temp2 .* (carg -1i*sin(arg));

% 	.002<k
%                      double kd = MAT.dot(k,d[i]);
%                      double temp = Math.sin(arg)-arg*Math.cos(arg);
%                             temp /= (2.0*Math.pow(Math.PI,2)*Math.pow(K,3));
	
%                             temp *= rho[i]*abc[i][0]*abc[i][1]*abc[i][2];
temp = sl.rho(e)*abc*(sargK(n)-argK(n).*cargK(n))./(2*(pi^2)*K3(n));
%                  
%                      signal[0] += temp * Math.cos(2.0 * Math.PI * kd);
%                      signal[1] -= temp * Math.sin(2.0 * Math.PI * kd);
	signalAdd2 = temp.*(carg-1i*sarg);
	end

% 		Find logical array
		idx0 = K(n)==0;
		 idx1 = logical((0<K(n)) .* (K(n)<=.002));
		 idx2 = (.002<K(n));
		 
% 		 Translage to indeces
		idx0s = find(idx0)+n(1)-1;
		idx1s = find(idx1)+n(1)-1;
		idx2s = find(idx2)+n(1)-1;
		 
		 if(length(signalAdd0)>1)
			 signal(idx0s) = signal(idx0s) + signalAdd0(idx0);
		 else
			signal(idx0s) = signal(idx0s) + signalAdd0;
		 end
		 
		 if(~isempty(idx1))
			signal(idx1s) = signal(idx1s) + signalAdd1(idx1);
		 end
		 
		 if(~isempty(idx2))
			signal(idx2s) = signal(idx2s) + signalAdd2(idx2);
		 end
	end
end
end