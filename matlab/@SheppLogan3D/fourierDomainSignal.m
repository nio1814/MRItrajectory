%     /**
%      * Returns complex k-space signal evaluated at ( kx, ky, kz)
%      */
function signal = fourierDomainSignal(sl, kx, ky, kz)
     
npts = length(kx(:));
k = [kx(:) ky(:) kz(:)]';
clear kx;
clear ky;
clear kz;

signal = zeros(1,npts);
for e=1:sl.NE
    clear idx0;
    clear idx1;
    clear idx2;
    clear signalAdd0;
    clear signalAdd1;
    clear signalAdd2;
    clear argK;
    
	abc = sl.a(e)*sl.b(e)*sl.c(e);

	K = sqrt(sum((diag([sl.a(e) sl.b(e) sl.c(e)])*sl.RT(:,:,e)*k).^2, 1));

    idx0 = find(K==0);
    idx1 = find((0<K) .* (K<=.002));
    idx2 = find(.002<K);
    
    if(~isempty(idx2))
        argK = 2*pi*K(idx2);
    end
    
    if(norm(sl.d(e,:))==0.0)
    % 	k=0
        if(~isempty(idx0))
                signalAdd0 = (4/3)*pi* sl.rho(e)*abc;
        end
        
 %   	k<=.002
        
        if(~isempty(idx1))
                signalAdd1 = sl.rho(e)*abc*4.1887902047863905 - 16.5366808961599*K(idx1).^2 + 23.315785507450016*K(idx1).^4;
        end

        if(~isempty(idx2))
   % 					 .002<k
            signalAdd2 = sl.rho(e)*abc*(sin(argK)-argK.*cos(argK))./(2*(pi^2)*K(idx2).^3);
        end
    else
% 	k=0
		kd = sl.d(e,:)*k;
		arg = 2.0 * pi * kd;
		carg = cos(arg);
		sarg = sin(arg);
        
        if(~isempty(idx0))
	

		signalAdd0  = (4/3)*pi* sl.rho(e)*abc * (carg(idx0) - 1i*sarg(idx0));
        end
        
        if(~isempty(idx1))
% 	k<=.002

		signalAdd1 = sl.rho(e)*abc*4.1887902047863905 - 16.5366808961599*K(idx1).^2 + 23.315785507450016*K(idx1).^4 .* (carg(idx1) -1i*sarg(idx1));
        end
        
        if(~isempty(idx2))
% 	.002<k

signalAdd2 = sl.rho(e)*abc*(sin(argK)-argK.*cos(argK))./(2*(pi^2)*K(idx2).^3).*(carg(idx2)-1i*sarg(idx2));
        end
    end

if(~isempty(idx0))
        signal(idx0) = signal(idx0) + signalAdd0;
end

if(~isempty(idx1))
    signal(idx1) = signal(idx1) + signalAdd1;
end
    signal(idx2) = signal(idx2) + signalAdd2;	 
end