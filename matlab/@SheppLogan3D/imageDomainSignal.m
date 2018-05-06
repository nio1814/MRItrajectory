% /**
%      * returning real value of the image intensity at (x,y,z).
%      *
%      */
%     public double ImageDomainSignal(double x, double y, double z){
function signal = imageDomainSignal(sl, x, y, z)

%         double[] r = {x,y,z};
npts = length(x(:));
r = [x(:) y(:) z(:)]';
        
%         double signal = 0.0; 
signal = zeros(1,npts);
        
%         double[] p = new double[3];
        
%         double sum = 0.0;
        
%         for(int i=0; i<this.NE; i++){ // loop through each of the ellipsoids
for e=1:sl.numEllipsoids
        
                 
%              p = MAT.dot(RT[i],new double[] {r[0]-d[i][0],r[1]-d[i][1],r[2]-d[i][2]});
	p = sl.RT(:,:,e)*(r-(sl.d(e,:)'*ones(1,npts)));
            
%              sum = Math.pow(p[0]/abc[i][0],2) + Math.pow(p[1]/abc[i][1],2) + Math.pow(p[2]/abc[i][2],2);
	sm = sum((diag(1./[sl.a(e) sl.b(e) sl.c(e)])*p).^2, 1);
        
%              signal += (sum<=1.0)?rho[i]:0;
	idx = sm<=1;
	signal(idx) = signal(idx) + sl.rho(e);
        
%         }
end
        
%         return signal;
%     }
		end