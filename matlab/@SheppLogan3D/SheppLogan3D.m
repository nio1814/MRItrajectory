% Shepp Logan 3D phantom
%  *
%  * Three-dimensional Shepp-Logan Phantom in both the Fourier and image domains.
%  *
%  * This is a class called SheppLogan3D. It can be used to generate Fourier domain 
%  * signal (or k-space) as well as image domain signal based on a 3-dimensional 
%  * analytical phantom in both domains. Please refer to the source code or 
%  * the article referenced below for further information.
%  * 
%  * <br>
%  * Based on 
%  * <br>
%  * Koay CG, Sarlls JE, &#214zarslan E. 
%  * Three Dimensional Analytical Magnetic Resonance Imaging Phantom in the Fourier Domain. Magn Reson Med. 58: 430-436 (2007)
%  * <br> 
%  * for further information.
%  * @see <a href=http://dx.doi.org/10.1002/mrm.21292>Ref</a>
%  * @author  Cheng Guan Koay
%  * @since 07/25/2007
%  *
%  * @version &#8722&#8734.
classdef SheppLogan3D
	properties
		numEllipsoids = 0;	% number of ellipsoids

		RT = [];
		d = [];

%         	The length of the principal axes.
		a = [];
		b = [];
		c = [];

%         	signal intensity.
		rho = [];
		
		% Rotation angles
		theta = [];
		phi = [];
		psi = [];
	end
    
    methods
		
	function obj = SheppLogan3D(E)
        slColA = 1;
		slCola = 2;
		slColb = 3;
		slColc = 4;
		slColx0 = 5;
		slColy0 = 6;
		slColz0 = 7;
		slColphi = 8;
		slColtheta = 9;
		slColpsi = 10;
		
		if(nargin==0)
			E =  [       0,       0,       0,     0.69,    0.92,     0.9,              0,      0,    0      2. ;
              0,       0,       0,   0.6624,   0.874,    0.88,              0,      0,    0,    -0.8 ;
          -0.22,      0.,   -0.25,     0.41,    0.16,    0.21, (3*pi)/5.,      0,    0,    -0.2 ;
           0.22,      0.,   -0.25,     0.31,    0.11,    0.22, (2*pi)/5.,      0,    0,    -0.2 ;
              0,    0.35,   -0.25,     0.21,    0.25,     0.5,              0,      0,    0,     0.2 ;
              0,     0.1,   -0.25,    0.046,   0.046,   0.046,              0,      0,    0,     0.2 ;
          -0.08,   -0.65,   -0.25,    0.046,   0.023,    0.02,              0,      0,    0,     0.1 ;
           0.06,   -0.65,   -0.25,    0.046,   0.023,    0.02,              0,      0,    0,     0.1 ;
           0.06,  -0.105,   0.625,    0.056,    0.04,     0.1,     pi/2.,      0,    0,     0.2 ;
             0.,     0.1,   0.625,    0.056,   0.056,     0.1,     pi/2.,      0,    0,    -0.2 ];
			
% 			Reformat E matrix
			 E = E(:, [10 4 5 6 1 2 3 8 7 9]);
		end
		
		obj.numEllipsoids = size(E,1);

		obj.d = E(:,[slColx0 slColy0 slColz0]);	% ellipsoid offset 
        

		obj.a = E(:,slCola);
		obj.b = E(:,slColb);
        obj.c = E(:,slColc);

		obj.theta = E(:, slColtheta);
		obj.phi = E(:, slColphi);
		obj.psi = E(:, slColpsi);
			for e=1:obj.numEllipsoids
			   cphi = cos(obj.phi(e));
			   sphi = sin(obj.phi(e));
			   ctheta = cos(obj.theta(e));
			   stheta = sin(obj.theta(e));
			   cpsi = cos(obj.psi(e));
			   spsi = sin(obj.psi(e));
   
				obj.RT(:,:,e) = [cpsi*cphi-ctheta*sphi*spsi   cpsi*sphi+ctheta*cphi*spsi  spsi*stheta;
			        -spsi*cphi-ctheta*sphi*cpsi  -spsi*sphi+ctheta*cphi*cpsi cpsi*stheta;
				    stheta*sphi                  -stheta*cphi                ctheta];        
			end
         
			obj.rho = E(:,slColA);
		end
	end
    
end
	
% Description
% 	Calculate the rotation matrix about the x-axis
% Inputs
% 	t:	Rotation angle (radians)
% Output
% 	R:	Rotation matrix
function R = Rx(t)
	R = [1, 0, 0; 0, cos(t), -sin(t); 0, sin(t), cos(t)];
         
end
  
% Description
% 	Calculate the rotation matrix about the y-axis
% Inputs
% 	t:	Rotation angle (radians)
% Output
% 	R:	Rotation matrix
function R = Ry(t)
	R = [cos(t), 0, sin(t); 0, 1, 0; -sin(t), 0, cos(t)];
         
end
     
% Description
% 	Calculate the rotation matrix about the z-axis
% Inputs
% 	t:	Rotation angle (radians)
% Output
% 	R:	Rotation matrix
function R = Rz(t)

	R = [cos(t), -sin(t), 0; sin(t), cos(t), 0; 0, 0, 1];
         
end