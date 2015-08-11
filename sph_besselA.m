function [x1] = sph_besselA(n,z) 
if abs(z) < 0.000000000001 x1 = ones(size(z)); 
else 
x1 = besselj(n+1./2,z).*sqrt(pi./(2*z));
end
return


