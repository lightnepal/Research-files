function mine (id11,id22)
id1 = str2num(id11)
id2 = str2num(id22)



load Alm22.mat;
N = 26;
dr = 50/26;%corrected 

[x y z] = ndgrid(-N/2:N/2,-N/2:N/2,-N/2:N/2);
 r = [x(:) y(:) z(:)];
 rx = dr*r(:,1);
 ry = dr*r(:,2);
 rz = dr*r(:,3);
 [phi,theta,r] = cart2sph(rx,ry,rz);
 theta = pi/2.0-theta;


lmax = 81;

dq = 1/(3* 50);

rmax=43.4;
rmin=0.00001;
npoints = 4000;
dr=(rmax-rmin)/npoints;
x1 = 2*pi*rmin:2*pi*dr:2*pi*rmax;

ii = id1

ij = id2 
 

for l = 0:lmax
    
q1 = dq * (ii -1);

q2 = dq * (ij -1);

bessarr1 = sph_besselA(l,q1*x1);

bessarr2 = sph_besselA(l,q2*x1);
g1(ii,l+1,1:(N+1)^3) = interp1(x1,bessarr1,2*pi*r);
g2(ij,l+1,1:(N+1)^3) = interp1(x1,bessarr2,2*pi*r);
end

for k = 1:(N+1)^3
for l= 0:lmax
    for m = -l:l
        h1(l+1,m+l+1,k) = compute_ylm(l,m,theta(k),phi(k));
    end
end
end



for ik = 1:27 %1:11  %15:20
            
            
            
   
    for l = 0:lmax
    q3 = (ik + 54 -1)*dq

    bessarr3 = sph_besselA(l,q3*x1);

    g3(ik + 54,l+1,1:(N+1)^3) = interp1(x1,bessarr3,2*pi*r);
    
    end
end

for ik = 1:27
            
    for k= 1 : 19683 % general k
    
    sum4=0;
    
    for l = 0:lmax
        
        for m = -l:l
                    
                    
                    
            
            sum4 = sum4 + 4 * pi * (((j)^-l * ((Alm22(ii,l+1, m+l+1)* g1(ii,l+1,k)) + (Alm22(ij,l+1, m+l+1)* g2(ij,l+1,k)) + (Alm22(ik+54,l+1, m+l+1) * g3(ik+54,l+1,k))) * (h1(l+1,m+l+1,k))) + ((conj((j)^-l)* (conj(Alm22(ii,l+1, m+l+1)* g1(ii,l+1,k)) + conj((Alm22(ij,l+1, m+l+1)* g2(ij,l+1,k))) + conj((Alm22(ik+54,l+1, m+l+1) * g3(ik+54,l+1,k))))) * conj(h1(l+1,m+l+1,k))))
            
            
        end
                    
                    
                    
    end
    sum4;
    M94(ik,k) = sum4;
    end
    
end
            
        

stri1 = int2str(id1)
stri2 = int2str(id2)
namefile = strcat('nn',stri1,'_',stri2,'.mat')

save(namefile,'M94');

end


