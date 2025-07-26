function tau=tau_calc(semiax,Mu,nu,eta,k,zi)
%eta:Viscosity
%cv=volume of crack
%k:permeability
%zi:grain size
%Mu:shear Modulus
%ar:aspect ratio
%nu:poisson ratio
%a:crack radius
%cf=Kf:fluid bulk modulus
% Kc=(pi*Mu*Kf*ar)/(2*(1-nu));
% Qc=(pi*Mu*ar)/(2*(1-nu));
% cv=(4*pi*(a^3)*ar)/3;
% 
% if ar>10^-2
% tau=(eta*cv*(1+Kc))/(6*k*zi*Qc);
% else

a=max(semiax) ;  
tau=(4*eta*(a^3)*(1-nu))/(9*((trace(k))/3)*zi*Mu);

end