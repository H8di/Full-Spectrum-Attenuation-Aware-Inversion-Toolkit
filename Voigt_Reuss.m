function [K_V, K_R, Mu_V, Mu_R]=Voigt_Reuss(fi, F, K, Mu, Kf)

 %F:volume fraction of constituation 
   %K:Bulk modulus of constituation 
     %Mu:Shear modulus of constituation 
       %fi:porosity in percent
         %Kf:fluid Bulk muduli
           %Muf:fluid shear moduli

%% Voigt_Ruess Computation Band Elastic Moduli

%Voigt
K_V=(fi*Kf)+((1-fi)*sum((F.*K)));
Mu_V=((1-fi)*sum((F.*Mu)));

%Rouess
K_R=((fi./Kf)+((1-fi).*sum(F./K))).^-1;
if fi~=0
    Mu_R=0;
else
Mu_R=(((1-fi)*sum(F./Mu))).^-1;
end

end


