function [K_ave, Mu_ave]=Hash_Sht(pt, F, K , Mu, Kf, Muf)

% F:volume fraction of constituation 
% K:Bulk modulus of constituation 
% Mu:Shear modulus of constituation 
% pt:porosity in percent
% Kf:fluid Bulk muduli
% Muf:fluid shear moduli

Z_upper=(max(Mu)/6)*(9*max(K)+8*max(Mu))/(max(K)+2*max(Mu));
K_upper=((((pt/(Kf+((4/3)*max(Mu))))+((1-pt)*(sum(F/(K+((4/3)*max(Mu))))))))^-1)-(4/3*max(Mu));

if pt>0
Z_lower=(Muf/6)*(9*Kf+8*Muf)/(Kf+2*Muf);
K_lower=( ((pt/Kf)+((1-pt)*sum(F/K))) )^-1;
elseif pt==0
Z_lower=(min(Mu)/6)*(9*min(K)+8*min(Mu))/(min(K)+2*min(Mu));
K_lower=((1-pt)*sum(F/K))^-1;  
end

Mu_upper=(((pt/Z_upper)+((1-pt)*sum(F/(Mu+Z_upper))))^-1)-Z_upper;
Mu_lower=(((pt/Z_lower)+((1-pt)*sum(F/(Mu+Z_lower))))^-1)-Z_lower;

K_ave=(K_upper+K_lower)/2;
Mu_ave=(Mu_upper+Mu_lower)/2;
 end
