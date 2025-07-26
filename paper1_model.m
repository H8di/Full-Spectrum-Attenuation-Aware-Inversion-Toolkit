clear, close all
clc


iso='yes'; %iso or aniso model 
k_matrix = 77e9;
u_matrix = 32e9;
kw = 2.25e9;
ko = 0.8e9;
kg = 0.01e9;
sw = 1;
so = 0;
sg = 0;
sat=sw+so+sg;
rho_mineral = 2.71e3;%kg/m3
rho = 0.8e3;%kg/m3
rhw = 1e3;%kg/m3
rhg = 1;%kg/m3
eta_w = 1e-3;%viscosity
eta_o = 60e-3;
eta_g = 0.02e-3;
rho_fluid = so*rho+sw*rhw+sg*rhg;
theta = 0;
phi = 0;
tau = 1e-5;
pt = 0.005;
fr = 0;
v_mineral = 1;
c_minerals = elastictotensor(k_matrix,u_matrix);
c_matrix = elastictotensor(k_matrix,u_matrix);
a_mineral = 1;%meter
b_mineral = 1;
c_mineral = 1;
v_pore = [0.9 0.1];
uf = 0;
a_pore = [1e-6,1e-6];
b_pore = [1e-6,1e-6];
c_pore = [1e-6,1e-9];
orientation_type_m = "simple_orientation";
theta_m = 0;
phi_m = 0;
psi_m = 0;
mu_theta_m = 0;
sigma_theta_m = 0;
mu_phi_m = 0;
sigma_phi_m = 0;
alpha_theta_m = 0;
beta_theta_m = 0;
alpha_phi_m = 0;
beta_phi_m = 0;
aM_theta_m = 0;
aM_phi_m = 0;
orientation_type_f = {"simple_orientation","chaotic"};
theta_f = [0,0];
phi_f = [0,0];
psi_f = [0,0];
mu_theta_f = [0,0];
sigma_theta_f = [0,0];
mu_phi_f = [0,0];
sigma_phi_f = [0,0];
alpha_theta_f = [0,0];
beta_theta_f = [0,0];
alpha_phi_f = [0,0];
beta_phi_f = [0,0];
aM_theta_f = [0,0];
aM_phi_f = [0,0];
communicative = 'yes';
[c_dry]=t_matrix(iso,[],[],[],[],[],[],[],[],pt,fr,sat,...
 v_mineral,c_minerals,c_matrix,a_mineral,b_mineral,c_mineral,...
 v_pore,0,uf,a_pore,b_pore,c_pore,...
 orientation_type_m,...
 theta_m,phi_m,psi_m,mu_theta_m,sigma_theta_m,mu_phi_m,...
 sigma_phi_m,alpha_theta_m,beta_theta_m,alpha_phi_m,...
 beta_phi_m,aM_theta_m,aM_phi_m,...
 orientation_type_f,...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f,'no');

c_eff = c_dry;

d = inf;
x0 = eye(3,3)*9.8692326671601e-16;
perm = eye(3,3).*t_matrix_perm(pt,x0,orientation_type_f,v_pore,[a_pore',b_pore',c_pore'],...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f);
%%%%sens on oil sauration%%%%
counter = 0;
omega = logspace(0,9,50);
so = linspace(0,.6,50);
N_omega = length(omega);
N_so = length(so);
for i = 1:N_so
 sg = 0;
 sw = 1 - so(i);
 kf = sw*kw+so(i)*ko+sg*kg;
 eta = eta_w*sw + eta_g*sg + eta_o*so(i);
 rho_fluid = so(i)*rho+sw*rhw+sg*rhg;
 Rho(i)= rho_mineral*(1-pt)+rho_fluid*pt;
 for j = 1:N_omega
 d = inf;

 while d>0.1
 c_old = c_eff;
 [c_eff]=...
 t_matrix(iso,c_eff,Rho(i),theta,phi,omega(j),tau,eta,perm,pt,fr,sat,...
 v_mineral,c_minerals,c_matrix,a_mineral,b_mineral,c_mineral,...
 v_pore,kf,uf,a_pore,b_pore,c_pore,...
 orientation_type_m,...
 theta_m,phi_m,psi_m,mu_theta_m,sigma_theta_m,mu_phi_m,...
 sigma_phi_m,alpha_theta_m,beta_theta_m,alpha_phi_m,...
 beta_phi_m,aM_theta_m,aM_phi_m,...
 orientation_type_f,...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f,communicative);
 d = norm(c_old-c_eff,"fro");
 end
y = Christoffel(c2dto4d(c_eff),Rho(i),theta,phi);
Vp1(i,j) = y(1,1);
Vs11(i,j) = y(2,1);
Vs21(i,j) = y(3,1);
 end
end
Qp=real(Vp1.^2)./imag(Vp1.^2);
alpha1 = -1000./Qp;
figure(1)
surf(so,omega,real(Vp1'))
xlabel("so")
ylabel("\omega, Hz")
zlabel("Vp, km/s")
set(gca,'yscale','log','FontName','Times New Roman','FontSize',14)
colormap jet
figure(2)
surf(so,omega,real((Rho.*Vp1)'))
xlabel("so")
ylabel("\omega, Hz")
zlabel("AI")
set(gca,'yscale','log','FontName','Times New Roman','FontSize',14)
colormap jet
figure(3)
surf(so,omega,real(Vs11'))
xlabel("so")
ylabel("\omega, Hz")
zlabel("Vs_1, km/s")
set(gca,'yscale','log','FontName','Times New Roman','FontSize',14)
colormap jet
figure(4)
surf(so,omega,real(Vs21'))
xlabel("so")
ylabel("\omega, Hz")
zlabel("Vs_2, km/s")
set(gca,'yscale','log','FontName','Times New Roman','FontSize',14)
colormap jet
figure(5)
surf(so,omega,alpha1')
xlabel("S_{oil}")
ylabel("\omega, Hz")
zlabel("\alpha")
set(gca,'yscale','log','FontName','Times New Roman','FontSize',14)
colormap jet

%%%%sens on gas sauration%%%%
k_matrix = 77e9;
u_matrix = 32e9;
kw = 2.25e9;
ko = 0.8e9;
kg = 0.01e9;
sw = 1;
so = 0;
sg = 0;
rho_mineral = 2.71e3;
rho = 0.8e3;
rhw = 1e3;
rhg = 1;
eta_w = 1e-3;
eta_o = 60e-3;
eta_g = 0.02e-3;
rho_fluid = so*rho+sw*rhw+sg*rhg;
theta = 0;
phi = 0;
tau = 1e-5;
pt = 0.005;
rho= rho_mineral*(1-pt)+rho_fluid*pt;
fr = 0;
v_mineral = 1;
c_minerals = elastictotensor(k_matrix,u_matrix);
c_matrix = elastictotensor(k_matrix,u_matrix);
a_mineral = 1;
b_mineral = 1;
c_mineral = 1;
v_pore = [0.9 0.1];
uf = 0;
a_pore = [1e-6,1e-6];
b_pore = [1e-6,1e-6];
c_pore = [1e-6,1e-9];
orientation_type_m = "simple_orientation";
theta_m = 0;
phi_m = 0;
psi_m = 0;
mu_theta_m = 0;
sigma_theta_m = 0;
mu_phi_m = 0;
sigma_phi_m = 0;
alpha_theta_m = 0;
beta_theta_m = 0;
alpha_phi_m = 0;
beta_phi_m = 0;
aM_theta_m = 0;
aM_phi_m = 0;
orientation_type_f = {"simple_orientation","chaotic"};
theta_f = [0,0];
phi_f = [0,0];
psi_f = [0,0];
mu_theta_f = [0,0];
sigma_theta_f = [0,0];
mu_phi_f = [0,0];
sigma_phi_f = [0,0];
alpha_theta_f = [0,0];
beta_theta_f = [0,0];
alpha_phi_f = [0,0];
beta_phi_f = [0,0];
aM_theta_f = [0,0];
aM_phi_f = [0,0];
communicative = 'yes';
[c_dry ]=...
 t_matrix(iso,[],[],[],[],[],[],[],[],pt,fr,sat,...
 v_mineral,c_minerals,c_matrix,a_mineral,b_mineral,c_mineral,...
 v_pore,0,uf,a_pore,b_pore,c_pore,...
 orientation_type_m,...
 theta_m,phi_m,psi_m,mu_theta_m,sigma_theta_m,mu_phi_m,...
 sigma_phi_m,alpha_theta_m,beta_theta_m,alpha_phi_m,...
 beta_phi_m,aM_theta_m,aM_phi_m,...
 orientation_type_f,...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f,'no');
c_eff = c_dry;
x0 = eye(3,3)*9.8692326671601e-16;
perm = eye(3,3).*t_matrix_perm(pt,x0,orientation_type_f,v_pore,[a_pore',b_pore',c_pore'],...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f);
counter = 0;
omega = logspace(0,9,50);
sg = linspace(0,.6,50);
N_omega = length(omega);
N_sg = length(sg);
for i = 1:N_sg
 so = 0;
 sw = 1 - sg(i);
 kf = 1/(sw/kw+0/ko+sg(i)/kg);
 eta = eta_w*sw + eta_w*sg(i) + eta_o*so;
 rho_fluid = so*rho+sw*rhw+sg(i)*rhg;
 rho= rho_mineral*(1-pt)+rho_fluid*pt;
 for j = 1:N_omega
 d = inf;
%  c_bk = LowFrequency_BrownKorringa(c_dry,c_matrix,kf,pt);
 while d>0.1
 c_old = c_eff;
 [c_eff]=...
 t_matrix(iso,c_eff,rho,theta,phi,omega(j),tau,eta,perm,pt,fr,sat,...
 v_mineral,c_minerals,c_matrix,a_mineral,b_mineral,c_mineral,...
 v_pore,kf,uf,a_pore,b_pore,c_pore,...
 orientation_type_m,...
 theta_m,phi_m,psi_m,mu_theta_m,sigma_theta_m,mu_phi_m,...
 sigma_phi_m,alpha_theta_m,beta_theta_m,alpha_phi_m,...
 beta_phi_m,aM_theta_m,aM_phi_m,...
 orientation_type_f,...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f,communicative);
 d = norm(c_old-c_eff,"fro");
 end
y = Christoffel(c2dto4d(c_eff),rho,theta,phi);
Vp2(i,j) = y(1,1);
Vs12(i,j) = y(2,1);
Vs22(i,j) = y(3,1);
 end
end
Qp=real(Vp2.^2)./imag(Vp2.^2);
alpha2 = -1000./Qp;
figure(6)
surf(sg,omega,real(Vp2'))
xlabel("sg")
ylabel("\omega, Hz")
zlabel("Vp, km/s")
set(gca,'yscale','log','FontName','Times New Roman','FontSize',14)
colormap jet
figure(7)
surf(sg,omega,real((rho.*Vp2)'))
xlabel("sg")
ylabel("\omega, Hz")
zlabel("Vp, km/s")
set(gca,'yscale','log','FontName','Times New Roman','FontSize',14)
colormap jet
figure(8)
surf(sg,omega,real(Vs12'))
xlabel("sg")
ylabel("\omega, Hz")
zlabel("Vs_1, km/s")
set(gca,'yscale','log','FontName','Times New Roman','FontSize',14)
colormap jet
figure(9)
surf(sg,omega,real(Vs22'))
xlabel("sg")
ylabel("\omega, Hz")
zlabel("Vs_2, km/s")
set(gca,'yscale','log','FontName','Times New Roman','FontSize',14)
colormap jet
figure(10)
surf(sg,omega,alpha2')
xlabel("sg")
ylabel("\omega, Hz")
zlabel("\alpha, km/s")
set(gca,'yscale','log','FontName','Times New Roman','FontSize',14)
colormap jet
%%%%%%%relaxation time%%%%%%%%%%
k_matrix = 77e9;
u_matrix = 32e9;
kw = 2.25e9;
ko = 0.8e9;
kg = 0.01e9;
sw = 1;
so = 0;
sg = 0;
rho_mineral = 2.71e3;
rho = 0.8e3;
rhw = 1e3;
rhg = 1;
eta_w = 1e-3;
eta_o = 60e-3;
eta_g = 0.02e-3;
rho_fluid = so*rho+sw*rhw+sg*rhg;
theta = 0;
phi = 0;
tau = 1e-5;
pt = 0.005;
rho= rho_mineral*(1-pt)+rho_fluid*pt;
fr = 0;
v_mineral = 1;
c_minerals = elastictotensor(k_matrix,u_matrix);
c_matrix = elastictotensor(k_matrix,u_matrix);
a_mineral = 1;
b_mineral = 1;
c_mineral = 1;
v_pore = [0.9 0.1];
uf = 0;
a_pore = [1e-6,1e-6];
b_pore = [1e-6,1e-6];
c_pore = [1e-6,1e-9];
orientation_type_m = "simple_orientation";
theta_m = 0;
phi_m = 0;
psi_m = 0;
mu_theta_m = 0;
sigma_theta_m = 0;
mu_phi_m = 0;
sigma_phi_m = 0;
alpha_theta_m = 0;
beta_theta_m = 0;
alpha_phi_m = 0;
beta_phi_m = 0;
aM_theta_m = 0;
aM_phi_m = 0;
orientation_type_f = {"simple_orientation","chaotic"};
theta_f = [0,0];
phi_f = [0,0];
psi_f = [0,0];
mu_theta_f = [0,0];
sigma_theta_f = [0,0];
mu_phi_f = [0,0];
sigma_phi_f = [0,0];
alpha_theta_f = [0,0];
beta_theta_f = [0,0];
alpha_phi_f = [0,0];
beta_phi_f = [0,0];
aM_theta_f = [0,0];
aM_phi_f = [0,0];
communicative = 'yes';
[c_dry]=...
 t_matrix(iso,[],[],[],[],[],[],[],[],pt,fr,sat,...
 v_mineral,c_minerals,c_matrix,a_mineral,b_mineral,c_mineral,...
 v_pore,0,uf,a_pore,b_pore,c_pore,...
 orientation_type_m,...
 theta_m,phi_m,psi_m,mu_theta_m,sigma_theta_m,mu_phi_m,...
 sigma_phi_m,alpha_theta_m,beta_theta_m,alpha_phi_m,...
 beta_phi_m,aM_theta_m,aM_phi_m,...
 orientation_type_f,...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f,'no');
c_eff = c_dry;
x0 = eye(3,3)*9.8692326671601e-16;
perm = eye(3,3).*t_matrix_perm(pt,x0,orientation_type_f,v_pore,[a_pore',b_pore',c_pore'],...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f);
counter = 0;
sg = 0;
so = 0.6*rand();
sw = 1 - so;
 rho_fluid = so*rho+sw*rhw+sg*rhg;
 rho= rho_mineral*(1-pt)+rho_fluid*pt;
kf = 1/(sw/kw+so/ko+0/kg);
eta = eta_w*sw + eta_w*sg + eta_o*so;
tau = logspace(-10,-3,50);
omega = logspace(0,9,50);
N_omega = length(omega);
N_sg = length(so);
for i = 1:numel(tau)
 for j = 1:N_omega
 d = inf;
%  c_bk = LowFrequency_BrownKorringa(c_dry,c_matrix,kf,pt);
 while d>0.1
 c_old = c_eff;
 [c_eff]=...
 t_matrix(iso,c_eff,rho,theta,phi,omega(j),tau(i),eta,perm,pt,fr,sat,...
 v_mineral,c_minerals,c_matrix,a_mineral,b_mineral,c_mineral,...
 v_pore,kf,uf,a_pore,b_pore,c_pore,...
 orientation_type_m,...
 theta_m,phi_m,psi_m,mu_theta_m,sigma_theta_m,mu_phi_m,...
 sigma_phi_m,alpha_theta_m,beta_theta_m,alpha_phi_m,...
 beta_phi_m,aM_theta_m,aM_phi_m,...
 orientation_type_f,...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f,communicative);
 d = norm(c_old-c_eff,"fro");
 end
y = Christoffel(c2dto4d(c_eff),rho,theta,phi);
Vp3(i,j) = y(1,1);
Vs13(i,j) = y(2,1);
Vs23(i,j) = y(3,1);
 end
end
Qp=real(Vp3.^2)./imag(Vp3.^2);
alpha3 = -1000./Qp;
figure(11)
surf(tau,omega,real(Vp3'))
xlabel("\tau")
ylabel("\omega, Hz") 
zlabel("Vp, km/s")
set(gca,'xscale','log','yscale','log','FontName','Times New Roman','FontSize',14)
colormap jet
figure(12)
surf(tau,omega,real(Vs13'))
xlabel("\tau")
ylabel("\omega, Hz")
zlabel("Vs1, km/s")
set(gca,'xscale','log','yscale','log','FontName','Times New Roman','FontSize',14)
colormap jet
figure(13)
surf(tau,omega,real(Vs23'))
xlabel("\tau")
ylabel("\omega, Hz")
zlabel("Vs2, km/s")
set(gca,'xscale','log','yscale','log','FontName','Times New Roman','FontSize',14)
colormap jet
figure(14)
surf(tau,omega,alpha3')
xlabel("\tau")
ylabel("\omega, Hz")
zlabel("\alpha, km/s")
set(gca,'xscale','log','yscale','log','FontName','Times New Roman','FontSize',14)
colormap jet
%%%%sens on aspect ratios%%%%
k_matrix = 77e9;
u_matrix = 32e9;
kw = 2.25e9;
ko = 0.8e9;
kg = 0.01e9;
sw = 1;
so = 0;
sg = 0;
rho_mineral = 2.71e3;
rho = 0.8e3;
rhw = 1e3;
rhg = 1;
eta_w = 1e-3;
eta_o = 60e-3;
eta_g = 0.02e-3;
rho_fluid = so*rho+sw*rhw+sg*rhg;
theta = 0;
phi = 0;
tau = 1e-5;
pt = 0.005;
rho= rho_mineral*(1-pt)+rho_fluid*pt;
fr = 0;
v_mineral = 1;
c_minerals = elastictotensor(k_matrix,u_matrix);
c_matrix = elastictotensor(k_matrix,u_matrix);
a_mineral = 1;
b_mineral = 1;
c_mineral = 1;
v_pore = [0.9 0.1];
uf = 0;
a_pore = [1e-6,1e-6];
b_pore = [1e-6,1e-6];
c_pore = [1e-6,1e-9];
orientation_type_m = "simple_orientation";
theta_m = 0;
phi_m = 0;
psi_m = 0;
mu_theta_m = 0;
sigma_theta_m = 0;
mu_phi_m = 0;
sigma_phi_m = 0;
alpha_theta_m = 0;
beta_theta_m = 0;
alpha_phi_m = 0;
beta_phi_m = 0;
aM_theta_m = 0;
aM_phi_m = 0;
orientation_type_f = {"simple_orientation","chaotic"};
theta_f = [0,0];
phi_f = [0,0];
psi_f = [0,0];
mu_theta_f = [0,0];
sigma_theta_f = [0,0];
mu_phi_f = [0,0];
sigma_phi_f = [0,0];
alpha_theta_f = [0,0];
beta_theta_f = [0,0];
alpha_phi_f = [0,0];
beta_phi_f = [0,0];
aM_theta_f = [0,0];
aM_phi_f = [0,0];
communicative = 'yes';
[c_dry]=...
 t_matrix(iso,[],[],[],[],[],[],[],[],pt,fr,sat,...
 v_mineral,c_minerals,c_matrix,a_mineral,b_mineral,c_mineral,...
 v_pore,0,uf,a_pore,b_pore,c_pore,...
 orientation_type_m,...
 theta_m,phi_m,psi_m,mu_theta_m,sigma_theta_m,mu_phi_m,...
 sigma_phi_m,alpha_theta_m,beta_theta_m,alpha_phi_m,...
 beta_phi_m,aM_theta_m,aM_phi_m,...
 orientation_type_f,...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f,'no');
c_eff = c_dry;
d = inf;
x0 = eye(3,3)*9.8692326671601e-16;
perm = eye(3,3).*t_matrix_perm(pt,x0,orientation_type_f,v_pore,[a_pore',b_pore',c_pore'],...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f);
counter = 0;
a1 = logspace(-6,-3,50);
a2 = logspace(-6,-3,50);
a3 = logspace(-6,-3,10);
a_pore = [1e-6,1e-6];
b_pore = [1e-6,1e-6];
c_pore = [1e-6,1e-9];
omega = logspace(0,9,50);
sg = 0.6*rand(); %linspace(0,.6,50);
N_omega = length(omega);
N_sg = length(so);
 so = 0;
 sw = 1 - sg;
 rho_fluid = so*rho+sw*rhw+sg*rhg;
 rho= rho_mineral*(1-pt)+rho_fluid*pt;
 kf = 1/(sw/kw+0/ko+sg/kg);
 eta = eta_w*sw + eta_w*sg + eta_o*so;
for i = 1:numel(a3)
 for l = 1:numel(a1)
 for k = 1:numel(a2)
 a_pore = [1e-6,a1(l)];
 b_pore = [1e-6,a2(k)];
 c_pore = [1e-6,a3(i)];
 perm = eye(3,3).*t_matrix_perm(pt,x0,orientation_type_f,v_pore,[a_pore',b_pore',c_pore'], ...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f);
 for j = 1:N_omega
 d = inf;
%  c_bk = LowFrequency_BrownKorringa(c_dry,c_matrix,kf,pt);
 while d>0.1
 c_old = c_eff;
 [c_eff]=...
 t_matrix(iso,c_eff,rho,theta,phi,omega(j),tau,eta,perm,pt,fr,sat,...
 v_mineral,c_minerals,c_matrix,a_mineral,b_mineral,c_mineral,...
 v_pore,kf,uf,a_pore,b_pore,c_pore,...
 orientation_type_m,...
 theta_m,phi_m,psi_m,mu_theta_m,sigma_theta_m,mu_phi_m,...
 sigma_phi_m,alpha_theta_m,beta_theta_m,alpha_phi_m,...
 beta_phi_m,aM_theta_m,aM_phi_m,...
 orientation_type_f,...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f,communicative);
 d = norm(c_old-c_eff,"fro");
 end
 y = Christoffel(c2dto4d(c_eff),rho,theta,phi);
 Vp4(i,l,k,j) = y(1,1);
 Vs14(i,l,k,j) = y(2,1);
 Vs24(i,l,k,j) = y(3,1);
 
 end
 end
 end
end
Qp=real(Vp4.^2)./imag(Vp4.^2);
alpha4 = -1000./Qp;
% save all2


%%%%%%%%%%%%
%    NEW
%%%%%%%%%%%%
%%%%sens on aspect ratios/permeability%%%%
%%%%sens on aspect ratios%%%%
k_matrix = 77e9;
u_matrix = 32e9;
kw = 2.25e9;
ko = 0.8e9;
kg = 0.01e9;
sw = 1;
so = 0;
sg = 0;
rho_mineral = 2.71e3;
rho = 0.8e3;
rhw = 1e3;
rhg = 1;
eta_w = 1e-3;
eta_o = 60e-3;
eta_g = 0.02e-3;
rho_fluid = so*rho+sw*rhw+sg*rhg;
theta = 0;
phi = 0;
tau = 1e-5;
pt = 0.005;
rho= rho_mineral*(1-pt)+rho_fluid*pt;
fr = 0;
v_mineral = 1;
c_minerals = elastictotensor(k_matrix,u_matrix);
c_matrix = elastictotensor(k_matrix,u_matrix);
a_mineral = 1;
b_mineral = 1;
c_mineral = 1;
v_pore = [0.9 0.1];
uf = 0;
a_pore = [1e-6,1e-6];
b_pore = [1e-6,1e-6];
c_pore = [1e-6,1e-9];
orientation_type_m = "simple_orientation";
theta_m = 0;
phi_m = 0;
psi_m = 0;
mu_theta_m = 0;
sigma_theta_m = 0;
mu_phi_m = 0;
sigma_phi_m = 0;
alpha_theta_m = 0;
beta_theta_m = 0;
alpha_phi_m = 0;
beta_phi_m = 0;
aM_theta_m = 0;
aM_phi_m = 0;
orientation_type_f = {"simple_orientation","chaotic"};
theta_f = [0,0];
phi_f = [0,0];
psi_f = [0,0];
mu_theta_f = [0,0];
sigma_theta_f = [0,0];
mu_phi_f = [0,0];
sigma_phi_f = [0,0];
alpha_theta_f = [0,0];
beta_theta_f = [0,0];
alpha_phi_f = [0,0];
beta_phi_f = [0,0];
aM_theta_f = [0,0];
aM_phi_f = [0,0];
communicative = 'yes';
d = inf;
x0 = eye(3,3)*9.8692326671601e-16;
perm = eye(3,3).*t_matrix_perm(pt,x0,orientation_type_f,v_pore,[a_pore',b_pore',c_pore'],...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f);
counter = 0;
N1 = 50;%50
N2 = 50;%50
N3 = 10;%10
a1 = logspace(-6,-3,N1);
a2 = logspace(-6,-3,N2);
a3 = logspace(-6,-3,N3);
a_pore = [1e-6,1e-6];
b_pore = [1e-6,1e-6];
c_pore = [1e-6,1e-9];
omega = logspace(0,9,50);%(0,9,55)
sg = 0.6*rand(); %linspace(0,.6,50);
N_omega = length(omega);
N_sg = length(so);
 so = 0;
 sw = 1 - sg;
 rho_fluid = so*rho+sw*rhw+sg*rhg;
 rho= rho_mineral*(1-pt)+rho_fluid*pt;
 kf = 1/(sw/kw+0/ko+sg/kg);
 eta = eta_w*sw + eta_w*sg + eta_o*so;
     [c_dry ]=...
 t_matrix(iso,[],[],[],[],[],[],[],[],pt,fr,sat,...
 v_mineral,c_minerals,c_matrix,a_mineral,b_mineral,c_mineral,...
 v_pore,0,uf,a_pore,b_pore,c_pore,...
 orientation_type_m,...
 theta_m,phi_m,psi_m,mu_theta_m,sigma_theta_m,mu_phi_m,...
 sigma_phi_m,alpha_theta_m,beta_theta_m,alpha_phi_m,...
 beta_phi_m,aM_theta_m,aM_phi_m,...
 orientation_type_f,...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f,'no');
for i = 1:N3
c_eff = c_dry;
 for l = 1:N1
 for k = 1:N2
 a_pore = [1e-6,a1(l)];
 b_pore = [1e-6,a2(k)];
 c_pore = [1e-6,a3(i)];
 perm = eye(3,3).*t_matrix_perm(pt,x0,orientation_type_f,v_pore,[a_pore',b_pore',c_pore'], ...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f);
 permeability(i,l,k) = perm(1,1);
 
 for j = 1:N_omega
 d = inf;

 while d>0.1
 c_old = c_eff;
 [c_eff]=...
 t_matrix(iso,c_eff,rho,theta,phi,omega(j),tau,eta,perm,pt,fr,sat,...
 v_mineral,c_minerals,c_matrix,a_mineral,b_mineral,c_mineral,...
 v_pore,kf,uf,a_pore,b_pore,c_pore,...
 orientation_type_m,...
 theta_m,phi_m,psi_m,mu_theta_m,sigma_theta_m,mu_phi_m,...
 sigma_phi_m,alpha_theta_m,beta_theta_m,alpha_phi_m,...
 beta_phi_m,aM_theta_m,aM_phi_m,...
 orientation_type_f,...
 theta_f,phi_f,psi_f,mu_theta_f,sigma_theta_f,mu_phi_f,...
 sigma_phi_f,alpha_theta_f,beta_theta_f,alpha_phi_f,...
 beta_phi_f,aM_theta_f,aM_phi_f,communicative);
 d = norm(c_old-c_eff,"fro");
 end
 y = Christoffel(c2dto4d(c_eff),rho,theta,phi);
 Vp4(i,l,k,j) = y(1,1);
 Vs14(i,l,k,j) = y(2,1);
 Vs24(i,l,k,j) = y(3,1);
 
 end
 end
 end
end
Qp=real(Vp4.^2)./imag(Vp4.^2);
alpha4 = -1000./Qp;
% save all3
