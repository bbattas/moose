function [kappa,g,L,m,deltax,deltat,order_par,bound_cond,delta_s,delta_s_kin,n_fold,n_fold_kin,ofset_angle,ofset_angle_kin,anis] = incl_depence_input

kappa = 2.100e07; %0.3; % energy gradient coefficient
g = [0 0.5097; 0 0]; %0.6796 is 1x % g-values
L = [0 2.0861e-06; 0 0]; %[0 0.8333; 0 0];-6 % kinetic coefficient in Allen Cahn equations
m = 1.305e7; %0.9375; % prefactor in homogeneous free energy
deltax = 0.5; %1; %1 %0.2; % grid spacing
deltat = 0.0002; %0.003; % time step
order_par = 2; % number of order parameters (eta) or number different grain types considered
bound_cond = 'pp'; % type of boundary conditions; pp : periodic (possible vlaues 'pp' (periodic), 'pn' (periodic along 1 direction and Neumann along the other), 'nn' (Neumann), 'dd' (Dirichlett))
delta_s = [0 0.0; 0 0]; % anisotropy strength factor grain boundary energy

delta_s_kin = [0 0; 0 0]; % anisotropy strength factor grain boundary mobility
n_fold = [0 2; 0 0]; % 4- or 6- fold symmetry in grainboundary energy
n_fold_kin = [0 2; 0 0]; % 4- or 6-fold symmetry in grain boundary mobility
ofset_angle_kin = [0 0; 0 0]; % offset angle for inclination dependence mobility
ofset_angle = [0 0; 0 0]; % offset  angle for inclination dependence grain boundary energy
anis = ['w' 'w'; 'w' 'w']; % type of anisotropy function for grain boundary energy (possible values 'w', 'c', 's')
%'w' : %f_sigma = (1+delta_sigma*cos(n*(theta+theta_ofset)))
%'c' : %cubes scaled : f_sigma = (1/1+delta_sigma)*(1+delta_sigma(abs(cos(theta+theta_ofset))+abs(sin(theta+theta_ofset))))
%'s' : cubes not scaled : f_sigma = (1+delta_sigma(abs(cos(theta+theta_ofset))+abs(sin(theta+theta_ofset))))
