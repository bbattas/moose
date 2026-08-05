function [kappa,g,L,m,deltax,deltat,order_par,bound_cond,delta_s,delta_s_kin,n_fold,n_fold_kin,ofset_angle,ofset_angle_kin,anis] = incl_depence_input

kappa = 0.3; %% see parameters sim 5
g = [0 sqrt(2)/3 sqrt(2)/3 sqrt(2)/3 sqrt(2)/3 sqrt(2)/3;  0 0 sqrt(2)/3 sqrt(2)/3 sqrt(2)/3 sqrt(2)/3; 0 0 0 sqrt(2)/3 sqrt(2)/3 sqrt(2)/3; ...
    0 0 0 0 sqrt(2)/3 sqrt(2)/3;  0 0 0 0 0 sqrt(2)/3;  0 0 0 0 0 0];
% gamma_l = 0.6328; gamma_h = 5.1364;
%gamma = [0 gamma_h gamma_l gamma_h gamma_l gamma_h;  0 0 gamma_h gamma_l gamma_h gamma_l; 0 0 0 gamma_h gamma_l gamma_h; ...
 %   0 0 0 0 gamma_h gamma_l;  0 0 0 0 0 gamma_h;  0 0 0 0 0 0];
%L_l = 0.4167; L_h = 1.1111;
L_l=0.833; L_h=0.833;
L = [0 L_h L_l L_h L_l L_h; 0 0 L_h L_l L_h L_l; ...
    0 0 0 L_h L_l L_h; 0 0 0 0 L_h L_l; 0 0 0 0 0 L_h; 0 0 0 0 0 0];
m = 0.9375;
deltax = 0.1;
deltat = 0.001;
order_par = 6;
bound_cond = 'pp';
delta_s = [0 0.05 0.05 0.05 0.05 0.05;0 0 0.05 0.05 0.05 0.05; 0 0 0 0.05 0.05 0.05; ...
    0 0 0 0 0.05 0.05; 0 0 0 0 0 0.05 ; 0 0  0 0 0 0];
% delta_s = zeros(6,6);

%delta_s_kin = 4*[0 0.15 0.15 0.15 0.15 0.15;0 0 0.15 0.15 0.15 0.15; 0 0 0 0.15 0.15 0.15; ...
 %   0 0 0 0 0.15 0.15; 0 0 0 0 0 0.15 ; 0 0  0 0 0 0];
%delta_s = 1*[0 0 0 0.15 0 0; 0 0 0 0 0.15 0; 0 0 0 0 0 0.15; ...
 %   0 0 0 0 0 0; 0 0 0 0 0 0 ; 0 0  0 0 0 0];
delta_s_kin = zeros(6,6);
%delta_s = [0 0;0 0];
n_fold = [0 2 2 2 2 2; 0 0 2 2 2 2; 0 0 0 2 2 2; 0 0 0 0 2 2 ; 0 0 0 0 0 2; 0 0 0 0 0 0];
n_fold_kin = [0 2 2 2 2 2; 0 0 2 2 2 2; 0 0 0 2 2 2; 0 0 0 0 2 2 ; 0 0 0 0 0 2; 0 0 0 0 0 0];
ofset_angle = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
%ofset_angle = [0 0 0 0 4*pi/6 0; 0 0 0 pi/6 0 5*pi/6; 0 0 0 0 2*pi/6 0; 0 0 0 0 0 3*pi/6; 0 0 0 0 0 0; 0 0 0 0 0 0];
%ofset_angle = zeros(6,6);
ofset_angle_kin = zeros(6,6);
anis = [0 'w' 'w' 'w' 'w' 'w'; 0 0 'w' 'w' 'w' 'w'; 0 0 0 'w' 'w' 'w'; 0 0 0 0 'w' 'w'; 0 0 0 0 0 'w'; 0 0 0 0 0 0];
