function start_incl_dependence
addpath .
% invoer = 'incl_dependence_single_grain_inputt';%'incl_dependence_bicr';
invoer = 'incl_dependence_bicr_m4';
n = 20000;%30000; %30000;
n_old = 0;

%   eta{1} = zeros(100);
% for ii = 1 : 100
%     for jj = 1 : 100
%         if (50-ii)^2 + (50-jj)^2 < 40^2 %%changed for the simulations with a
%             eta{1}(ii,jj) = 1;
%         end
%     end
% end

% --- domain & initial condition (replace the nested loops) ---
N  = 161;          % 161 points -> 160 intervals ("elements")
cx = 81; cy = 81;  % center index (1-based). 81 gives symmetric -80..80
% cx = 321; cy = 321;  % center index (1-based). 81 gives symmetric -80..80
R  = 60;           % radius in grid points (physical R = R * deltax)

eta{1} = zeros(N, N);
[I, J] = ndgrid(1:N, 1:N);

mask = ( (I - cx).^2 + (J - cy).^2 ) < R^2;
eta{1}(mask) = 1;
eta{2} = 1 - eta{1};


% eta{2} = 1- eta{1};
grain_area = [];
%load incl_dependence_single_grain_input_kin_acnn20000
% Was saving at 2500 for 30000

[eta,sumetasqu,grain_area] = anisotropic_grgr_incl_in_gamma(invoer,n,n_old,eta,1000,grain_area);
