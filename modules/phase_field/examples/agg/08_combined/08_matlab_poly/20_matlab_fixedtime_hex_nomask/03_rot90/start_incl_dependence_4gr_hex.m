function start_incl_dependence
addpath .
addpath ../../00_sub
invoer = 'incl_dependence_4gr_hex';
n = 5000;
n_old = 0;

% load_moose_ic.m
% Load MOOSE Exodus initial condition and populate eta cell array.
% Usage: run this before calling anisotropic_grgr_incl_in_gamma(...)

ic_file = '18_pf_fixedtime_hex_nomask_1.5708rot_ic_4gr_t0.mat';   % <-- set to your generated filename
data = load(ic_file);

fprintf('Loaded IC from: %s\n', data.source_file);
fprintf('Exodus time: %g,  step: %d\n', data.exodus_time, data.exodus_step);
fprintf('Grid size: %d x %d\n', data.nx, data.ny);

eta = cell(4, 1);
eta{1} = data.gr0;
eta{2} = data.gr1;
eta{3} = data.gr2;
eta{4} = data.gr3;
% eta{5} = data.gr4;
% eta{6} = data.gr5;

% Sanity check
for ii = 1:4
    fprintf('eta{%d}: size=[%d %d], min=%.4f, max=%.4f\n', ...
        ii, size(eta{ii},1), size(eta{ii},2), min(min(eta{ii})), max(max(eta{ii})));
end


grain_area = [];

% RUN
[eta,sumetasqu,grain_area] = anisotropic_grgr_incl_in_gamma(invoer,n,n_old,eta,50,grain_area);
