% Compiling the gr0/gr1 data to csv for reading into python to compare with
% PF results with plotting

outfolder = 'results';
if ~exist(outfolder, 'dir')
    mkdir(outfolder);
end

% 5s results
matlab_to_csv('m2/incl_dependence_bicr_m2nn5000.mat', 'results/matlab_m2_05s.csv')
matlab_to_csv('m4/incl_dependence_bicr_m4nn5000.mat', 'results/matlab_m4_05s.csv')
matlab_to_csv('m6/incl_dependence_bicr_m6nn5000.mat', 'results/matlab_m6_05s.csv')

% 10 results
matlab_to_csv('m2/incl_dependence_bicr_m2nn10000.mat', 'results/matlab_m2_10s.csv')
matlab_to_csv('m4/incl_dependence_bicr_m4nn10000.mat', 'results/matlab_m4_10s.csv')
matlab_to_csv('m6/incl_dependence_bicr_m6nn10000.mat', 'results/matlab_m6_10s.csv')

% 15s results
matlab_to_csv('m2/incl_dependence_bicr_m2nn15000.mat', 'results/matlab_m2_15s.csv')
matlab_to_csv('m4/incl_dependence_bicr_m4nn15000.mat', 'results/matlab_m4_15s.csv')
matlab_to_csv('m6/incl_dependence_bicr_m6nn15000.mat', 'results/matlab_m6_15s.csv')

% 20s results
matlab_to_csv('m2/incl_dependence_bicr_m2einde20000.mat', 'results/matlab_m2_20s.csv')
matlab_to_csv('m4/incl_dependence_bicr_m4einde20000.mat', 'results/matlab_m4_20s.csv')
matlab_to_csv('m6/incl_dependence_bicr_m6einde20000.mat', 'results/matlab_m6_20s.csv')