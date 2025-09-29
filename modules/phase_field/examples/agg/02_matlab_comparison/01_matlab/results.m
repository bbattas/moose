% Compiling the gr0/gr1 data to csv for reading into python to compare with
% PF results with plotting

outfolder = 'results';
if ~exist(outfolder, 'dir')
    mkdir(outfolder);
end

% 5s results
matlab_to_csv('m2/incl_dependence_bicr_m2nn50000.mat', 'results/matlab_m2_05s.csv')
matlab_to_csv('m4/incl_dependence_bicr_m4nn50000.mat', 'results/matlab_m4_05s.csv')
matlab_to_csv('m6/incl_dependence_bicr_m6nn50000.mat', 'results/matlab_m6_05s.csv')

% 10 results
matlab_to_csv('m2/incl_dependence_bicr_m2nn100000.mat', 'results/matlab_m2_10s.csv')
matlab_to_csv('m4/incl_dependence_bicr_m4nn100000.mat', 'results/matlab_m4_10s.csv')
matlab_to_csv('m6/incl_dependence_bicr_m6nn100000.mat', 'results/matlab_m6_10s.csv')

% 15s results
matlab_to_csv('m2/incl_dependence_bicr_m2nn150000.mat', 'results/matlab_m2_15s.csv')
matlab_to_csv('m4/incl_dependence_bicr_m4nn150000.mat', 'results/matlab_m4_15s.csv')
matlab_to_csv('m6/incl_dependence_bicr_m6nn150000.mat', 'results/matlab_m6_15s.csv')

% 20s results
matlab_to_csv('m2/incl_dependence_bicr_m2nn200000.mat', 'results/matlab_m2_20s.csv')
matlab_to_csv('m4/incl_dependence_bicr_m4nn200000.mat', 'results/matlab_m4_20s.csv')
matlab_to_csv('m6/incl_dependence_bicr_m6nn200000.mat', 'results/matlab_m6_20s.csv')

% 25s results
matlab_to_csv('m2/incl_dependence_bicr_m2nn250000.mat', 'results/matlab_m2_25s.csv')
matlab_to_csv('m4/incl_dependence_bicr_m4nn250000.mat', 'results/matlab_m4_25s.csv')
matlab_to_csv('m6/incl_dependence_bicr_m6nn250000.mat', 'results/matlab_m6_25s.csv')

% 30 results
matlab_to_csv('m2/incl_dependence_bicr_m2nn300000.mat', 'results/matlab_m2_30s.csv')
matlab_to_csv('m4/incl_dependence_bicr_m4nn300000.mat', 'results/matlab_m4_30s.csv')
matlab_to_csv('m6/incl_dependence_bicr_m6nn300000.mat', 'results/matlab_m6_30s.csv')

% 40s results
matlab_to_csv('m2/incl_dependence_bicr_m2einde400000.mat', 'results/matlab_m2_40s.csv')
matlab_to_csv('m4/incl_dependence_bicr_m4einde400000.mat', 'results/matlab_m4_40s.csv')
matlab_to_csv('m6/incl_dependence_bicr_m6einde400000.mat', 'results/matlab_m6_40s.csv')
