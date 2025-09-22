function matlab_to_csv(matfile, outfile)
% export_eta_to_csv.m
% Usage:
%   matlab_to_csv('incl_dependence_bicr_m2nn1000.mat', 'bicrystal_eta.csv')
    if nargin < 2
        error('Usage: export_eta_to_csv(<matfile>, <outfile.csv>)');
    end

    % Load .mat
    S = load(matfile);

    % Expect eta as a cell array: eta{1}, eta{2}
    if ~isfield(S, 'eta') || ~iscell(S.eta) || numel(S.eta) < 2
        error('File does not contain eta{1} and eta{2}.');
    end
    eta1 = S.eta{1};
    eta2 = S.eta{2};

    % Grid info (from your description)
    N = 161;
    deltax = 1.0;

    % Build coordinates 0..160
    x = (0:N-1) * deltax;
    y = (0:N-1) * deltax;

    % Use ndgrid to match eta?s row/col layout (size N x N)
    [X, Y] = ndgrid(x, y);

    % Flatten to columns (MATLAB is column-major so (: ) is consistent)
    Xc    = X(:);
    Yc    = Y(:);
    eta1c = eta1(:);
    eta2c = eta2(:);

    % Make a table and write CSV
    T = table(Xc, Yc, eta1c, eta2c, 'VariableNames', {'x','y','eta1','eta2'});
    writetable(T, outfile);

    fprintf('Wrote %d rows to %s\n', height(T), outfile);
end
