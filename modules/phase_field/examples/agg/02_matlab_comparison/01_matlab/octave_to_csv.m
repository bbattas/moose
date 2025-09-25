function octave_to_csv(matfile, outfile)
  % matlab_to_csv_nontable.m
  % Usage:
  %   matlab_to_csv_nontable('incl_dependence_bicr_m2nn1000.mat', 'bicrystal_eta.csv')

      if nargin < 2
          error('Usage: matlab_to_csv_nontable(<matfile>, <outfile.csv>)');
      end

      % Load .mat
      S = load(matfile);

      % Expect eta as a cell array: eta{1}, eta{2}
      if ~isfield(S, 'eta') || ~iscell(S.eta) || numel(S.eta) < 2
          error('File does not contain eta{1} and eta{2}.');
      end
      eta1 = S.eta{1};
      eta2 = S.eta{2};

      % Grid info (assuming 161x161 as before)
      N = size(eta1, 1);  % safer than hardcoding 161
      deltax = 0.2;%1.0;

      % Build coordinate arrays
      x = (0:N-1) * deltax;
      y = (0:N-1) * deltax;
      [X, Y] = ndgrid(x, y);

      % Flatten and concatenate into numeric matrix
      data = [X(:), Y(:), eta1(:), eta2(:)];

      % Write header first
      fid = fopen(outfile, 'w');
      if fid == -1
          error('Could not open %s for writing.', outfile);
      end
      fprintf(fid, 'x,y,eta1,eta2\n');   % header line

      % Write data rows (comma separated)
      fmt = '%.6f,%.6f,%.6f,%.6f\n';
      fprintf(fid, fmt, data.');
      fclose(fid);

      fprintf('Wrote %d rows to %s\n', size(data, 1), outfile);
  end
