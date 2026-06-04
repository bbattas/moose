% =========================================================
% extract_contours.m
% Extracts contours of sumetasqu = 0.75 from phase field
% .mat files and saves results to JSON for use in Python.
% =========================================================

%% 1. CONFIGURATION
data_dir       = '.';
file_pattern   = 'incl_dependence_6grnn*';
contour_level  = 0.75;
dx             = 0.2;
output_file    = 'contours_matlab_6gr.json';
debug_fig_file = 'debug_contour_step1.png';

%% 2. DISCOVER AND SORT FILES BY STEP NUMBER
listing = dir(fullfile(data_dir, file_pattern));

if isempty(listing)
    error('No files found matching pattern: %s', file_pattern);
end

n_files    = numel(listing);
step_nums  = zeros(n_files, 1);
for i = 1:n_files
    tok = regexp(listing(i).name, 'nn(\d+)$', 'tokens');
    if isempty(tok)
        error('Could not parse step number from filename: %s', listing(i).name);
    end
    step_nums(i) = str2double(tok{1}{1});
end

[step_nums_sorted, sort_idx] = sort(step_nums);
listing_sorted = listing(sort_idx);

fprintf('Found %d files. Step range: %d to %d\n', ...
    n_files, step_nums_sorted(1), step_nums_sorted(end));

%% 3. LOOP OVER FILES AND EXTRACT CONTOURS
all_contours = cell(n_files, 1);

for i = 1:n_files
    filepath   = fullfile(data_dir, listing_sorted(i).name);
    step_num   = step_nums_sorted(i);

    fprintf('Processing step %d: %s\n', step_num, listing_sorted(i).name);

    data       = load(filepath);
    sumetasqu  = data.sumetasqu;

    nx = size(sumetasqu, 1);
    ny = size(sumetasqu, 2);
    x  = (0:nx-1) * dx;
    y  = (0:ny-1) * dx;

    C = contourc(x, y, sumetasqu', [contour_level, contour_level]);

    segments = {};
    k = 1;
    while k < size(C, 2)
        num_pts = C(2, k);
        seg_x   = C(1, k+1 : k+num_pts);
        seg_y   = C(2, k+1 : k+num_pts);
        segments{end+1} = [seg_x(:), seg_y(:)];
        k = k + num_pts + 1;
    end

    all_contours{i} = struct('step_number', step_num, 'segments', {segments});
end

%% 4. WRITE JSON MANUALLY (Octave-compatible)
fid = fopen(output_file, 'w');
if fid == -1
    error('Could not open output file for writing: %s', output_file);
end

fprintf(fid, '[\n');
for i = 1:n_files
    step_num = all_contours{i}.step_number;
    segs     = all_contours{i}.segments;
    n_segs   = numel(segs);

    fprintf(fid, '  {\n');
    fprintf(fid, '    "step_number": %d,\n', step_num);
    fprintf(fid, '    "segments": [\n');

    for s = 1:n_segs
        seg   = segs{s};
        n_pts = size(seg, 1);

        fprintf(fid, '      [\n');
        for p = 1:n_pts
            if p < n_pts
                fprintf(fid, '        [%.6f, %.6f],\n', seg(p,1), seg(p,2));
            else
                fprintf(fid, '        [%.6f, %.6f]\n',  seg(p,1), seg(p,2));
            end
        end

        if s < n_segs
            fprintf(fid, '      ],\n');
        else
            fprintf(fid, '      ]\n');
        end
    end

    fprintf(fid, '    ]\n');
    if i < n_files
        fprintf(fid, '  },\n');
    else
        fprintf(fid, '  }\n');
    end
end
fprintf(fid, ']\n');
fclose(fid);

fprintf('Contours written to %s\n', output_file);

%% 5. DEBUG PLOT (first file only)
first_filepath  = fullfile(data_dir, listing_sorted(1).name);
first_step      = step_nums_sorted(1);
first_data      = load(first_filepath);
first_sumetasqu = first_data.sumetasqu;

nx_d = size(first_sumetasqu, 1);
ny_d = size(first_sumetasqu, 2);
x_d  = (0:nx_d-1) * dx;
y_d  = (0:ny_d-1) * dx;

fig = figure('Visible', 'off');
hold on;

first_segs = all_contours{1}.segments;
for s = 1:numel(first_segs)
    plot(first_segs{s}(:,1), first_segs{s}(:,2), 'b-', 'LineWidth', 1.5);
end

xlim([x_d(1), x_d(end)]);
ylim([y_d(1), y_d(end)]);
axis equal tight
xlabel('x');
ylabel('y');
title(sprintf('sumetasqu = %.2f contour(s), step %d', contour_level, first_step));

saveas(fig, debug_fig_file);
close(fig);

fprintf('Debug plot saved to %s\n', debug_fig_file);
