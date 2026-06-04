function plot_frames()

% Get all nn files
files = dir('incl_dependence_4gr_hexnn*');

% Extract trailing step numbers and sort numerically
steps = arrayfun(@(f) str2double(regexp(f.name, '\d+$', 'match', 'once')), files);
[steps, idx] = sort(steps);
files = files(idx);

for i = 1:length(files)
    % Load the file
    load(files(i).name);

    % Compute frame number
    frame = steps(i) / 2500;

    % Plot
    figure('Visible', 'off');
    imagesc(microstructure')
    colorbar
    colormap('viridis')
    axis equal tight xy
    title(sprintf('Step %d (t = %d)', frame, steps(i)))

    % Save
    saveas(gcf, sprintf('frame_%03d.png', frame));
    close(gcf);

    fprintf('Saved frame %d / %d (step %d)\n', i, length(files), steps(i));
end

end
