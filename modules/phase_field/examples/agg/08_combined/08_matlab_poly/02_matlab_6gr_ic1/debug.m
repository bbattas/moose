load('incl_dependence_6greinde60000')
whos

% Microstructure
imagesc(microstructure')
colorbar
colormap('viridis')
axis equal tight xy

% All 6 order parameters
figure;
for i = 1:6
    subplot(2, 3, i);
    imagesc(eta{i}');      % <-- transpose here
    colorbar;
    title(['eta{' num2str(i) '}']);
    axis equal tight xy;
end


% imagesc(sumetasqu)
% colorbar
% colormap('jet')

% surf(eta{1})


% for i = 1:6
%   printf('eta{%d}: size = %s, class = %s\n', i, mat2str(size(eta{i})), class(eta{i}));
% end

% eta{1}
% eta{1}(1:5, 1:5)
