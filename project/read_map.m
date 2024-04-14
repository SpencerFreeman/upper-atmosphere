function zbar = read_map(xbar, map)

llabar = ecef2lla(xbar(1:3, :)')';

% tic
zbar = interp2(...
    map(:, :, 1), map(:, :, 2), map(:, :, 3), ...
    360 + llabar(2, :), llabar(1, :));
% toc

% X = map(:, 1, 1);
% Y = map(1, :, 2);
% V = map(:, :, 3);
% 
% tic 
% for i = 1:size(zbar, 2)
%     zbar(i) = lininterp2(X, Y, V, 360 + llabar(2, i), llabar(1, i));
% end
% toc

end