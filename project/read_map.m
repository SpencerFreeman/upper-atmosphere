function zbar = read_map(xbar, map)

llabar = ecef2lla(xbar(1:3, :)')';

zbar = interp2(...
    map(:, :, 1), map(:, :, 2), map(:, :, 3), ...
    360 + llabar(2, :), llabar(1, :));

end