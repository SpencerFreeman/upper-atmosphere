function [x_truth, lats, lons, alts] = generate_trajectory(lla0, llaf, v, n)

lats = linspace(lla0(1), llaf(1), n);
lons = linspace(lla0(2), llaf(2), n);
alts = linspace(lla0(3), llaf(3), n);

x_truth = lla2ecef([lats(:), lons(:), alts(:)])';


end