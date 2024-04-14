function [x_truth, t_truth, n] = generate_trajectory(lla0, llaf, speed, Ts)

G = 6.6743e-11;
m_earth = 5.9722e24; % kg;

groundrange = distance(lla0(1), lla0(2), llaf(1), llaf(2), wgs84Ellipsoid("m"));

tend = groundrange / speed; % s

n = round(tend/Ts); % number of samples

lats = linspace(lla0(1), llaf(1), n);
lons = linspace(lla0(2), llaf(2), n);
alts = linspace(lla0(3), llaf(3), n);

r_truth = lla2ecef([lats(:), lons(:), alts(:)])';

v_truth = diff(r_truth')' / Ts; % m/s
v_truth = [v_truth(:, 1), v_truth];

rmag = vecnorm(r_truth); % m
rhat = r_truth ./ rmag;

% a_truth = G*m_earth ./ rmag.^2 .* -rhat; % m/s^2

a_truth = -rhat * speed^2 ./ rmag; % m/s^2

x_truth = [r_truth; v_truth; a_truth];

t_truth = Ts * (0:(n - 1)); % s

end