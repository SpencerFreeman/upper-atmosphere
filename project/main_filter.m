clear;clc;close all

%% load data + plot reference map
addpath('functions')
addpath('data')
load('EMAG2_V3_Blacksburg-Roanoke', 'data')

% text = fileread('data/FRD20240418.json');
text = fileread('data/FRD20220203.json'); % February Solar Storm (killed Starlink satellites)
temporal_data = jsondecode(text);

n_temporal_data = length(temporal_data.datetime);
t_temporal_data = NaT(n_temporal_data, 1, 'TimeZone', 'Etc/UTC');
for i = 1:n_temporal_data
    t_temporal_data(i) = datetime( ...
        temporal_data.datetime{i}(1:(end - 5)), ...
        'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss', ...
        'TimeZone', 'Etc/UTC');
end
t_temporal_data.TimeZone = "America/New_York";

[XYZ,H,D,I,F] = igrfmagm( ...
    temporal_data.x_info.altitude, ...
    temporal_data.x_info.latitude, temporal_data.x_info.longitude, ...
    decyear( ...
        t_temporal_data(1).Year, ...
        t_temporal_data(1).Month, ...
        t_temporal_data(1).Day));

h0 = figure;
h0.WindowStyle = 'Docked';
plot(t_temporal_data, temporal_data.S - F)
grid on
ylabel('Magnetic Field Strength (nT)')
xlabel('Time (s)')
text()

lat_range = [ 36.758719  37.556273];
lon_range = [-80.874399 -79.464463];

[map, h] = plot_mag_countour(data, lat_range, lon_range);

%% generate truth + measurements
lla0 = [36.90165855141354, -79.70578451436336, 4e3]; % deg, deg
% lla0 = [37.1979516376812, -79.5787821351352, 4e3]; % deg, deg
% lla0 = [37.0000, -79.9198, 4e3]; % deg, deg
llaf = [37.22891412982679, -80.43067879145124, 4e3]; % deg, deg

Ts = 1;%1/10; % s
[xs_truth, ts_truth, n] = generate_trajectory(lla0, llaf, 50, Ts); %

zs_truth = read_map(xs_truth, map); % using truth measurements, perfect sensor

%% filter
use_mag = true;

R = 2e1;%1e-5; % measurement noise, very small for perfect sensor, nT
R2 = 1e1;

qtilda = 1e0; % m^2/s^5
Q = process_noise(qtilda, Ts);
I3 = eye(3);
z3 = zeros(3);
Phi = [...
    I3,  I3*Ts, .5*I3*Ts^2; ...
    z3, I3,    I3*Ts; ...
    z3, z3,   I3];

% x0 = [xs_truth(1:3, 1); zeros(6, 1)]; % m, m/s, m/s^2
x0 = xs_truth(:, 1) + ...
    [randn(3, 1) * 10; randn(3, 1) * 2; randn(3, 1) * 0]; % initialize with truth

% x0 = [913107.851200612; -5027503.00541662; 3811061.22541831; -39.8270679230392; 8.35545013041825; 22.4064255614077; -5.61802231030889e-05; 0.000309317309728480; -0.000234476519221953];

init_pos_err = 1e3; % m
large_v = 5; % m/s
large_a = .05; % m/s^2
p0 = zeros(9); % covariance matrix
for i = 1:3
    i3 = i + 3;
    i6 = i + 6;
    p0(i, i) = init_pos_err^2;
    p0(i3, i3) = large_v^2;
    p0(i6, i6) = large_a^2;
end

% create map reading function handle
f = @(x) read_map(x, map);
fa = @(x) norm(x(1:3));

% loop through measurements -----------------------------
tic
xhat = x0;
phat = p0;
NIS = 0;
zbar = 0;
xs = nan(length(x0), n);
NISs = nan(1, n);
for i = 1:n % %%%% loop measurements %%%%%

    if isnan(zbar); i = i - 1; break; end % break if off the map

    xs(:, i) = xhat;
    NISs(i) = NIS;

    % propagate state and covariance
    M = Phi*phat*Phi' + Q; % covariance before update
    xbar = Phi*xhat; % estimate before update

    % Update filter state and covariance matrix using the measurement
    if use_mag
        z = zs_truth(i) + randn * .1; % current magnetometer reading, nT

        zbar = read_map(xbar, map); % consult the map! nT

        resid = z - zbar; % nT

        H = numerical_jacobian(f, 1, xbar); % state-measurement linearization, nT/m

        Sj = H*M*H' + R;
        K = M*H'*inv(Sj); % calculate optimal gain
        IKH = eye(9) - K*H;
        phat = IKH*M*IKH' + K*R*K'; % update covariance (Joseph form)

        xhat = xbar + K*resid; % update estimate

        S = R + H*phat*H'; % residual covariance
        NIS = resid'*(S\resid); % normalized innovation (residual) squared
    else
        xhat = xbar;
        phat = M;
    end

    % process altimeter measurement
    [xhat, phat] = update_altimeter(xhat, phat, xs_truth, i, fa, R2);

end %%%% loop measurements %%%%%
toc

lla_estimate = ecef2lla(xs(1:3, :)')';

r_error = vecnorm(xs_truth(1:3, 1:i) - xs(1:3, :), 2, 1);
v_error = vecnorm(xs_truth(4:6, 1:i) - xs(4:6, :), 2, 1);
a_error = vecnorm(xs_truth(7:9, 1:i) - xs(7:9, :), 2, 1);

%% plotting

fprintf('Mean Error Stats ---------------\n')
fprintf('\tPosition:     %f (m)\n', mean(r_error))
fprintf('\tVelocity:     %f (m/s)\n', mean(v_error))
fprintf('\tAcceleration: %f (m/s^2)\n', mean(a_error))

h2 = figure;
h2.WindowStyle = 'Docked';
plot(ts_truth(1:i), NISs)
grid on

figure(h)
plot(llaf(2), llaf(1), 'xk')
plot(linspace(lla0(2), llaf(2), n), linspace(lla0(1), llaf(1), n), '--k')
plot(lla_estimate(2, :), lla_estimate(1, :), 'm')










