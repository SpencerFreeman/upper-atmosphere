clear;clc;close all

%% load data + plot reference map
load('EMAG2_V3_Blacksburg-Roanoke', 'data')

lat_range = [ 36.758719  37.556273];
lon_range = [-80.874399 -79.464463];

[Xq, Yq, Vq] = plot_mag_countour(data, lat_range, lon_range);

%% generate truth + measurements
lla0 = [36.90165855141354, -79.70578451436336, 4e3]; % deg, deg
llaf = [37.22891412982679, -80.43067879145124, 4e3]; % deg, deg

plot(llaf(2), llaf(1), 'xk')

n = 1000;
[xs_truth, lats, lons, alts] = generate_trajectory(lla0, llaf, 10, n); % m\

plot(linspace(lla0(2), llaf(2), n), linspace(lla0(1), llaf(1), n), '--k')

zs_truth = griddata(Xq(:), Yq(:), Vq(:), 360 + lons, lats);


%% filter

Ts = 1/10; % s

R = 1e-5; % measurement noise, very small for perfect sensor, T

qtilda = 1e0; % m^2/s^5

Q = process_noise(qtilda, Ts);

I3 = eye(3);
z3 = zeros(3);
Phi = [...
    I,  I*Ts, .5*I*Ts^2; ...
    z3, I,    I*Ts; ...
    z3, z3,   I];


x0_pos = lla2ecef([lla0, 4e3]); % m

x0 = [x0_pos(:); zeros(6, 1)]; % m, m/s, m/s^2

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

% loop through measurements -----------------------------
xhat = x0;
phat = p0;
for i = 1:n % %%%% loop measurements %%%%%
    % propagate state and covariance
    M = Phi*phat*Phi' + Q; % covariance before update
    xbar = Phi*xhat; % estimate before update

    % Update filter state and covariance matrix using the measurement
    z = zs(i);

    resid = z - H*xbar;

    Sj = H*M*H' + R;
    K = M*H'*inv(Sj); % calculate optimal gain
    IKH = eye(9) - K*H;
    phat = IKH*M*IKH' + K*R*K'; % update covariance (Joseph form)

    xhat = xbar + K*resid; % update estimate

    S = R + H*phat*H'; % residual covariance
    NIS = resid'*(S\resid); % normalized innovation (residual) squared
end %%%% loop measurements %%%%%










