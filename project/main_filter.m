clear;clc;close all

%% load data + build map
addpath('functions')
addpath('data')
load('data/EMAG2_V3_Blacksburg-Roanoke', 'data')
load('data/EMAG2_V3_NGS-Fredericksburg', 'data_frd')
load('test.mat', 'f_temp', 't_temporal_frd_14day_avg')

% generate map via interpolation ------------------------------------------
lat_range = [ 36.758719  37.556273]; nlat = 70;
lon_range = [-80.874399 -79.464463]; nlon = 75;

map = build_map(data, lat_range, nlat, lon_range, nlon);

%% load temporal mag data --------------------------------------------------
% accessed from web: https://imag-data.bgs.ac.uk/GIN_V1/GINForms2
% text = fileread('data/FRD20240418.json');
% text = fileread('data/FRD20220202.json'); % February 2, 2022
text = fileread('data/FRD20220203.json'); % February 3, 2022 Solar Storm (killed Starlink satellites)
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

[XYZ,H,D,I,H_core_frd] = igrfmagm( ...
    temporal_data.x_info.altitude, ...
    temporal_data.x_info.latitude, temporal_data.x_info.longitude, ...
    decyear( ...
        t_temporal_data(1).Year, ...
        t_temporal_data(1).Month, ...
        t_temporal_data(1).Day));

H_crust_frd = data_frd.UpCont; % nT

%% load 14 day mag data 
% temporal_data_14day = jsondecode(fileread('data/FRD20220116_14days_minute.json')); % January 2022
% 
% for i = 1:(60*24)
%     t_temporal_frd_14day_avg(i) = datetime( ...
%         temporal_data_14day.datetime{i}(1:(end - 5)), ...
%         'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss', ...
%         'TimeZone', 'Etc/UTC');
% end
% t_temporal_frd_14day_avg.TimeZone = "America/New_York";
% 
% 
% H_mag = vecnorm([temporal_data_14day.X, temporal_data_14day.Y, temporal_data_14day.Z]')';
% 
% H_temporal_frd_14day = H_mag - H_core_frd - H_crust_frd; % nT
% t_temporal_frd_14day = (1:length(H_temporal_frd_14day)) * 60; % s
% 
% temp = [];
% 
% for i = 1:14 % 14 days...
%     inds = (i - 1) * 60*24 + (1:60*24);
%     temp = [temp, H_temporal_frd_14day(inds)];
% end
% H_temporal_frd_14day_avg = mean(temp')';
% t_temporal_frd_14day_avg_s = (1:length(H_temporal_frd_14day_avg)) * 60; % s
% 
% inds = ~isnan(H_temporal_frd_14day_avg);
% 
% f_temp = fit( ...
%     t_temporal_frd_14day_avg_s(inds)', H_temporal_frd_14day_avg(inds), 'sin6');


%% plotting
% figure
% hold on
% grid on
% plot(t_temporal_frd_14day_avg, temp)
% plot(t_temporal_frd_14day_avg, H_temporal_frd_14day_avg, 'r', 'Linewidth', 1.5)
% plot(t_temporal_frd_14day_avg, f_temp(t_temporal_frd_14day_avg_s), 'b', 'Linewidth', 1.5)
% xtickformat('HH:mm')
% xlabel('Local Time')
% ylabel('Magnetic Field Strength (nT)')
% 
% empt = repmat({''}, 1, 14);
% legend(empt{:}, 'Averaged Data', 'Model Fit', 'Location', 'Southwest')

% figure
% plot(t_temporal_frd_14day, H_temporal_frd_14day)
% grid on
% hold on

% save('test.mat', 'f_temp', 't_temporal_frd_14day_avg')



%% generate truth + measurements
lla0 = [36.90165855141354, -79.70578451436336, 4e3]; % deg, deg
% lla0 = [37.1979516376812, -79.5787821351352, 4e3]; % deg, deg
% lla0 = [37.0000, -79.9198, 4e3]; % deg, deg
llaf = [37.22891412982679, -80.43067879145124, 4e3]; % deg, deg

Ts = 1;%1/10; % s
[xs_truth, ts_truth, n] = generate_trajectory(lla0, llaf, 50, Ts); %

% t_temporal_data(i_temporal:(i_temporal + n - 1));

i_temporal = 32401; % 4am i think

seconds_offset = seconds(t_temporal_data(i_temporal) - t_temporal_data(1)); % s

H_temporal_frd = ...
    temporal_data.S(i_temporal:(i_temporal + n - 1)) - H_core_frd - H_crust_frd; % nT


zs_truth = read_map(xs_truth, map); % using truth measurements, perfect sensor




%% filter
use_mag = true;

R = 1e2;%1e-5; % measurement noise, very small for perfect sensor, nT
R2 = 1e1; % alt, m
R3 = 1 * eye(3); % heading

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
    [randn(3, 1) * 10; randn(3, 1) * 10; randn(3, 1) * .1]; % initialize with truth

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
fh = @(x) x(4:6) / norm(x(4:6)); % current heading reading, (direction vector)


% loop through measurements -----------------------------
tic
xhat    = x0;
phat    = p0;
NIS     = 0;
zbar    = 0;
grad    = 0;
xs      = nan(length(x0), n);
phats   = nan(length(p0(:)), n);
NISs    = nan(1, n);
grads   = nan(1, n);

for i = 1:n % %%%% loop measurements %%%%%

    iend = i; if isnan(zbar); iend = i - 1; break; end % break if off the map

    xs(:, i)        = xhat;
    phats(:, i)     = phat(:);
    NISs(i)         = NIS;
    grads(i)        = grad;

    % propagate state and covariance
    M = Phi*phat*Phi' + Q; % covariance before update
    xbar = Phi*xhat; % estimate before update

    % Update filter state and covariance matrix using the measurement
    if use_mag
        z = zs_truth(i) + randn*1 + H_temporal_frd(i)*1; % current magnetometer reading, nT

%         z = z - f_temp(seconds_offset + ts_truth(i)); % remove daily variation

        zbar = read_map(xbar, map); % consult the map! nT

        resid = z - zbar; % nT

        H = numerical_jacobian(f, 1, xbar); % state-measurement linearization, nT/m

        grad = norm(H(1:3)); % nT/m

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

    % process heading measurement
    [xhat, phat] = update_heading(xhat, phat, xs_truth, i, fh, R3);

end %%%% loop measurements %%%%%
toc

xx_cov_pos = sqrt(sum(phats([1,  11, 21], 1:iend))); % 1-sigma pos, m
xx_cov_vel = sqrt(sum(phats([31, 41, 51], 1:iend))); % 1-sigma vel, m/s
xx_cov_acc = sqrt(sum(phats([61, 71, 81], 1:iend))); % 1-sigma acc, m/s^2

lla_estimate = ecef2lla(xs(1:3, 1:iend)')';

r_error = vecnorm(xs_truth(1:3, 1:iend) - xs(1:3, 1:iend), 2, 1);
v_error = vecnorm(xs_truth(4:6, 1:iend) - xs(4:6, 1:iend), 2, 1);
a_error = vecnorm(xs_truth(7:9, 1:iend) - xs(7:9, 1:iend), 2, 1);

%% plotting
% clc;close all

fprintf('Mean Error Stats ---------------\n')
fprintf('\tPosition:     %f (m)\n',         mean(r_error))
fprintf('\tVelocity:     %f (m/s)\n',       mean(v_error))
fprintf('\tAcceleration: %f (m/s^2)\n',     mean(a_error))

% h0 = figure;
% h0.WindowStyle = 'Docked';
% plot(t_temporal_data, temporal_data.S - H_core_frd - data_frd.UpCont)
% % plot(t_temporal_data, temporal_data.S)
% grid on
% hold on
% plot(t_temporal_data(1:60:end), f_temp(t_temporal_frd_14day_avg))
% title('24 Hour Magnetic Field (NGS FRD)')
% ylabel('Magnetic Field Strength (nT)')
% xlabel('Time (Eastern Time Zone)')
% % legend('Truth', 'Estimated', 'Location', 'Southwest')
% 
h2 = figure;
h2.WindowStyle = 'Docked';
plot(ts_truth(1:iend), NISs(1:iend))
grid on
xlabel('Time (s)')
ylabel('NIS')
% 
% h3 = figure;
% h3.WindowStyle = 'Docked';
% plot(ts_truth(1:iend), grads(1:iend))
% grid on
% xlabel('Time (s)')
% ylabel('Gradient (nT/m)')
% 
% h4 = figure;
% h4.WindowStyle = 'Docked';
% plot(ts_truth(1:iend), H_temporal_frd(1:iend))
% grid on
% xlabel('Time (s)')
% ylabel('Temporal Effects (nT)')

h1 = plot_mag_countour(map);

figure(h1)
% plot(llaf(2), llaf(1), 'xk')
plot(linspace(lla0(2), llaf(2), n), linspace(lla0(1), llaf(1), n), '--k')
% plot(lla_estimate(2, :), lla_estimate(1, :), 'm')
% legend('', '', '', '', 'Truth', 'Estimate', 'Location', 'Southwest')
grid on
% 
% plot_filter( ...
%     ts_truth(1:iend), r_error, v_error, a_error, xx_cov_pos, xx_cov_vel, xx_cov_acc, ...
%     NISs(1:iend), qtilda, true)







