function plot_filter( ...
    time, ...
    xx_kalmanD_err_pos, ...
    xx_kalmanD_err_vel, ...
    xx_kalmanD_err_acc, ...
    xx_covD_pos, ...
    xx_covD_vel, ...
    xx_covD_acc, ...
    NIS, ...
    qtilda, ...
    plot_log)

linew = 1.5;
nx = 4;
ny = 1;
figure('Name', ...
    ['Filter Performance, qtilda = ', num2str(qtilda/(1000^2), '%e'), ' (km^2/s^5)'], ...
    'units','normalized', 'outerposition', [0 0 1 1]);


% position error discrete -------------------------------------------------
subplot(nx, ny, 1)
if plot_log
    semilogy(time, xx_kalmanD_err_pos, ...
        time, xx_covD_pos, 'LineWidth', linew)
else
    plot(time, xx_kalmanD_err_pos, ...
        time, xx_covD_pos, 'LineWidth', linew)
end
grid on
title('Filter Performance')
legend('True Error', 'Covariance')
ylabel('Position Error (m)')

% velocity error discrete -------------------------------------------------
subplot(nx, ny, 2)
if plot_log
    semilogy(time, xx_kalmanD_err_vel, ...
        time, xx_covD_vel, 'LineWidth', linew)
else
    plot(time, xx_kalmanD_err_vel, ...
        time, xx_covD_vel, 'LineWidth', linew)
end
grid on
legend('True Error', 'Covariance')
ylabel('Velocity Error (m/s)')

% accel error discrete ----------------------------------------------------
subplot(nx, ny, 3)
if plot_log
    semilogy(time, xx_kalmanD_err_acc, ...
        time, xx_covD_acc, 'LineWidth', linew)
else
    plot(time, xx_kalmanD_err_acc, ...
        time, xx_covD_acc, 'LineWidth', linew)
end
grid on
legend('True Error', 'Covariance')
ylabel('Acceleration Error (m/s^2)')

% NIS ---------------------------------------------------------------------
subplot(nx, ny, 4)
if plot_log
    semilogy(time, NIS, 'LineWidth', linew)
else
    plot(time, NIS, 'LineWidth', linew)
end
grid on
xlabel('Time (s)')
ylabel('NIS')

end



















