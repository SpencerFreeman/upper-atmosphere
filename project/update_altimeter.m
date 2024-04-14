function [xhat, phat] = update_altimeter(xbar, M, xs_truth, i, fa, R2)

% Update filter state and covariance matrix using the measurement (altitude)
z = norm(xs_truth(:, i)); % current altimeter reading, m

zbar = norm(xbar(1:3)); % predicted altitude, m

resid = z - zbar; % m

H = numerical_jacobian(fa, 1, xbar); % state-measurement linearization, m/m

Sj = H*M*H' + R2;
K = M*H'*inv(Sj); % calculate optimal gain
IKH = eye(9) - K*H;
phat = IKH*M*IKH' + K*R2*K'; % update covariance (Joseph form)

xhat = xbar + K*resid; % update estimate

S = R2 + H*phat*H'; % residual covariance
NIS = resid'*(S\resid); % normalized innovation (residual) squared

end