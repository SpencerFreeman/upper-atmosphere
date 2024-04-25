function [xhat, phat] = update_altimeter(xbar, M, xs_truth, i, f, R)

% Update filter state and covariance matrix using the measurement (altitude)
z = norm(xs_truth(1:3, i)); % current altimeter reading, m

zbar = norm(xbar(1:3)); % predicted altitude, m

resid = z - zbar; % m

% H = numerical_jacobian(fa, 1, xbar); % state-measurement linearization, m/m

H = xbar' / norm(xbar);

Sj = H*M*H' + R;
K = M*H'*inv(Sj); % calculate optimal gain
IKH = eye(9) - K*H;
phat = IKH*M*IKH' + K*R*K'; % update covariance (Joseph form)

xhat = xbar + K*resid; % update estimate

S = R + H*phat*H'; % residual covariance
NIS = resid'*(S\resid); % normalized innovation (residual) squared

end