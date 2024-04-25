function [xhat, phat] = update_heading(xbar, M, xs_truth, i, f, R)

% Update filter state and covariance matrix using the measurement (heading)
z = xs_truth(4:6, i) / norm(xs_truth(4:6, i)); % current heading reading, (direction vector)

zbar = xbar(4:6) / norm(xbar(4:6)); % predicted heading

resid = z - zbar; % m

H = numerical_jacobian(f, 3, xbar); % state-measurement linearization, m/m

Sj = H*M*H' + R;
K = M*H'*inv(Sj); % calculate optimal gain
IKH = eye(9) - K*H;
phat = IKH*M*IKH' + K*R*K'; % update covariance (Joseph form)

xhat = xbar + K*resid; % update estimate

S = R + H*phat*H'; % residual covariance
NIS = resid'*(S\resid); % normalized innovation (residual) squared

end