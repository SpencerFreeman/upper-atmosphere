% This function is the objectve function which minimizes the difference
% between the actual qE data and the qE fitted from the chapman production
% function.

function [ErrorNorm] = ChapmanObjectiveFunction(Zenith, StarVars, h, qE)

qE_star_max = StarVars(1);
h_star_max = StarVars(2);
H = StarVars(3);

% Evaluate Chapman Production function at Initial Guess
for z = 1:length(h)
    tao(z) = secd(Zenith)*exp(-(h(z) - h_star_max)/H); % Optical Depth
    ChapmanFit_qE(z,1) = qE_star_max*exp(1 - (h(z) - h_star_max)/H - tao(z)); % Energy Deposition Rate per Unit Volume (Eq 3.35 in text)
end

Error = abs(qE - ChapmanFit_qE);
ErrorNorm = norm(Error,2);
