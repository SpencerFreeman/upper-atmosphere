% This function solves for the qE_star and h_star of a Chapman Production
% Function to fit a curve to real data
function [qE_star, h_star, H] = GetChapmanFit(Zenith, StarVars_hat, h, qE)

options = optimoptions(...
    'fsolve', ...
    'FunctionTolerance', 1E-5, ...
    'OptimalityTolerance', 1E-5, ...
    'MaxIterations', 1000, ...
    'display', 'none', ...
    'Diagnostics', 'off');

ObjectiveFunction = @(StarVars) ChapmanObjectiveFunction(Zenith, StarVars, h, qE);
[solution, ~] = fsolve(ObjectiveFunction, StarVars_hat, options);

%note, this is going to give you some error warnings that you can ignore for our purposes here
qE_star = solution(1);
h_star = solution(2);
H = solution(3);


end