function varargout = numerical_jacobian(f, nf, varargin)
% NUMERICAL_JACOBIAN Numerical Jacobian routine.
% -------------------------------------------------------------------------
% Author: Spencer Freeman, E21
% Date: 12/12/22
%
%   Input:
%       f - function handle, must accept column vector input
%       nf - function output size
%       varargin - [x01, x02, ... x0n] where x0n is the n'th input column
%       vector in the order that they are accepted by the function handle.
%
%   Output:
%       varargout - [M1, M2, ... Mn] where Mn is the n'th Jacobian matrix
%       with respec to the input vectors.
%
% Linearizes the nonlinear input equation y = f(x). The solution takes the
% form dy = M*dx. The output is the matrix M (or Mn if decomposed into
% variables, which is(are) the numerical Jacobian(s) of the input function
% f(x).
% 
% Example Usage:
% 
%     f = @(x) ...
%         ...
%         [x(1)^2 + 3*x(2) + x(3); ...
%         sin(x(3))];
%     
%     x0 =  [2, 0, -pi];
%     
%     C = numerical_jacobian(f, 2, x0)
%     
%     C_truth = ...
%         [   2*x0(1),  3,  0; ...
%             0,          0,  cos(x0(3))]
%     
%     %% alternate split output
%     x0_1 = x0(1:2);
%     x0_2 = x0(3);
%     
%     [A, B] = numerical_jacobian(f, 2, x0_1, x0_2)
%     
%     A_truth = ...
%         [   2*x0_1(1),  3; ...
%             0,          0]
%     
%     B_truth = ...
%         [   1; ...
%             cos(x0_2(1))]
% 
% For more information, see: <a href="matlab: 
% web('https://mathworld.wolfram.com/Jacobian.html')">Jacobian -- from Wolfram MathWorld</a>.
%
% -------------------------------------------------------------------------

nargs = length(varargin);
x0 = [];
nxi = zeros(1, nargs);
for k = 1:nargs
    x0 = [x0; varargin{k}(:)]; % unpack inputs
    nxi(k) = length(x0);
end
nxi = [0, nxi];
nx = length(x0);

del = 1;%10e-12; % numerical perturbation size
xPlus = repmat(x0, 1, nx);
xMinus = xPlus;
dels = del*eye(nx);
xPlus = xPlus + dels; % perturb up
xMinus = xMinus - dels; % perturb down

fPlus = zeros(nf, nx);
fMinus = fPlus;
for j = 1:nx
    fPlusTemp = f(xPlus(:, j)); % evaluate function up
    fPlus(:, j) = fPlusTemp(:); % enforce dimensions
    
    fMinusTemp = f(xMinus(:, j)); % evaluate function down
    fMinus(:, j) = fMinusTemp(:); % enforce dimensions
end

M = (fPlus - fMinus)/del/2; % central difference

varargout = cell(1, nargs);
for k = 1:nargs
    varargout{k} = M(:, nxi(k) + 1:nxi(k + 1)); % package outputs
end

end


