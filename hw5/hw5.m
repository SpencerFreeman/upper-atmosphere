clear;clc;close all


% 1 ---------------------------------

% Boltzmann's constant 
kb = 1.385e-23; % J/K

% Permittivity of free space
ep = 8.85e-12; % F/m

% Charge of an electron 
e_charge = -1.6e-19; % C

% Mass of an electron
e_mass = 9.11e-31; % kg


B = [2,  30,  0]*10^-3; % T

% An electron in this region has a velocity of
v = [400, 200, 300]; % m/s

vp = v - dot(v, B/norm(B)) * B/norm(B); % m/s

rg = e_mass*norm(vp) / abs(e_charge)/norm(B); % m (𝑚𝑣_⊥)/|𝑞|𝐵

rg_nm = rg / 10^-9; % nm

fprintf('1: %f\n', rg_nm)

% 2 ---------------------------------

e_ndensity = 1.3e10; % electrons/m^3

T = 1650; % K

lam = sqrt(ep * kb * T / e_ndensity / e_charge^2); % m

fprintf('2: %f\n', lam*100)

% 3 -------------------------------------

e_ndensity2 = 1.3e10; % electrons per m3.

n_ndensity = 2e16; % neutrals per m3.

sza = 160; % degrees.

ilr = 3e8; % 1/m^3/s, per m3 per second.

beta = ilr / n_ndensity / e_ndensity;

fprintf('3: %f\n', beta)

% 4 -----------------------------------------

G = 6.67430e-11;

earth_mass = 5.97e24; % kg

earth_r = 6378e3; % m

alt = 400e3; % m

r = alt + earth_r; % m

theta = 80; % deg

B0 = 3.1e-5; % T

B2 = [-2 * B0 * earth_r / r^3 * cosd(theta); ...
          -B0 * earth_r / r^3 * sind(theta); ...
                                         0]; % T

Fg = G*earth_mass*e_mass/r^2 * [-1; 0; 0]; % N

Bdir = B2/norm(B2);

Fp = Fg - dot(Fg, Bdir) * Bdir; 

vcg = 1/e_charge * cross(Fp, B2) / (norm(B2)^2);

fprintf('4: %e\n', norm(vcg))

%%

B2d = B2/norm(B2);
Fpd = Fp/norm(Fp);
Fgd = Fg/norm(Fg);

% figure
% plot3([0 B2d(1)], [0 B2d(2)], [0 B2d(3)])
% hold on
% plot3([0 Fpd(1)], [0 Fpd(2)], [0 Fpd(3)])
% plot3([0 Fgd(1)], [0 Fgd(2)], [0 Fgd(3)])
% 
% grid on
% axis equal
% xlabel('r')
% ylabel('theta')
% zlabel('psi')



















