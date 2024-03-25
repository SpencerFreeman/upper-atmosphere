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


rg = e_mass*dot(v, B/norm(B)) / abs(e_charge)/norm(B); % m (𝑚𝑣_⊥)/|𝑞|𝐵

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

earth_r = 6378; % km

alt = 400; % km

r = ((alt + earth_r)*1000)^3; % m

theta = 80; % deg

B0 = 3.1e-5; % T

B2 = [-2 * B0 * earth_r*1000 / r * cosd(theta); ...
    -B0 * earth_r*1000 / r * sind(theta)]; % T

Fg = -G*earth_mass*e_mass/r^2 * [1; 0; 0]; % N

% Fp = Fp * 

% vcg = cross(Fp, B2) / norm(B)^2;



























