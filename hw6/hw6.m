clear;clc;close all


% Boltzmann's constant 
kb = 1.385e-23; % J/K

% Permittivity of free space
ep = 8.85e-12; % F/m

% Charge of an electron 
e_charge = -1.6e-19; % C

% Mass of an electron
e_mass = 9.11e-31; % kg

%% 2
% You are developing a radio system to transmit over the horizon. You have 
% a working transmitter, that is fixed at 3.00x108 Hz, and you are trying 
% to see if this can work for your application.
% 
% Using IRI, you have determined the typical maximum electron densities in 
% the D, E and F-layers of the ionosphere around your location to be D=5.6x109, 
% E=4.3x1011 and F=4.9x1014 m-3.
% 
% Determine the following:
% 
% The minimum angle with respect to vertical that you will need to transmit 
% to ensure it is refracted back to the ground, in degrees.

f = 3.00e8; % Hz

n_D = 5.6e9; % # / m^3 
n_E = 4.3e11; % # / m^3 
n_F = 4.9e14; % # / m^3
ns = [n_D, n_E, n_F];

% 𝜔_𝑝=((𝑛_𝑒 𝑞_𝑒^2)/(𝑚_𝑒 𝜀_0 ))^(1⁄2), radians/s

wp_D = sqrt(n_D *  e_charge^2 / (e_mass * ep)); % rad/s
wp_E = sqrt(n_E *  e_charge^2 / (e_mass * ep)); % rad/s
wp_F = sqrt(n_F *  e_charge^2 / (e_mass * ep)); % rad/s
wps = [wp_D, wp_E, wp_F];

[wp_max, ind] = max(wps); % rad/s
MUF = f * 2*pi; % rad/s

theta_min = acosd(wp_max / MUF); % deg

fprintf('2: %f (deg)\n', theta_min)

%% 3
% You set up the system described above and are able to communicate with a 
% receiver at a remote location across the horizon. However, on the occasion 
% the attenuation coefficient (alpha) is 70%, which is too much for your 
% design parameters. Your system requires the attenuation coefficient alpha 
% to be no more than 0.5. If nothing else changed, what transmit frequency 
% would you need to achieve this attenuation, in MHz. 

alpha_max = .5; % max attenuation coefficient

k = 1.16e-3;

alpha_cur = .7; % current

v_collision = alpha_cur * MUF^2 / k / ns(ind);

% alpha=(1.16×10^(−3) 𝑛_𝑒 𝜈_𝑐𝑜𝑙𝑙𝑖𝑠𝑖𝑜𝑛)/(𝜔_𝑟𝑎𝑑𝑖𝑜^2 ), dB/km

w_radio = sqrt(k * ns(ind) * v_collision / alpha_max); % rad/s

fprintf('3: %f (MHz)\n', w_radio / (2*pi) / 1e6)













