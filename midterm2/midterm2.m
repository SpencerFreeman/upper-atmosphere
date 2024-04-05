clear;clc;close all



%% 2) A He+ ion is in near Earth space, at the geomagnetic equator. It has 
% a mass of 6.64x10-27 kg, and a charge of +1.60x10-19. It has a kinetic 
% energy of 3.2 eV (where 1eV = 1.6x10- 19 J). The particle’s motion is 
% directed at an angle of 70 degrees away from the local geomagnetic field 
% vector. The local geomagnetic field has a magnitude of 4700 nT, and a 
% radius of curvature of 7300 km. Compute the following:
fprintf('\n')

r2cyc = 1/(2*pi); % cycles/radian

eV2J = 1.6e-19; % J/eV

He_mass = 6.64e-27; % kg
He_c = 1.60e-19; % C
He_u = 3.2 *  eV2J; % J 

dpsi = 70 * pi/180; % The particles motion directed away from the local geomagnetic field vector, rad

B_mag = 4700e-9; % local geomagnetic field, T
r = 7300; % radius of curvature, km

% a) The velocity of the particle parallel to the magnetic field, in m/s 
% (hint, use its kinetic energy to find its total speed first)

v_mag = sqrt(2*He_u/He_mass); % m/s
v_B_para = v_mag * cos(dpsi); % m/s

fprintf('2a: %e, m/s\n', v_B_para)

% b) The velocity of the particle perpendicular to the magnetic field, in m/s.

v_B_perp = v_mag * sin(dpsi); % m/s

fprintf('2b: %e, m/s\n', v_B_perp)

% c) The gyroradius of the particle, in m.

r_l = He_mass * v_B_perp / He_c / B_mag; % m

fprintf('2c: %f, m\n', r_l)

% d) The gyroperiod of the particle, in s. (hint, don’t forget 2-pi)

w_c = v_B_perp / r_l * r2cyc; % Hz
T = 1 / w_c; % s

fprintf('2d: %f, s\n', T)

% e) The magnitude of the gradient plus curvature drift of the particle, in m/s.

% 𝑉_∇𝐵+𝑉_(∇×𝐵)=𝑚/𝑞  (𝑅 ⃗_𝑐×𝐵 ⃗_ )/(𝑅_𝑐^2 𝐵^2 ) (𝑣_∥^2+1/2 𝑣_⊥^2 )

v_grad_curv = He_mass / He_c * (v_B_para^2 + .5*v_B_perp^2); % m/s

fprintf('2e: %f, m/s <--- CHECK THIS\n', v_grad_curv)

% f) The magnetic field strength, in nT, at which this particle will mirror.

B_m = B_mag * (v_mag / v_B_perp)^2; % T

fprintf('2f: %f, nT\n', B_m * 10^9)


%% 3) In the E-region ionosphere, intense auroral electrons are creating a 
% Chapman production layer. The altitude of peak production is at 125 km, 
% with a production rate of 3.2*10^8 ions m-3 s-1. The recombination coefficient 
% (alpha) in this region is 0.0079.
fprintf('\n')
clear

% a) Estimate the ion density in this region at the time of the intense auroral electrons.

% b) A few moments later, the auroral electrons move away from this region. 
% If it is currently nighttime (i.e. sun has set), determine the rate of 
% change of ion density 3 seconds after the auroral electrons disappear. 
% (hint, if you solve this numerically, use a small time-step to ensure 
% dt * dn/dt is << n).

























