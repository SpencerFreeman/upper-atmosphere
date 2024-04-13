clear;clc;close all

% Solution Output:
% Spencer Freeman, 4/12/2024
% 
% 2a: 4247.345825, m/s
% 2b: 11669.486746, m/s
% 2c: 103.039085, m
% 2d: 0.055479, s
% 2e: 0.104177, m/s
% 2f: 5322.629358, nT
% 
% 3a: 201261.842171, ions/something...
% 3b: 42.185130, 42.204214, 42.185249, ions/something... / s
% 
% 4: -3.781406e-12, Pa/km

%% 2) A He+ ion is in near Earth space, at the geomagnetic equator. It has 
% a mass of 6.64x10-27 kg, and a charge of +1.60x10-19. It has a kinetic 
% energy of 3.2 eV (where 1eV = 1.6x10- 19 J). The particle’s motion is 
% directed at an angle of 70 degrees away from the local geomagnetic field 
% vector. The local geomagnetic field has a magnitude of 4700 nT, and a 
% radius of curvature of 7300 km. Compute the following:
fprintf('Solution Output:\n')
fprintf('Spencer Freeman, 4/12/2024\n\n')

r2cyc = 1/(2*pi); % cycles/radian

eV2J = 1.6e-19; % J/eV

He_mass = 6.64e-27; % kg
He_c = 1.60e-19; % C
He_u = 3.2 *  eV2J; % J 

dpsi = 70 * pi/180; % The particles motion directed away from the local geomagnetic field vector, rad

B_mag = 4700e-9; % local geomagnetic field, T
Rc_mag = 7300e3; % radius of curvature, m

% a) The velocity of the particle parallel to the magnetic field, in m/s 
% (hint, use its kinetic energy to find its total speed first)

v_mag = sqrt(2*He_u/He_mass); % m/s
v_B_para = v_mag * cos(dpsi); % m/s

fprintf('2a: %f, m/s\n', v_B_para)

% b) The velocity of the particle perpendicular to the magnetic field, in m/s.

v_B_perp = v_mag * sin(dpsi); % m/s

fprintf('2b: %f, m/s\n', v_B_perp)

% c) The gyroradius of the particle, in m.

r_l = He_mass * v_B_perp / He_c / B_mag; % m

fprintf('2c: %f, m\n', r_l)

% d) The gyroperiod of the particle, in s. (hint, don’t forget 2-pi)

w_c = v_B_perp / r_l * r2cyc; % Hz
T = 1 / w_c; % s

fprintf('2d: %f, s\n', T)

% e) The magnitude of the gradient plus curvature drift of the particle, in m/s.

% 𝑉_∇𝐵+𝑉_(∇×𝐵)=𝑚/𝑞  (𝑅 ⃗_𝑐×𝐵 ⃗_ )/(𝑅_𝑐^2 𝐵^2 ) (𝑣_∥^2+1/2 𝑣_⊥^2 )

B = B_mag * [0; 1; 0]; % B field vector (at geomagnetic equator), T
Rc = Rc_mag * [1; 0; 0]; % m

v_grad_curv = norm(He_mass / He_c * cross(Rc, B) / (Rc_mag^2 * B_mag^2) * (v_B_para^2 + .5*v_B_perp^2)); % m/s

fprintf('2e: %f, m/s\n', v_grad_curv)

% f) The magnetic field strength, in nT, at which this particle will mirror.

B_m = B_mag * (v_mag / v_B_perp)^2; % T

fprintf('2f: %f, nT\n', B_m * 10^9)


%% 3) In the E-region ionosphere, intense auroral electrons are creating a 
% Chapman production layer. The altitude of peak production is at 125 km, 
% with a production rate of 3.2*10^8 ions m-3 s-1. The recombination coefficient 
% (alpha) in this region is 0.0079.
fprintf('\n')
clear

q = 3.2*10^8; % ions / m^3 / s
alpha = 0.0079;

% a) Estimate the ion density in this region at the time of the intense auroral electrons.

% 𝑛_𝑒 (𝑧,𝑡)≈√(𝑞(𝑧,𝑡)/𝛼)

n_e0 = sqrt(q / alpha); % ions / 

fprintf('3a: %f, ions/something...\n', n_e0)

% b) A few moments later, the auroral electrons move away from this region. 
% If it is currently nighttime (i.e. sun has set), determine the rate of 
% change of ion density 3 seconds after the auroral electrons disappear. 
% (hint, if you solve this numerically, use a small time-step to ensure 
% dt * dn/dt is << n).

% (𝜕𝑛_𝑒)/𝜕𝑡=𝑞(𝑧,𝑡)−𝛼𝑛_𝑒^2

n_e = n_e0;

tend = 3; % s
dt = 1e-6; % s
ts = linspace(0, tend, round(tend/dt)); % s


for i = 1:length(ts)

    n_e_dot = -alpha * n_e^2;

    n_e = n_e + n_e_dot * dt;

    n_es(i) = n_e;

end

[t,nes] = ode45(@(t, ne) -alpha * ne^2, [0 tend], n_e0);
n_e_ode = nes(end);

n_e_tend = 1 / (alpha * tend + 1 / n_e0);

% figure
% plot(ts, n_es)
% grid on

fprintf('3b: %f, %f, %f, ions/something... / s\n', n_e, n_e_ode, n_e_tend)




%% 4) You have obtained observations of the O+ density in the upper F-region 
% near the magnetic pole. At this altitude, the ion temperature is 
% approximately constant with height, and the electron temperature is 1.5x 
% the ion temperature. 

% Find the vertical gradient in the pressure associated with the O+ plasma 
% in this region. You may ignore any other sources of pressure, and you may 
% assume that O+ acts as an ideal gas when computing pressure (see Part 1, 
% Lesson 1 slides if needed). 
fprintf('\n')
clear

n_O = 17312751370; % O+ density, m^-3

T_O = 1606; % Ion temperature, K
T_e = T_O * 1.5; % electron temperature, K
h = 584; % Altitude of observation, km

Re = 6378; % Radius of Earth, km
mass_O = 2.66e-26; % Mass of O+, kg
G = 6.67e-11; % Gravitational constant G 
mass_earth = 5.97e24; % Mass of Earth, kg
kb = 1.38e-23; % Boltzmann's constant

% use iterative method to calculate density
km2m = 1000; % m/km
m2km = 1/km2m; % km/m
delz = .1*km2m; % m

z2 = h:(delz*m2km):(h + 100); % km

rho0 = n_O * mass_O; % kg/m^3
rho = rho0;

for i = 2:length(z2)
    zi = z2(i);

    r = Re + zi; % km
    g = G*mass_earth/(r*km2m)^2; % m/s^2
    H = kb*T_O/mass_O/g;
    rho(i) = rho(i - 1)*exp(-delz/H); % kg/m^3
end
N = rho/mass_O; % number of ions / m^3

% ideal gas law
p = N*kb*T_O; % Pa

% forward difference
dpdz = diff(p(1:2)) / diff(z2(1:2)); % Pa/km

figure
plot(p, z2)
grid on
xlabel('Pressure (Pa)')
ylabel('alt (km)')

fprintf('4: %e, Pa/km\n', dpdz)





