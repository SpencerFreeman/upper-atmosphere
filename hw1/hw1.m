clear;clc;close all



data = importfile1("Homework_1_input_SLE_012224.xlsx", "Input", [2, Inf]);


comp = importfile2("Homework_1_input_SLE_012224.xlsx", "Input", [2, Inf]);

% 1
% honor code

% 2
% Mean distance between particles d=1/∛𝑛
alt = data.HeightKm; % km
d = (data.OCm3 + data.N2Cm3 + data.O2Cm3 + data.ArCm3) .^ (-1/3); % cm

m2cm = 10^2; % cm/m
m2nm = 10^9; % nm/m
m2km = 10^-3; % km/m
cm32m2 = (1/m2cm)^3; % m^3/cm^3

d_130k = interp1(alt, d, 130) / m2cm * m2nm % nm

% 3
d_400k = interp1(alt, d, 400) / m2cm * m2nm % nm

% 4
kb = 1.38e-23; % m^2 * kg / s^2 / K
fv = @(v, m, T) sqrt(2/pi) * (m/kb/T)^(3/2) * v.^2 .*exp(-m*v.^2 /2/kb/T);

v = 0:30:2500; % m/s
T = 298.5; % K
m_He = 6.64647907*10^-27; % kg
m_Ne = 3.35092*10^-26; % kg
m_Ar = comp{4, 4}; % kg
m_Xe = 2.18017*10^-25; % kg

fv_He = fv(v, m_He, T);
fv_Ne = fv(v, m_Ne, T);
fv_Ar = fv(v, m_Ar, T);
fv_Xe = fv(v, m_Xe, T);

% plot(v, fv_He, v, fv_Ne, v, fv_Ar, v, fv_Xe);
% grid on
% xlabel('Speed (m/s)')
% ylabel('PDF (s/m)')

T_400k = interp1(alt, data.TemperatureK, 400); % K
m_O = comp{1, 4}; % kg
v_int = 1000:1:2500; % m/s

fv_O = fv(v_int, m_O, T_400k);

prob_O = trapz(v_int, fv_O) * 100 % probability, %

% 5
fv_Ar = fv(v_int, m_Ar, T_400k);

prob_Ar = trapz(v_int, fv_Ar) * 100 % probability, %

% 6
% 𝑙_1,2=𝑐_1/𝜈_1,2 =1/(𝜎_1,2 𝑛_2 √(1+(𝑚_1 𝑇_2 )∕(𝑚_2 𝑇_1 ) ))
% l_1,1 = 1/(𝜎_1,1 n √2)

n = interp1(alt, data.N2Cm3, 300); 
r = comp{2, 2}; % m
sig_11 = 4*pi*r^2; % m

l_11 = 1/(sig_11*n*sqrt(2)) * m2km % m

% 7
rho_ball = 1.8e-6; % kg/m^2
% [~, ~, ~, rho_145k] = atmosisa(145e3);

m_N2 = comp{2, 4}; % kg
m_O2 = comp{3, 4}; % kg
rho_145k = ...
    (interp1(alt, data.OCm3,  145) * m_O + ...
    interp1(alt, data.N2Cm3, 145) * m_N2 + ...
    interp1(alt, data.O2Cm3, 145) * m_O2 + ...
    interp1(alt, data.ArCm3, 145) * m_Ar) / cm32m2; % kg/m^3

r_b = 3*rho_ball/rho_145k % m

% 8
ep = .589;
r_O = comp{1, 2}; % m
sig_11_O = 4*pi*r_O^2; % m
T_500k = interp1(alt, data.TemperatureK, 500); % K
eta = ep*sqrt(m_O*kb*T_500k)/sig_11_O
















