clear;clc;close all

%% soultion output
% Spencer Freeman, 3/1/24

% Q11: 5.026548e-19 (m^2)
% Q12: 51.430651 (1/s)
% Q13: 890.316469 (m/s)
% Q14: 7.024940 (m)
% Q15: 333986.018339 
% Q2: 191.340879 (km) 
% 
% Q31: Richardson Number [-0.341711] is >0 due to >0 N2, 
% thus likely unstable since Ri < 1/4. 
% 
% Q32: Richardson Number [-0.079082] is >0 due to >0 N2, 
% thus likely unstable since Ri < 1/4. 


%% 1
% 1

r_O = 1e-10; % m

r_N2 = 3e-10; % m

sig_ON2 = pi*(r_O + r_N2)^2; % m^2

fprintf('Q11: %d (m^2)\n', sig_ON2)

% 2
kb = 1.38e-23; % m2 kg s-2 K-1

m1 = 2.66e-26; % kg
m2 = 4.65e-26; % kg

n2 = 8.9e16; % m-3

T1 = 600; % K
T2 = 700; % K

T12 = (m2*T1 + m1*T2)/(m1 + m2);
m12 = m1*m2/(m1 + m2);

v = sig_ON2*n2*sqrt(8*kb*T12/pi/m12); % 1/s

fprintf('Q12: %f (1/s)\n', v)

% 3

cbar_O = sqrt(8*kb*T1/pi/m1); % m/s

fprintf('Q13: %f (m/s)\n', cbar_O)

% 4

sig_N2 = pi*(r_N2 + r_N2)^2; % m^2

l_N2 = 1/(sig_N2*n2*sqrt(2)); % m

fprintf('Q14: %f (m)\n', l_N2)

% 5
sig_O = pi*(r_O + r_O)^2; % m^2

n1 = 1.5e16; % m-3

l_O = 1/(sig_O*n1*sqrt(2)); % m

D12 = cbar_O*l_O;

fprintf('Q15: %f \n', D12)

%% 2
% 1
clear;close all

kb = 1.38e-23; % m2 kg s-2 K-1

r_N2 = 3e-10; % m
m_N2 = 4.65e-26; % kg

sig_N2 = pi*(r_N2 + r_N2)^2; % m^2



Tinf = 700; % K
Tz0 = 120; % K
s = 0.01/1000; % m-1
g = 9; % m/s^2

delz = 1000; % m

zs = 120e3:(delz):300e3; % m
z0 = 120e3; % m

rho0 = 2e17 * m_N2; % kg/m^3

rho = rho0;

for i = 2:length(zs)

    zi = zs(i);
    T = Tinf - (Tinf - Tz0)*exp(-s*(zi - z0));

    H = kb*T/m_N2/g;

    rho(i) = rho(i - 1)*exp(-delz/H); % kg/m^3

    n = rho(i)/m_N2; % n/m^3

    Kn(i) = 1/sig_N2/n/sqrt(2) * m_N2*g/kb/T; 


end

plot(Kn, zs/1000)
grid on
xlabel('Kn')
ylabel('Alt (km)')

z_kn_1 = interp1(Kn, zs, 1)/1000; % km

fprintf('Q2: %f (km) \n', z_kn_1)


%% 3
clear

g = 9.5; % m/s^2

cp = 1.4e3; % m2s-1K-1

T = 190; % K


% region 1

dTdz = -9/1000; % K/m
dudz = -18/1000; % m/s / m

dvdz = 0; % entirely east in both directions.

N2 = g/T*(dTdz + g/cp);

Ri = N2/(dudz^2 + dvdz^2);

fprintf('\nQ31: Richardson Number [%f] is >0 due to >0 N2, \nthus likely unstable since Ri < 1/4. \n\n', Ri)

% region 2

aK = 1e-3; % W/Km

dudz = sqrt(0.0014);%  m/s/m

kl = 1.45e-2;% W/km
kh = 1.38e-2; % W/km

Tl = (kl/aK)^2; % K
Th = (kh/aK)^2; % K

dTdz = (Th - Tl)/(5000); % K/m

T = (Tl + Th)/2; % K

N2_2 = g/T*(dTdz + g/cp);
Ri_2 = N2/(dudz^2 + dvdz^2);

fprintf('Q32: Richardson Number [%f] is >0 due to >0 N2, \nthus likely unstable since Ri < 1/4. \n', Ri_2)









