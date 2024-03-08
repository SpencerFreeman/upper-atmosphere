clear;clc;close all

m2cm = 10^2; % cm/m
m2nm = 10^9; % nm/m
m2km = 10^-3; % km/m
km2m = 10^3; % m/km
cm32m3 = (1/m2cm)^3; % m^3/cm^3
m32cm3 = m2cm^3; % cm^3/m^3

data = importfile1("C:\Users\spenc\Desktop\school shit 2\Spring_2024\upper_atmosphere\hw2\Homework_2_input_SLE_020524.xlsx", "Sheet1", [2, Inf]);

m_Ar = data.VarName8(1); % kg
m_N2 = 4.65E-26; % kg
kb = data.VarName8(4); % kg
G = data.VarName8(2); % 
m_earth = data.VarName8(3); % kg
r_earth = 6378; % km

figure
plot(data.TemperatureK, data.AltitudeKm, 'o-')
grid on
xlabel('Temperature (K)')
ylabel('Alt (km)')


%% 1

n = data.N2Cm3(1:3); % /cm^3
z = data.AltitudeKm2(1:3); % km

ind = 2;

dndz = (n(ind + 1) - n(ind - 1)) / (z(ind + 1) - z(ind - 1)); % n/km

Hn = (abs(dndz)/n(ind))^-1; % km
H = Hn;

r = r_earth + z(ind); % km
g = G*m_earth/(r*km2m)^2; % m/s^2

disp('Q1 answer:')
Tinf = H*km2m*m_N2*g/kb % ANSWER, K

%% 2

zt  = 300; % m
z0 = 120; % m
Tz0 = interp1( data.AltitudeKm, data.TemperatureK, z0); % K
k = (data.TemperatureK(end) - data.TemperatureK(end - 1))/ ...
    (data.AltitudeKm(end) - data.AltitudeKm(end - 1)); % K/km

z = 120:5:300; % km

disp('Q2 answer:')
s = k/(Tinf - Tz0) % ANSWER

tz = Tinf - (Tinf - Tz0)*exp(-s*(z - z0));

dtdz = s*(Tinf - Tz0)*exp(-s*(z - z0));

figure
plot(tz, z)
grid on
xlabel('Temp (K)')
ylabel('Alt (km)')

figure
plot(dtdz, z)
grid on
xlabel('dTdz (K)')
ylabel('Alt (km)')


%% 3

zi = 210; %km

disp('Q3 answer:')
tz = Tinf - (Tinf - Tz0)*exp(-s*(zi - z0)) % ANSWER, K



%% 4

delz = 1*km2m; % m

z2 = 120:(delz*m2km):210; % km

rho0 = data.ArCm3(1) * m_Ar / cm32m3; % kg/m^3

rho = rho0;

for i = 2:length(z2)

    zi = z2(i);

    T = Tinf - (Tinf - Tz0)*exp(-s*(zi - z0));

    r = r_earth + zi; % km

    g = G*m_earth/(r*km2m)^2; % m/s^2

    H = kb*T/m_Ar/g;

    rho(i) = rho(i - 1)*exp(-delz/H); % kg/m^3

end

n = rho/m32cm3/m_Ar; % n/cm^3

figure
plot(n, z2)
grid on
xlabel('rho (n/cm^3)')
ylabel('alt (km)')

disp('Q4 answer:')
n(end) % ANSWER, n/cm^3


%% 5

rho0 = interp1(z2, rho, 160);
rho = rho0;

z3 = 160:(delz*m2km):260; % km
ntot = 0;
for i = 2:length(z3)

    zi = z3(i);

    T = Tinf - (Tinf - Tz0)*exp(-s*(zi - z0));

    r = r_earth + zi; % km

    g = G*m_earth/(r*km2m)^2; % m/s^2

    H = kb*T/m_Ar/g;

    rho(i) = rho(i - 1)*exp(-delz/H); % kg/m^3

    n = rho(i)/m_Ar; % n/m^3

    ntot = ntot + n*delz*1*1; 

end

disp('Q5 answer:')
ntot % ANSWER





%% 6+7

clear

% Quick commands to open the file in Matlab

latitude = ncread("ICON_L2-4_FUV_Day-Limb_2021-04-10_v05r000.NC","ICON_L24_Latitude");
altitude = ncread("ICON_L2-4_FUV_Day-Limb_2021-04-10_v05r000.NC","ICON_L24_Altitude");
longitude = ncread("ICON_L2-4_FUV_Day-Limb_2021-04-10_v05r000.NC","ICON_L24_Longitude");
O = ncread("ICON_L2-4_FUV_Day-Limb_2021-04-10_v05r000.NC","ICON_L24_O");
O2 = ncread("ICON_L2-4_FUV_Day-Limb_2021-04-10_v05r000.NC","ICON_L24_O2");
N2 = ncread("ICON_L2-4_FUV_Day-Limb_2021-04-10_v05r000.NC","ICON_L24_N2");
temp = ncread("ICON_L2-4_FUV_Day-Limb_2021-04-10_v05r000.NC","ICON_L24_Temperature");

n = N2(:, 3); % /cm^3
z = altitude(:, 3); % km

zi = 156.4; % km

ind = find(abs(z - zi) < .001);

dndz = (n(ind + 1) - n(ind - 1)) / (z(ind + 1) - z(ind - 1)); % n/km

disp('Q6 answer:')
Hn = (abs(dndz)/n(ind))^-1 % km


%% 7

r_earth = 6378; % km
G = 6.67E-11;
m_earth = 5.97E+24; % kg
kb = 1.38E-23;

T = interp1(z, temp(:, 3), zi); % K
m_N2 = 4.65E-26; % kg
r = r_earth + zi; % km
g = G*m_earth/(r*1000)^2; % m/s^2

H = kb*T/m_N2/g/1000; % km

disp('Q7 answer:')
d = (Hn - H)/H * 100 % ANSWER, percent

















