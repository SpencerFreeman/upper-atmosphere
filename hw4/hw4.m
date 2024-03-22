clear;clc;close all


ICON_L31_Meridional_Wind = ncread("ICON_L3-1_MIGHTI_2022-06-09-v01r000.NC","ICON_L31_Meridional_Wind");
ICON_L31_Zonal_Wind = ncread("ICON_L3-1_MIGHTI_2022-06-09-v01r000.NC","ICON_L31_Zonal_Wind");
ICON_L31_Temperature = ncread("ICON_L3-1_MIGHTI_2022-06-09-v01r000.NC","ICON_L31_Temperature");
ICON_L31_Altitude = ncread("ICON_L3-1_MIGHTI_2022-06-09-v01r000.NC","ICON_L31_Altitude");



km2m = 1000; % m/km

Re = 6378; % km
G = 6.67430e-11;
m_earth = 5.9722e24; % kg

cp = 1.4e3;


indx = [2 3];
indy = 40:49;

mwind = ICON_L31_Meridional_Wind(indx, indy); % m/s
zwind = ICON_L31_Zonal_Wind(indx, indy); % m/s
T = ICON_L31_Temperature(indx, indy); % K
alt = ICON_L31_Altitude(indx); % km

dz = (alt(2) - alt(1)) * km2m; % m

g = G*m_earth/((mean(alt) + Re)*km2m)^2; % m/s^2

dTdz = (T(2, :) - T(1, :)) / dz; % K/m

N2 = g./T(2, :).*(dTdz + g/cp);

dudz = (mwind(2, :) - mwind(1, :)) / dz; % m/s / m
dvdz = (zwind(2, :) - zwind(1, :)) / dz; % m/s / m

Ri = N2./(dudz.^2 + dvdz.^2);

[~, ind] = min(Ri);

indy(ind)

% 2

[~, ind] = max(Ri);

indy(ind)


%% 3
clear

data = importfile("C:\Users\spencer.freeman\Desktop\projects\school\upper-atmosphere\hw4\HWK4_Mars_Pressure_data_SLE_031124.xlsx", "Sheet1", [1, 8877]);

t = data.VarName1;
p = data.VarName2;

wavemodel = 'A0 + A1*cos(1*t*(2*pi) - phi1) ';%fit equation
f = fittype(wavemodel, 'independent', 't', 'dependent', 'p','coefficients', {'A0','A1', 'phi1'});%parameters to fit
start_points = [600, 1,1]; % Initial guess
[fitresult1, gof] = fit(t, p, f, 'Start', start_points);%do the fit


fitresult1.A1

wavemodel = 'A0 + A1*cos(1*t*(2*pi)/.25 - phi1) ';%fit equation
f = fittype(wavemodel, 'independent', 't', 'dependent', 'p','coefficients', {'A0','A1', 'phi1'});%parameters to fit

[fitresult2, gof] = fit(t, p, f, 'Start', start_points);%do the fit

fitresult2.A1

fitresult2.A1/fitresult1.A1



















