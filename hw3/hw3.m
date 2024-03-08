clear;clc;close all


data = importfile("C:\Users\spenc\Desktop\school shit 2\Spring_2024\upper_atmosphere\hw3\Homework_3_input_SLE_020624.xlsx", "Sheet1", [2, Inf]);


%% 1


n_O = data.Ocm3; % n/cm^3
alt_msis = data.MSISAltkm; % km

ind = find(alt_msis == 165);

dndz = (n_O(ind + 1) - n_O(ind - 1)) / (alt_msis(ind + 1) - alt_msis(ind - 1)); 

Hn = (abs(dndz)/n_O(ind))^-1; % km

fprintf('Q1: %f (km)\n', Hn)

%% 2

sza = data.ICONDiskSolarZenithAngleDegrees;
inds = ~isnan(sza);
sza = sza(inds);

tvc_brightness = data.BrightnessR(inds);

temp = sortrows([sza, tvc_brightness]);

sza = temp(:, 1);
tvc_brightness = temp(:, 2);

figure
plot(sza, tvc_brightness, 'o')
grid on
xlabel('SZA (deg)')
ylabel('Brightness (Q)')



f = polyfit(sza, tvc_brightness, 1);

cs = cosd(sza);
fk = polyfit(cs, tvc_brightness, 1);


hold on


plot([0 sza(end)], [f(2), f(2) + f(1)*sza(end)], '-x')

sz_test = linspace(0, sza(end), 100);
plot(sz_test, fk(2) + fk(1)*cosd(sz_test), '-')
plot(0, fk(2) + fk(1)*cosd(0), 'xy')


k_tvc_brightness = tvc_brightness./cosd(sza);

figure
plot(cs, tvc_brightness)
hold on
plot(cs,  fk(2) + fk(1)*cs)
grid on
xlabel('SZA (deg)')
ylabel('Brightness (Q)')


fprintf('Q2: %f (R)\n', fk(2) + fk(1)*cosd(0))


%% 3

alt_icon = data.ICONLimbAltituideKm; % km

inds = ~isnan(alt_icon);

alt_icon = alt_icon(inds);

icon_brightness = data.BrightnessR1(inds); % q

figure
plot(icon_brightness, alt_icon, 'o')
grid on
xlabel('Brightness (Q)')
ylabel('Altitude (km)')

[q_max, ind_q_max] = max(icon_brightness);

fprintf('Q3: %f (R)\n', q_max)


%% 4

chi = 0; % rad

z = 200; % km

zmax = alt_icon(ind_q_max); % km

q_200km = q_max*exp(1 - (z - zmax)/Hn - sec(chi)*exp(-(z - zmax)/Hn));

hold on
plot(q_200km, z, 'xr')


zs = linspace(min(alt_icon), max(alt_icon), 100); %km
qs = q_max*exp(1 - (zs - zmax)/Hn - sec(chi)*exp(-(zs - zmax)/Hn));

plot(qs, zs, 'r')

fprintf('Q4: %f (R)\n', q_200km)


%% 5

% Load Data from Excel sheet
% Data = readmatrix('YOURFILEHERE.xlsx');
h = alt_icon; %Data(:,1); % [km] height
qE = icon_brightness; %Data(:,2); % [1/vol/s] energy deposition
Zenith = chi; % [deg] solar zenith angle

% Find Initial Guess from data of star values for solver
qE_star_max_hat = q_max; % max energy deposition rate
h_star_max_hat = zmax; % Max absorption height 
H = Hn; % [km] scale height 

% Evaluate Chapman Production function at Initial Guess
for z = 1:length(h)
    tao(z) = secd(Zenith)*exp(-(h(z) - h_star_max_hat)/H); % Optical Depth
    ChapmanInitGuess_qE(z,1) = qE_star_max_hat*exp(1 - (h(z) - h_star_max_hat)/H - tao(z)); % Energy Dchapmaneposition Rate per Unit Volume (Eq 3.35 in text)
end

% Initial Guess Vector
StarVars_hat = [qE_star_max_hat; h_star_max_hat; H];

disp('-------------------Fitted Values--------------------')
[qE_star_max, h_star_max, H] = GetChapmanFit(Zenith, StarVars_hat, h, qE);

% [ErrorNorm] = ChapmanObjectiveFunction(Zenith, StarVars_hat, h, qE)

qs_star = qE_star_max*exp(1 - (zs - h_star_max)/H - sec(chi)*exp(-(zs - h_star_max)/H));

plot(qs, zs, 'b--')

fprintf('Q5: %f (km)\n', H)


%% 6

sza_ref = 120; % deg

Re = 6378; % km

alt_dark = Re/cosd(sza_ref - 90) - Re; % km

fprintf('Q6: %f (km)\n', alt_dark)


























