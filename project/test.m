clear;clc;close all


raw_data = import_emag2_file("data\EMAG2_V3_20170530\EMAG2_V3_20170530.csv", [1, Inf]);


%%

lat_range = [ 36.758719  37.556273];
lon_range = [-80.874399 -79.464463];

inds = raw_data.LON > (360 + lon_range(1)) & raw_data.LON < (360 + lon_range(2)) & ...
    raw_data.LAT > lat_range(1) & raw_data.LAT < lat_range(2);

data = raw_data(inds, :);

% save('data/EMAG2_V3_Blacksburg-Roanoke', 'data')

%%

% 38.21666, -77.3833 --> pretty close to NGS FRD

range = .01; % deg

lat_range = 38.21116410524037 + range*[-1, 1]; % deg
lon_range = -77.37559992408342 + range*[-1, 1]; % deg

inds = raw_data.LON > (360 + lon_range(1)) & raw_data.LON < (360 + lon_range(2)) & ...
    raw_data.LAT > lat_range(1) & raw_data.LAT < lat_range(2);

data_frd = raw_data(inds, :);

save('data/EMAG2_V3_NGS-Fredericksburg', 'data_frd')


%% plotting
clc; close all

load('EMAG2_V3_Blacksburg-Roanoke', 'data')

xq = linspace(360 + lon_range(1), 360 + lon_range(2), 75);
yq = linspace(lat_range(1), lat_range(2), 70);
x  = data.LON;
y  = data.LAT;
v  = data.UpCont;
[Xq, Yq, Vq] = griddata(x, y, v, xq, yq', "linear");

figure('WindowStyle', 'Docked')
contourf(wrapTo180(Xq), Yq, Vq)%, "ShowText",true,"LabelFormat","%0.0f (nT)")
hold on

% plot(wrapTo180(x), y, '.k')
% plot(wrapTo180(Xq(1, :)), Yq(:, 1), '.k')

clr = [1 0 0];

blacksburg_lla = [37.230139960377606, -80.4293122107635];
plot(blacksburg_lla(2), blacksburg_lla(1), '*', 'Color',clr )
text(blacksburg_lla(2), blacksburg_lla(1) + .1, 'Blacksburg', 'Color',clr)

roanoke_lla = [37.2810435498791, -79.95796154021077];
plot(roanoke_lla(2), roanoke_lla(1), '*', 'Color',clr)
text(roanoke_lla(2), roanoke_lla(1) + .1, 'Roanoke', 'Color',clr)

clim([min(v) max(v)])
c = colorbar;
c.Label.String = '(nT)';

lon_plot = linspace(lon_range(1), lon_range(2), 5)';
for i = 1:length(lon_plot)
    temp = degrees2dms(lon_plot(i));
    lblsx{i} = [ num2str(temp(1), '%.0f'), '\circ', ...
                num2str(temp(2), '%.0f'), '''', ...
                num2str(temp(3), '%.0f'), ''''''];
end
xticks(lon_plot)
xticklabels(lblsx)

lat_plot = linspace(lat_range(1), lat_range(2), 3)';
for i = 1:length(lat_plot)
    temp = degrees2dms(lat_plot(i));
    lblsy{i} = [ num2str(temp(1), '%.0f'), '\circ', ...
                num2str(temp(2), '%.0f'), '''', ...
                num2str(temp(3), '%.0f'), ''''''];
end
yticks(lat_plot)
yticklabels(lblsy)

axis equal
title('NOAA EMAG2 V3 Magnetic Anomaly')
xlabel('Longitude')
ylabel('Latitude')
zlabel('Magnetic Anomaly (nT) ')




%%
clear;clc;close all

syms x y z vx vy vz

r = [x;y;z];
v = [vx;vy;vz];


% jacobian(norm(r), r)

Hv = jacobian(v / norm(v), v)



% i = 200;
% data(i, :)
% V = data(i, :).UpCont
% Vi = griddata(x, y, v, data(i, :).LON, data(i, :).LAT, "linear")

% plot(data(i, :).LON-360,data(i, :).LAT,'*k')














