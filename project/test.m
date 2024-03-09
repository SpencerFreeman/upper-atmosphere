clear;clc;close all


data = import_emag2_file("C:\Users\spenc\Desktop\school shit 2\Spring_2024\upper_atmosphere\project\data\EMAG2_V3_20170530\EMAG2_V3_20170530.csv", [1, Inf]);


%%

lat_range = [36 38];
lon_range = [-81 -78];

inds = data.LON > (360 + lon_range(1)) & data.LON < (360 + lon_range(2)) & ...
    data.LAT > lat_range(1) & data.LAT < lat_range(2);

data_sub = data(inds, :);


%%
clc; close all

xq = linspace(360 + lon_range(1), 360 + lon_range(2), 75);
yq = linspace(lat_range(1), lat_range(2), 70);
x  = data_sub.LON;
y  = data_sub.LAT;
v  = data_sub.UpCont;

[Xq, Yq, Vq] = griddata(x, y, v, xq, yq', "linear");

figure
contourf(Xq, Yq, Vq)
clim([min(v) max(v)])
c = contourcbar;
c.Label.String = '(nT)';
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
zlabel('Magnetic Anomaly (nT)')














