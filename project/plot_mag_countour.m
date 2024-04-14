function [map, h] = plot_mag_countour(data, lat_range, lon_range)

xq = linspace(360 + lon_range(1), 360 + lon_range(2), 75);
yq = linspace(lat_range(1), lat_range(2), 70);
x  = data.LON;
y  = data.LAT;
v  = data.UpCont;
[Xq, Yq, Vq] = griddata(x, y, v, xq, yq', "linear");
map(:, :, 1) = Xq;
map(:, :, 2) = Yq;
map(:, :, 3) = Vq;

h = figure('WindowStyle', 'Docked');
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

end