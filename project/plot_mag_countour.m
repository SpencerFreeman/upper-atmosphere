function h = plot_mag_countour(map)

Xq = map(:, :, 1); % lon
Yq = map(:, :, 2); % lat
Vq = map(:, :, 3); % mag



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

clim([min(Vq(:)) max(Vq(:))])
c = colorbar;
c.Label.String = '(nT)';

lon_plot = linspace(Xq(1, 1) - 360, Xq(1, end) - 360, 5)';
for i = 1:length(lon_plot)
    temp = degrees2dms(lon_plot(i));
    lblsx{i} = [ num2str(temp(1), '%.0f'), '\circ', ...
                num2str(temp(2), '%.0f'), '''', ...
                num2str(temp(3), '%.0f'), ''''''];
end
xticks(lon_plot)
xticklabels(lblsx)

lat_plot = linspace(Yq(1, 1), Yq(end, 1), 3)';
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