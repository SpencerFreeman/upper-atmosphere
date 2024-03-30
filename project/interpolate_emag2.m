function z = interpolate_emag2(lat, lon, data, lat_range, lon_range)

xq = linspace(360 + lon_range(1), 360 + lon_range(2), 75);
yq = linspace(lat_range(1), lat_range(2), 70);
x  = data.LON;
y  = data.LAT;
v  = data.UpCont;
[Xq, Yq, Vq] = griddata(x, y, v, xq, yq', "linear");



end