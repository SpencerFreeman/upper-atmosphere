function map = build_map(data, lat_range, nlat, lon_range, nlon)

xq = linspace(360 + lon_range(1), 360 + lon_range(2), nlon);
yq = linspace(lat_range(1), lat_range(2), nlat);
x  = data.LON;
y  = data.LAT;
v  = data.UpCont;
[Xq, Yq, Vq] = griddata(x, y, v, xq, yq', "linear");
map(:, :, 1) = Xq;
map(:, :, 2) = Yq;
map(:, :, 3) = Vq;

end