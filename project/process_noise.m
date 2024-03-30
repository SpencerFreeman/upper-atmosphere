function Q = process_noise(qtilda, T)

I = eye(3);
Q = zeros(9,9);
Q(1:3,1:3) = .05*T^5*I;
Q(1:3,4:6) = .125*T^4*I;
Q(1:3,7:9) = T^3/6*I;
Q(4:6,1:3) = .125*T^4*I;
Q(4:6,4:6) = T^3/3*I;
Q(4:6,7:9) = .5*T^2*I;
Q(7:9,1:3) = T^3/6*I;
Q(7:9,4:6) = .5*T^2*I;
Q(7:9,7:9) = T*I;
Q = Q*qtilda;

end