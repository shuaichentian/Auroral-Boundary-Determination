function dH = Dirac(phi, epsilon)
dH = 0.5 / epsilon * ( 1 + cos(phi * pi / epsilon)); % eqn. (24)
dH(phi >= epsilon) = 0;
dH(phi <= -epsilon) = 0;
end