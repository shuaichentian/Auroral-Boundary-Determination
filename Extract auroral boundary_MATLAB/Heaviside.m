function H=Heaviside(phi, epsilon)

H = 0.5 * (1 + phi / epsilon + 1/pi * sin(phi * pi / epsilon));  % eqn. (23)
H(phi > epsilon) = 1;
H(phi < -epsilon) = 0;
