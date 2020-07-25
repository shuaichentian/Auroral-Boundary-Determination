function [phi1, phi2] = evolution(phi1, phi2, I, mask, alpha, beta, gama, lamda, sigma, deltaT, epsilon)

sigma2 = 2 * sigma^2;  % 2 * sigma^2

hphi1 = Heaviside(phi1, epsilon);
hphi2 = Heaviside(phi2, epsilon);
dhphi1 = Dirac(phi1, epsilon);
dhphi2 = Dirac(phi2, epsilon);
kappa1 = kappa(-phi1);
kappa2 = kappa(-phi2);
dist = phi2 - phi1;
% miu = calcMiu(phi1, phi2, hphi1, hphi2);    % eqn.(9)
miu = 15;    % eqn.s(9)
[c1, c2, cb] = Mean_block(hphi1, hphi2, I, mask); % 曲线内部和外部区域均值
radius = 9;
[u1, v1] = localMean(I, hphi1, radius, mask);
[u2, v2] = localMean(I, hphi2, radius, mask);



Eglobal1 = dhphi1 .* ((I - c1).^2 - hphi2 .* (I - c2).^2); % eqn. (13)
maxEglobal = max(max(abs(Eglobal1)));
Eglobal1 = Eglobal1 / maxEglobal;
Elocal1 = computeLocalTerm(I, dhphi1, u1, v1,radius); % local energy
Ereg1 = dhphi1 .* kappa1; % eqn. (14)
Eshape1 = dhphi1 .* (1 - 2 .* hphi2) .* (dist - miu).^2 ./ sigma2 - 2 * (dist - miu) .* (hphi1 + hphi2 - 2 * hphi1 .* hphi2) / sigma2; % eqn. (15)
% Eshape1 = dhphi1 .* (1 - 2 .* hphi2) .* (dist - miu).^2 ./ sigma2; % eqn. (15)
% Eext1 = dhphi1 .* (hphi2 - 1) .* g; % eqn. (16)

Eglobal2 = dhphi2 .* ((1 - hphi1) .* (I - c2).^2 - (I - cb).^2); % eqn. (17)
maxEglobal = max(max(abs(Eglobal2)));
Eglobal2 = Eglobal2 / maxEglobal;
Elocal2 = computeLocalTerm(I, dhphi2, u2, v2, radius);
Ereg2 = dhphi2 .* kappa2; % eqn. (18)
Eshape2 = dhphi2 .* (1 - 2 .* hphi1) .* (dist - miu).^2 ./ sigma2 + 2 * (dist - miu) .* (hphi1 + hphi2 - 2 * hphi1 .* hphi2) / sigma2; % eqn. (19)
% Eshape2 = dhphi2 .* (1 - 2 .* hphi1) .* (dist - miu).^2 ./ sigma2; % eqn. (19)
% Eext2 = dhphi2 .* (hphi1 - 1) .* g; % eqn. (20)

phi1 = phi1 - deltaT * (Eglobal1 + alpha * Ereg1 + beta * Eshape1 + gama * Elocal1);  % eqn. (12)
phi2 = phi2 - deltaT * (Eglobal2 + alpha * Ereg2 + beta * Eshape2 + gama * Elocal2);