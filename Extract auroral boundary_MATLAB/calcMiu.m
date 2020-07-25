function miu = calcMiu(phi1, phi2, hphi1, hphi2)
intersect = (1 - hphi1) .* hphi2 + (1 - hphi2) .* hphi1;
miu = sum(sum(abs(phi1 - phi2) .* intersect)) / sum(intersect(:));