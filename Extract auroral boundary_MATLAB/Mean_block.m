function [c1,c2,cb]=Mean_block(hphi1, hphi2, I, mask)
% eqn. (22)

c1 = sum(sum(hphi1 .* I)) / sum(hphi1(:));
intersect = (1 - hphi1) .* hphi2;
c2 = sum(sum(intersect .* I)) / sum(intersect(:));
cb = sum(sum((1 - hphi2) .* I)) / sum(sum((1 - hphi2) .* mask));



