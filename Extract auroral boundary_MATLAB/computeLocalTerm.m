function F = computeLocalTerm( I,dhphi,u,v,radius )
[m, n] = size(I);
% F = zeros(size(I));
% for i = radius + 1 : m - radius
%     for j = radius + 1 : n - radius
%         
%         subDhphi = dhphi(i-radius:i+radius,j-radius:j+radius);
%         subI = I(i-radius:i+radius,j-radius:j+radius);
% %         subu = u(i-radius:i+radius,j-radius:j+radius);
% %         subv = v(i-radius:i+radius,j-radius:j+radius);
%         
%         F(i, j) = sum(sum(((subI - u(i, j)).^2 - (subI - v(i, j)).^2) .* subDhphi));
%         
%     end
% end

F = -(u-v).*(2.* I-u-v);

F = F .* dhphi;

maxF = max(max(abs(F)));
F = F / maxF;

end

