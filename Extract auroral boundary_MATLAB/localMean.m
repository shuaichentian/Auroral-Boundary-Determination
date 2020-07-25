function [u,v] = localMean(I, hphi, radius, mask)

[m, n] = size(I);
u = zeros(m, n);
v = zeros(m, n);

for i = radius + 1 : m - radius
    for j = radius + 1 : n - radius
        
        if(hphi(i, j) > 0 && hphi(i, j) < 1)       
            subHphi1 = hphi(i-radius:i+radius,j-radius:j+radius);
            subI = I(i-radius:i+radius,j-radius:j+radius);
            subMusk = mask(i-radius:i+radius,j-radius:j+radius);
            if sum(subHphi1(:)) > 0
                u(i, j) = sum(sum(subHphi1 .* subI .* subMusk)) / sum(sum(subHphi1  .* subMusk));
                v(i, j) = sum(sum((1-subHphi1) .* subI .* subMusk)) / sum(sum(1-subHphi1  .* subMusk));
            end
        end
    end
end

end