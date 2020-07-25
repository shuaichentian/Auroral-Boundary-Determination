%% 计算定量评价指标
function [Pd, Pf, Pmp] = MeasureAccuracy(Pa, Pm, Iauto, Imanual)

Na = size(Pa, 2);
Nm = size(Pm, 2);

da = zeros(Na, 1);
dm = zeros(Nm, 1);
for i = 1 : Na
    x = Pa(1, i) - Pm(1, :);
    y = Pa(2, i) - Pm(2, :);
    d = sqrt(x.^2 + y.^2);
    da(i) = min(d);
end

for i = 1 : Nm
    x = Pm(1, i) - Pa(1, :);
    y = Pm(2, i) - Pa(2, :);
    d = sqrt(x.^2 + y.^2);
    dm(i) = min(d);
end

%% Pixel deviation(Pd)
Pd = (sum(da) + sum(dm)) / (Nm + Na);

%% The percentage of distant pixels(Pf)
Tfar = 7;
Fa = length(find(da > Tfar));
Fm = length(find(dm > Tfar));
Pf = (Fa + Fm) / (Nm + Na) * 100;


%% percentage of mislabelled piexels(Pmp)
Pmp = sum(sum(abs(Imanual - Iauto))) / (sum(Imanual(:)) + sum(Iauto(:))) * 100;
