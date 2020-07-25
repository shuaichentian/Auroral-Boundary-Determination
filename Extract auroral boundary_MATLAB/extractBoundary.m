function Points = extractBoundary(Im)

% �������߷���ȡ��������߽�

[h, w] = size(Im);
% imshow(Im);
[yy, xx] = find(Im > 0);
centerX = fix(mean(xx));
centerY = fix(mean(yy));
% [centerX, centerY] = ginput(1);
% centerX = fix(centerX);
% centerY = fix(centerY);
flag = zeros(size(Im));

% �ϱ߽�
% disp('�ϱ߽�');
y = 1;
for x = 1 : w
    [ x_coord, y_coord ] = bresenham_line(x, y, centerX, centerY);
    ind = sub2ind([h, w], y_coord, x_coord);
    subIm = Im(ind);
    k = find(subIm > 0, 1);
    flag(x_coord(k), y_coord(k)) = 1;
end

% �±߽�
% disp('�±߽�');
y = h;
for x = 1 : w
    [ x_coord, y_coord ] = bresenham_line(x, y, centerX, centerY);
    ind = sub2ind([h, w], y_coord, x_coord);
    subIm = Im(ind);
    k = find(subIm > 0, 1);
    flag(x_coord(k), y_coord(k)) = 1;
end

% ��߽�
% disp('��߽�');
x = 1;
for y = 1 : h
    [ x_coord, y_coord ] = bresenham_line(x, y, centerX, centerY);
    ind = sub2ind([h, w], y_coord, x_coord);
    subIm = Im(ind);
    k = find(subIm > 0, 1);
    flag(x_coord(k), y_coord(k)) = 1;
end

% �ұ߽�
% disp('�ұ߽�');
x = w;
for y = 1 : h
    [ x_coord, y_coord ] = bresenham_line(x, y, centerX, centerY);
    ind = sub2ind([h, w], y_coord, x_coord);
    subIm = Im(ind);
    k = find(subIm > 0, 1);
    flag(x_coord(k), y_coord(k)) = 1;
end

[rows, cols] = find(flag>0);
% hold on;
% plot(rows, cols, '.r');
% figure; imshow(flag);

Points = [rows, cols]';
% save([filename '.mat'], 'Points');

% hold on; 
% [C, h] = contour(Im, [0.5 0.5], '-r');
% hold off;
% Points = C(:,2:end);
% save([filename '.mat'], 'Points');
% bw = im2bw(Im);
% figure;imshow(bw);
% bw2 = bwperim(bw, 8);
% figure;imshow(bw2);
% 
% 
% bw = roipoly(Im);
% 
% imwrite(bw, [filename '_bound.bmp']);