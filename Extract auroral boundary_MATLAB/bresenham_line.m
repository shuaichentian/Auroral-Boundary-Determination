function [ xx, yy ] = bresenham_line( Px0, Py0, Px1, Py1 )
% point = input('Give coord[ x0 y0 x1 y1]: ');
if(abs(Py1 - Py0) > abs(Px1 - Px0))         % If the line is steep
    x0 = Py0; y0 = Px0; x1 = Py1;y1 = Px1;  % then it would be converted to
    token = 1;                              % non steep by changing coordinate
else
    x0 = Px0; y0 = Py0; x1 = Px1; y1 = Py1;
    token = 0;
end

reverse = 0;
if(x0 >x1)
    temp1 = x0; x0 = x1; x1 = temp1;
    temp2 = y0; y0 = y1; y1 = temp2;
    reverse = 1;
end
dx = abs(x1 - x0) ;                              % Distance to travel in x-direction
dy = abs(y1 - y0);                               % Distance to travel in y-direction
sx = sign(x1 - x0);                              % sx indicates direction of travel in X-dir
sy = sign(y1 - y0);                              % Ensures positive slope line

x = x0; y = y0;                                  % Initialization of line
param = 2*dy - dx ;                              % Initialization of error parameter
for i = 0:dx-1                                   % FOR loop to travel along X
    x_coord(i+1) = x;                            % Saving in matrix form for plot
    y_coord(i+1) = y;
%     if (token ==0)                               % Plotting in point form
%         plot(x,y,'r*');                          % For steep line coordinate is again
%     else                                         % converted for actual line drawing.
%         plot(y,x,'r*');
%     end
    param = param + 2*dy;                        % parameter value is modified
    if (param >0)                                % if parameter value is exceeded
        y = y +1*sy;                             % then y coordinate is increased
        param = param - 2*(dx );                 % and parameter value is decreased
        
    end
    x = x + 1*sx;                                % X-coordinate is increased for next point
end

if (token == 0)
    xx = x_coord;
    yy = y_coord;
else
    xx = y_coord;
    yy = x_coord;
end

if (reverse == 1)
    xx = flip(xx);
    yy = flip(yy);
end

