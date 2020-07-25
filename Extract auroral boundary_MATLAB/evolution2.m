function [phi1, phi2] = evolution2(i,phi1, phi2, I,img_log, mask, rad, alpha, beta, gama, lamda, lamda2, sigma, deltaT, epsilon)
% EVOLUTION2 水平集演化
% 输入：
%      phi1: 内水平集函数
%      phi2: 外水平集函数
%      I   : 输入图像
%      mask: FOV 掩模

%-- 计算phi1和phi2的heviside函数
hphi1 = Heaviside(phi1, epsilon);
hphi2 = Heaviside(phi2, epsilon);

dist = phi2 - phi1;
% miu = calcMiu(phi1, phi2, hphi1, hphi2);    % eqn.(9)
miu = 13; 
sigma2 = 2 * sigma^2;  % 2 * sigma^2

% dhphi1 = Dirac(phi1, epsilon);
% dhphi2 = Dirac(phi2, epsilon);

%-- 获取窄带坐标
[idx1, x1, y1] = get_narrowband(phi1);
[idx2, x2, y2] = get_narrowband(phi2);

 %-- 计算三个区域的平均灰度
[c1, c2, cb] = Mean_block(hphi1, hphi2, I, mask);


%-- phi1
% globalTerm1 = dhphi1(idx1) .*((I(idx1) - c1).^2 - hphi2(idx1)  .* (I(idx1) - c2).^2);
globalTerm1 = (I(idx1) - c1).^2 - hphi2(idx1)  .* (I(idx1) - c2).^2;
globalTerm1 = globalTerm1 / (max(abs(globalTerm1)) + eps);

%logitsTerm
logits_term1 = img_log(idx1);

localTerm1 = get_localTerm(I, phi1, mask, idx1, x1, y1, rad);
% localTerm1 = localTerm1 .* dhphi1(idx1);
localTerm1 = localTerm1 / (max(abs(localTerm1)) + eps);

% shapeTerm1 = (1 - 2 .* hphi2(idx1)) .* (dist(idx1) - miu).^2 ./ sigma2 - 2 * (dist(idx1) - miu) .* (hphi1(idx1) + hphi2(idx1) - 2 * hphi1(idx1) .* hphi2(idx1)) / sigma2;
shapeTerm1 = - 2 * (dist(idx1) - miu) / sigma2;
shapeTerm1 = shapeTerm1 / (max(abs(shapeTerm1)) + eps);


% regTerm1 = dhphi1(idx1) .* get_curvature(-phi1,idx1,x1,y1);  
regTerm1 = get_curvature(-phi1,idx1,x1,y1);  
regTerm1 = regTerm1 / (max(abs(regTerm1)) + eps);

dphidt1 =   gama * localTerm1+ beta * shapeTerm1 + alpha * regTerm1 + lamda2*logits_term1+lamda * globalTerm1;
% dt1 = .45/(min(dphidt1)+eps);
dt1 = T(i,idx1,img_log,deltaT);
phi1(idx1) = phi1(idx1) - dt1 .* dphidt1;  
 %-- Keep SDF smooth
 phi1 = sussman(phi1, .5);

%--phi2
% globalTerm2 = dhphi2(idx2) .* ((1 - hphi1(idx2)) .* (I(idx2) - c2).^2 - (I(idx2) - cb).^2);
globalTerm2 = (1 - hphi1(idx2)) .* (I(idx2) - c2).^2 - (I(idx2) - cb).^2; 
globalTerm2 = globalTerm2 / (max(abs(globalTerm2)) + eps);

%logitsTerm
logits_term2=img_log(idx2);

localTerm2 = get_localTerm(I, phi2, mask, idx2, x2, y2, rad);
% localTerm2 = localTerm2 .*  dhphi2(idx2);
localTerm2 = localTerm2 / (max(abs(localTerm2)) + eps);

% shapeTerm2 = (1 - 2 .* hphi1(idx2)) .* (dist(idx2) - miu).^2 ./ sigma2 + 2 * (dist(idx2) - miu) .* (hphi1(idx2) + hphi2(idx2) - 2 * hphi1(idx2) .* hphi2(idx2)) / sigma2;
shapeTerm2 = 2 * (dist(idx2) - miu) / sigma2;
shapeTerm2 = shapeTerm2 / (max(abs(shapeTerm2)) + eps);

% regTerm2 = dhphi2(idx2) .* get_curvature(-phi2,idx2,x2,y2);  
regTerm2 = get_curvature(-phi2,idx2,x2,y2);  
regTerm2 = regTerm2 / (max(abs(regTerm2)) + eps);

dphidt2 =  gama * localTerm2 + beta * shapeTerm2 + alpha * regTerm2-lamda2*logits_term2+lamda * globalTerm2;
% dt2 = .45/(max(dphidt2)+eps);
dt2 = T(i,idx2,img_log,deltaT);
phi2(idx2) = phi2(idx2) - dt2 .* dphidt2;  
 %-- Keep SDF smooth
 phi2= sussman(phi2, .5);
 
    
function t=T(i,idx,img_log,deltaT)
if i<160
    t=deltaT;
    %t=( max(( -0.499.*(img_log(idx)+0.1).*(img_log(idx)+0.9).*6.25), 0 )+0.001 );
else
    t=deltaT;
end



%-- get the curves's narrow band
function [idx, x, y] = get_narrowband(phi)
idx = find(phi <= 1.3 & phi >= -1.5)';
[y, x] = ind2sub(size(phi),idx);


%-- compute local stats
function localTerm = get_localTerm(I, phi,mask, idx, x, y, rad)
%-- get windows for localized statistics
[dimy, dimx] = size(I);
xneg = x-rad; xpos = x+rad;      %get subscripts for local regions
yneg = y-rad; ypos = y+rad;
xneg(xneg<1)=1; yneg(yneg<1)=1;  %check bounds
xpos(xpos>dimx)=dimx; ypos(ypos>dimy)=dimy;

%-- initialize u,v
u=zeros(size(idx)); v=zeros(size(idx));
Ain=zeros(size(idx)); Aout=zeros(size(idx)); 
%-- compute local stats
for i = 1:numel(idx)  % for every point in the narrow band
    img = I(yneg(i):ypos(i),xneg(i):xpos(i)); %sub image
    P = phi(yneg(i):ypos(i),xneg(i):xpos(i)); %sub phi
    M = mask(yneg(i):ypos(i),xneg(i):xpos(i)); %sub mask
    PM = P .* M;
%     upts = find(P >= 0 & M > 0);            %local interior
    upts = find(PM > 0);            %local interior

    Ain(i) = length(upts)+eps;
    u(i) = sum(img(upts))/Ain(i);
    
%     vpts = find(P < 0 & M > 0);             %local exterior
    vpts = find(PM < 0);             %local exterior

    Aout(i) = length(vpts)+eps;
    v(i) = sum(img(vpts))/Aout(i);
end

% localTerm = -(u - v) .* (2 .* I(idx) - u - v);
localTerm = -((u-v).*((I(idx)-u)./Ain+(I(idx)-v)./Aout));

%-- compute curvature along SDF
function curvature = get_curvature(phi,idx,x,y)
    [dimy, dimx] = size(phi);        

    %-- get subscripts of neighbors
    ym1 = y-1; xm1 = x-1; yp1 = y+1; xp1 = x+1;

    %-- bounds checking  
    ym1(ym1<1) = 1; xm1(xm1<1) = 1;              
    yp1(yp1>dimy)=dimy; xp1(xp1>dimx) = dimx;    

    %-- get indexes for 8 neighbors
    idup = sub2ind(size(phi),yp1,x);    
    iddn = sub2ind(size(phi),ym1,x);
    idlt = sub2ind(size(phi),y,xm1);
    idrt = sub2ind(size(phi),y,xp1);
    idul = sub2ind(size(phi),yp1,xm1);
    idur = sub2ind(size(phi),yp1,xp1);
    iddl = sub2ind(size(phi),ym1,xm1);
    iddr = sub2ind(size(phi),ym1,xp1);
    
    %-- get central derivatives of SDF at x,y
    phi_x  = -phi(idlt)+phi(idrt);
    phi_y  = -phi(iddn)+phi(idup);
    phi_xx = phi(idlt)-2*phi(idx)+phi(idrt);
    phi_yy = phi(iddn)-2*phi(idx)+phi(idup);
    phi_xy = +0.25*phi(iddl)+0.25*phi(idur)...
             -0.25*phi(iddr)-0.25*phi(idul);  %!!!!!!
    phi_x2 = phi_x.^2;
    phi_y2 = phi_y.^2;
    
    %-- compute curvature (Kappa)
    curvature = ((phi_x2.*phi_yy + phi_y2.*phi_xx - 2*phi_x.*phi_y.*phi_xy)./...
              (phi_x2 + phi_y2 +eps).^(3/2)).*(phi_x2 + phi_y2).^(1/2);        
  

%-- level set re-initialization by the sussman method
function D = sussman(D, dt)
  % forward/backward differences
  a = D - shiftR(D); % backward
  b = shiftL(D) - D; % forward
  c = D - shiftD(D); % backward
  d = shiftU(D) - D; % forward
  
  a_p = a;  a_n = a; % a+ and a-
  b_p = b;  b_n = b;
  c_p = c;  c_n = c;
  d_p = d;  d_n = d;
  
  a_p(a < 0) = 0;
  a_n(a > 0) = 0;
  b_p(b < 0) = 0;
  b_n(b > 0) = 0;
  c_p(c < 0) = 0;
  c_n(c > 0) = 0;
  d_p(d < 0) = 0;
  d_n(d > 0) = 0;
  
  dD = zeros(size(D));
  D_neg_ind = find(D < 0);
  D_pos_ind = find(D > 0);
  dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                       + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
  dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                       + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;
  
  D = D - dt .* sussman_sign(D) .* dD;
  
%-- whole matrix derivatives
function shift = shiftD(M)
  shift = shiftR(M')';

function shift = shiftL(M)
  shift = [ M(:,2:size(M,2)) M(:,size(M,2)) ];

function shift = shiftR(M)
  shift = [ M(:,1) M(:,1:size(M,2)-1) ];

function shift = shiftU(M)
  shift = shiftL(M')';
  
function S = sussman_sign(D)
  S = D ./ sqrt(D.^2 + 1);    

  

