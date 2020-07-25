    imageDate='31-Jul-1997_09_53_45'
    %~strcmp([imageDate ],filename1)
    %if ~strcmp([imageDate ],filename1)
     %   continue
    %end

    disp('loading image');
    % load_image = 'images/test.bmp';
    load_image = ['img/' imageDate '.bmp'];
    mask_image = 'mask.bmp'; % fov mask
    image = imread(load_image);
    image_o=image;
    [m, n] = size(image); 
    image = histeq(image);
    mask = imread(mask_image);
    mask = double(mask);
    loadimg=['logits_img/' imageDate '.bmp'];
    img_log=imread(loadimg);
    img_log1=img_log
    img_log = double(img_log);
    img_log=img_log/255.0-0.5;
    flag=load('flag.mat');
    flag=struct2cell(flag);
    flag=cell2mat(flag);
    flag1=flag;
    flag2=flag;
    flag1(flag>30)=0;
    flag1(flag<30 | flag==30)=1;
    flag2(flag>70)=0;   %120 before changed
    flag2(flag<70 | flag==70)=1;
    phi1 = SDF(flag1);
    phi2 = SDF(flag2);
    %flag1 = 1 - flag1 * 2;
    %flag2 = imread(['images/010_012021_outer.bmp']);
    %flag2 = im2bw(flag2);
%flag2 = flag2 * 2 - 1;

    fMaxDist = m * n;

% figure;imshow(image);
% h = fspecial('average', 5);
    h = fspecial('gaussian', 3, 1);
    image = imfilter(image, h);
    figure;imshow(image);
    Img = mat2gray(image);
% figure;imshow(uint8(Img));


%%%%%%%% Parameters setting %%%%%%%
    rad = 13; % the radius of local information
    alpha = 5; % controls the smoothness of segmenting shapes
    beta = 1.2; % weights the significance of the shape prior knowledge
    gama = 1; % weights the significance of the local information   '''1 before changed'''
    lamda = 2.0; %  is a coefficient of the global region energy which helps speed up the curve evolution
    sigma = 30; % the variance of shape distance  %20 before cahnged
    deltaT = 0.5; % time step
    epsilon = 10; % narrowband width
    lamda2=1;  % the weights of the logits term.  2 before changed



%%%%%%%% curve initialization %%%%%%%
% flag  = roipoly(mat2gray(Img));
% flag1 = flag * 2 - 1;
%phi1 = -SDF(flag1);
% flag  = roipoly(mat2gray(Img));
% flag2 = flag * 2 - 1;
%phi2 = SDF(flag2);

    tic;
    for i = 1 : 1000
    
    %重新初始化SDF
    if mod(i-1, 10) == 0
%         flag = ones(m, n);
%         flag(phi1 < 0) = -1;
%         phi1 = sdf_mex(single(flag),single(fMaxDist),int32(1),int32(1));
%         flag = ones(m, n);
%         flag(phi2 < 0) = -1;
%         phi2 = sdf_mex(single(flag),single(fMaxDist),int32(1),int32(1));
        
            figure(3);
            imagesc(img_log1,[0, 255]); axis off; axis equal; colormap(gray);
            title(num2str(i));
            hold on;
            contour(phi1, [0,0], 'r', 'LineWidth', 3);
            contour(phi2, [0,0], 'g', 'LineWidth', 3);
            %contour(phi1, [0,0], 'b', 'LineWidth', 3);
            hold off;
            pause(0.001);
    end
    
%     [phi1, phi2] = evolution(phi1, phi2, Img, mask, alpha, beta, gama, lamda, sigma, deltaT, epsilon);
        [phi1, phi2] = evolution2(i,phi1, phi2, Img,img_log, mask, rad, alpha, beta, gama, lamda, lamda2, sigma, deltaT, epsilon);

    %显示轮廓曲线
%     figure(4);mesh(phi1);
%     figure(5);mesh(phi2);


    end
    toc;

    Iauto = zeros(size(image));
    Iauto(phi1 <= 0 & phi2 >=0) = 1;
%Pa = extractBoundary(Iauto);
%load(['结果/' imageDate '_manual.mat']);
%Pm = Points;
%Imanual = imread(['结果/' imageDate '_manual.bmp']);
%Imanual = im2bw(Imanual);
%[Pd, Pf, Pmp] = MeasureAccuracy(Pa, Pm, Iauto, Imanual);
%disp('Pd  Pf  Pmp');
%disp([Pd, Pf, Pmp]);
    imwrite(Iauto, ['LIDLSM/' imageDate '_LIDLSM.bmp']);