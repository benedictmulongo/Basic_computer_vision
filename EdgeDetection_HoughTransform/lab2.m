%addpath lab_library/;
% Add the directory and its subdirectories 
%addpath(genpath('lab_library/'))

function [] = lab2(p,q)
%     Question1()
%     Question2()
%     Question3()
%     Question5()
%     Question7()
%     Question7_best('tools')
%     Question7_best('house')
    Question8()

end 

function [curves] = lines(rho, theta,minimum,maximum, N)
X = linspace(minimum,maximum,N);
Y = (rho - X*cos(theta))./sin(theta);
curves = [X;Y];
curves =  round([X;Y], 0);
end

function [] = Question8()
% ntheta = 180;
% nrho = 100;
ntheta = 256;
nrho = 256;
gradmagnthreshold = 10;
verbose = 'all';
nlines = 8;
scale = 1;
pic = triangle128;
fprintf('Size of Image');
disp(size(pic))
% pic = houghtest256;
[linepar, acc] = houghedgeline(pic, scale, gradmagnthreshold, nrho, ntheta, nlines, verbose);
end

function [linepar, acc] = houghedgeline(pic, scale, gradmagnthreshold, nrho, ntheta, nlines, verbose)

[curves, magnitude] =  extractedgeHough(pic, scale, gradmagnthreshold, 'same');
% [linepar, line_edges,acc] = houghline(curves', magnitude, nrho, ntheta, gradmagnthreshold, nlines, verbose);
[linepar, acc] = houghline(curves, magnitude, nrho, ntheta, gradmagnthreshold, nlines, verbose);
figure
overlaycurves(pic, curves)
fprintf('Size of curves');
disp(size(curves))
fprintf('Size of line_edgesline_edgesline_edgesline_edges');
figure 
overlaycurves(pic, linepar)
end

function [linepar, line_edges, Accumulator] = houghline(curves, magnitude, nrho, ntheta, threshold, nlines, verbose)
fprintf('Size of magnitude');
disp(size(magnitude))

%%%%% HERE
disp(curves)
[M,N] = size(magnitude);
% [x_r,y_r] = 0.5 *[M,N];
x = 0.5 * [M N];
x_r = x(1);
y_r = x(2);
d_theta = pi / ntheta;
d_rho = sqrt(M^2 + N^2) / nrho ;
j_zero = nrho / 2;
% Accumulator = zeros(ntheta,nrho);
Accumulator = zeros(nrho,ntheta);
for k = 1:length(curves)
    uv = curves(k,:);
    r_xy = curves(k,:);
    uv = uv - [x_r, y_r];
    for i = 1:ntheta
        theta = d_theta * i;
        r = uv(1)*cos(theta) + uv(2)*sin(theta);
        j = j_zero + round(r / d_rho, 0) + 1;
%         if (round(r_xy(1), 0) <= M) && (round(r_xy(2), 0) <= N)
% %             Accumulator(i,j) = Accumulator(i,j) + magnitude(round(r_xy(1), 0),round(r_xy(2), 0));
%             Accumulator(j,i) = Accumulator(j,i) + magnitude(round(r_xy(1), 0),round(r_xy(2), 0));
%         end
  
        if (j >= 1) && (j <= nrho)
%             Accumulator(i,j) = Accumulator(i,j) + 1;
            Accumulator(j,i) = Accumulator(j,i) + 1;
        end
    end
end


[pos, value] = locmax8(Accumulator);
[dummy, indexvector] = sort(value);
nmaxima = size(value, 1);
for idx = 1:nlines
    rhoidxacc = pos(indexvector(nmaxima - idx + 1), 1);
    thetaidxacc = pos(indexvector(nmaxima - idx + 1), 2);
    angles = thetaidxacc * d_theta;
    rho = (rhoidxacc - j_zero)*d_rho;
    linepar(:,idx) = [rho; angles];
    
    x0 = 0;
    y0 = (rho - x0 * cos(angles))./sin(angles);
    dx = d_rho.^2;
    dy = (rho - dx * cos(angles))./sin(angles);
    %Given in lab description - visualizing results
    outcurves(1, 4*(idx-1) + 1) = 0;          
    outcurves(2, 4*(idx-1) + 1) = 3;
    outcurves(2, 4*(idx-1) + 2) = x0 - dx;
    outcurves(1, 4*(idx-1) + 2) = y0 - dy;
    outcurves(2, 4*(idx-1) + 3) = x0;
    outcurves(1, 4*(idx-1) + 3) = y0;
    outcurves(2, 4*(idx-1) + 4) = x0+dx;
    outcurves(1, 4*(idx-1) + 4) = y0+dy;
    
    line_edges = lines(rho, angles,5,max(M,N), 400);
end
linepar = outcurves;
% Accumulator
figure
imagesc(Accumulator);
% imagesc(Accumulator);


end

function [curves, Lv] =  extractedgeHough(inpic, scale, threshold, shape)

[pixels, Lv] = Lv1(discgaussfft(inpic, scale));
magnitude = (Lv  > threshold ) - 0.5; 
Lvvv = LvvvTilder(discgaussfft(inpic, scale), shape);
Lvv = LvvTilder(discgaussfft(inpic, scale), shape);
Lvvv_mask = (Lvvv < 0) - 0.5;
curves = zerocrosscurves(Lvv, Lvvv_mask);
fprintf('Extracted curve size:');
disp(size(curves))
curves = thresholdcurves(curves, magnitude);

end 

function [] = Question1()
[dxtools, dytools] = difference_operators();
end 

function [] = Question2()
    inpic = godthem256;
    showgrey(inpic)
    shape = 'same';
    threshold = 10;
    scale = 1;
    K = 9;
    inpic = discgaussfft(inpic, 1);
    [pixels, gradmagntools] = Lv(inpic, shape, threshold, K);
    figure;
    showgrey(gradmagntools);
    title("Magnitude ")
    figure;
    [pixels, gradmagntools] = Lv(discgaussfft(inpic, scale), shape, threshold, K);
    figure;
    showgrey(gradmagntools);
    title("Magnitude Gaussian filter ")
end 

function [] = Question3()
    img = few256;
    img = godthem256;
    scale = 64;
    threshold = 4;
    subplot(2, 3, 1)
    showgrey(img);
    r = [0.0001, 1.0, 4.0, 16.0,64.0];
     
    for i = 1:5
        subplot(2, 3, i+1)
        contour(LvvTilder(discgaussfft(img, r(i) ), 'same'), [0 0])
        axis('image')
        axis('ij')
        title(sprintf('scale = %d',r(i) ))
    end
end 

function [] = Question5()
    img = few256;
%     img = godthem256;
    subplot(2, 3, 1)
    showgrey(img);
    r = [0.0001, 1.0, 4.0, 16.0,64.0];
     
    for i = 1:5
        subplot(2, 3, i+1)
        contour(LvvTilder(discgaussfft(img, r(i) ), 'same'), [0 0])
        axis('image')
        axis('ij')
        title(sprintf('scale = %d',r(i) ))
    end
    
    figure
    subplot(2, 3, 1)
    showgrey(img);
    for i = 1:5
        subplot(2, 3, i + 1)
        showgrey(LvvvTilder(discgaussfft(img, r(i)), 'same') < 0)
        title(sprintf('scale = %d',r(i) ))
    end
    
end 
function [] = Question7_best(type)
img = few256;
img = godthem256;
r = [0.0001, 1.0, 4.0, 16.0,64.0];
K = 9;
N = 10;
[pixels, gradmagntools] = Lv1(img);
minimum = min(gradmagntools(:));
maximum = max(gradmagntools(:)) / K;
threshold = linspace(minimum,maximum,N);

% best_index = [9,13,14,18,33,42];
best_index = [9,14,18,33];

if type(1) == 'h'
    best_index = [23,13,18,32];
else
    best_index = [9,14,18,33];
end
for i = 1:length(best_index)
    row = floor(best_index(i)/length(r));
    col = best_index(i) - row*length(r);
    row = row + 1;
    subplot(2, length(best_index)/2, i)
%     figure;
    extractedge(img, r(col), threshold(row), 'same');
end

end 
function [] = Question7()

img = few256;
img = godthem256;
r = [0.0001, 1.0, 4.0, 16.0,64.0];
K = 9;
N = 10;
[pixels, gradmagntools] = Lv1(img);
minimum = min(gradmagntools(:));
maximum = max(gradmagntools(:)) / K;
threshold = linspace(minimum,maximum,N);
count = 0;
for i = 1:N
    for j = 1:length(r)
        count = count + 1;
        figure;
        extractedge(img, r(j), threshold(i), 'same');
    end
end
disp(count)


end 


function [curves] =  extractedge(inpic, scale, threshold, shape)

[pixels, Lv] = Lv1(discgaussfft(inpic, scale));
magnitude = (Lv  > threshold ) - 0.5; 
Lvvv = LvvvTilder(discgaussfft(inpic, scale), shape);
Lvv = LvvTilder(discgaussfft(inpic, scale), shape);
Lvvv_mask = (Lvvv < 0) - 0.5;
curves = zerocrosscurves(Lvv, Lvvv_mask);
curves = thresholdcurves(curves, magnitude);
overlaycurves(inpic, curves)

end 


function [gradmagntools] = point_wise_thresholding(dxtoolsconv,dytoolsconv,  threshold)
gradmagntools = sqrt(dxtoolsconv .^2 + dytoolsconv .^2);
showgrey((gradmagntools - threshold) > 0)
end

function [pixels_vv] = LvvvTilder(inpic, shape)
dxmask = [0 0 0 0 0;
          0 0 0 0 0;
          0 -0.5 0 0.5 0;
          0 0 0 0 0;
          0 0 0 0 0];

dymask = dxmask';

dxxmask = [0 0 0 0 0;
          0 0 0 0 0;
          0 1 -2 1 0;
          0 0 0 0 0;
          0 0 0 0 0];
dyymask = dxxmask';

dxymask = conv2(dxmask, dymask, 'same');
dxxxmask = conv2(dxmask, dxxmask, 'same');
dxxymask = conv2(dxxmask, dymask, 'same');
dxyymask = conv2(dxmask, dyymask, 'same');
dyyymask = conv2(dymask, dyymask, 'same');

Lx = filter2(dxmask, inpic, shape);
Ly = filter2(dymask, inpic, shape);

Lxx = filter2(dxxmask, inpic, shape);
Lyy = filter2(dyymask, inpic, shape);
Lxy = filter2(dxymask, inpic, shape);

Lxxx = filter2(dxxxmask, inpic, shape);
Lxxy = filter2(dxxymask, inpic, shape);
Lxyy = filter2(dxyymask, inpic, shape); 
Lyyy = filter2(dyyymask, inpic, shape); 
pixels_vv = (Lx.^3).*Lxxx + 3*(Lx.^2).*(Ly).*(Lxxy) + 3*(Lx).*(Ly.^2).*(Lxyy) + (Ly.^3).*Lyyy;
end

function [pixels_vv] = LvvTilder(inpic, shape)

% dxmask = [0 0 0; -0.5 0 0.5; 0 0 0];
% dymask = [0 -0.5 0; 0 0 0; 0 0.5 0];
% dxxmask = [0 0 0; 1 -2 1; 0 0 0];
% dyymask = [0 1 0; 0 -2 0; 0 1 0];

dxmask = [0 0 0 0 0;
          0 0 0 0 0;
          0 -0.5 0 0.5 0;
          0 0 0 0 0;
          0 0 0 0 0];

dymask = dxmask';

dxxmask = [0 0 0 0 0;
          0 0 0 0 0;
          0 1 -2 1 0;
          0 0 0 0 0;
          0 0 0 0 0];
dyymask = dxxmask';

Lx = filter2(dxmask, inpic, shape);
Ly = filter2(dymask, inpic, shape);


Lxx = filter2(dxxmask, inpic, shape);
Lyy = filter2(dyymask, inpic, shape);


dxymask = conv2(dxmask, dymask, 'same');
dxxxmask = conv2(dxmask, dxxmask, 'same');
dxxymask = conv2(dxxmask, dymask, 'same');

Lxy = filter2(dxymask, inpic, shape);

pixels_vv = (Lx.^2).*Lxx + 2*(Lx).*(Ly).*(Lxy) + (Ly.^2).*Lyy;

end

function [] = LvvTilde(N )

Lx = [0 0 0; -0.5 0 0.5; 0 0 0];
Ly = [0 -0.5 0; 0 0 0; 0 0.5 0];

Lxx = [0 0 0; 1 -2 1; 0 0 0];
Lyy = [0 1 0; 0 -2 0; 0 1 0];
Lxy = conv2(Lx, Ly, 'same');

dxxxmask = conv2(Lx, Lxx, 'same');
dxxymask = conv2(Lxx, Ly, 'same');

%N = 5;
[x, y] = meshgrid(-N:N, -N:N);

op1 = conv2(dxxxmask, x.^3, 'valid')
op2 = conv2(Lxx, x.^3, 'valid')
op3 = conv2(dxxymask, (x.^2).* y, 'valid')

figure;
showgrey(op1);
figure;
showgrey(op2);
figure;
showgrey(op3);


end

function [pixels, gradmagntools] = Lv1(inpic)

dymask = fspecial('sobel');
dxmask = dymask';
shape = 'same';
Lx = filter2(dxmask, inpic, shape);
Ly = filter2(dymask, inpic, shape);
pixels = Lx.^2 + Ly.^2;
gradmagntools = sqrt(pixels);

end 

function [pixels, gradmagntools] = Lv(inpic, shape, N,K)
subplot(2, N/2, 1)
showgrey(inpic)
if (nargin < 2)
    shape = 'same';
end
[dxmask, dymask] = sobels();
Lx = filter2(dxmask, inpic, shape);
Ly = filter2(dymask, inpic, shape);
pixels = Lx.^2 + Ly.^2;
gradmagntools = sqrt(pixels);
minimum = min(gradmagntools(:));
maximum = max(gradmagntools(:)) / K;
fprintf('Min %d \n',min(gradmagntools(:)));
fprintf('Max %d \n',max(gradmagntools(:)));
threshold = linspace(minimum,maximum,N);
for i=2:N
%     figure;
%     showgrey((gradmagntools - threshold(i)) > 0)
    subplot(2, N/2, i)
    showgrey((gradmagntools - threshold(i)) > 0)
    title(sprintf('threshold = %d', round(threshold(i))))
end 

end 

function [dxtools, dytools] = difference_operators()

tools = few256;
% subplot(2, N, i)
% showgrey(img)

% figure;
subplot(1, 3, 1)
showgrey(tools)
title('Original image')
[deltax, deltay] = sobels();
dxtools = conv2(tools, deltax, 'valid');
dytools = conv2(tools, deltay, 'valid');
% figure;
subplot(1, 3, 2)
showgrey(dxtools)
title('Sobel.X')
% figure;
subplot(1, 3, 3)
showgrey(dytools)
title('Sobel.Y')
fprintf('Size of Image');
disp(size(tools))
fprintf('Size of convolution.x ');
disp(size(dxtools))
fprintf('Size of convolution.x ');
disp(size(dytools))
end


function [s_x, s_y] = sobels()
s_temp = [-1 0 1; -2 0 2; -1 0 1];
s_x = (1/32)*[-3 0 3; -10 0 10; -3 0 3];
s_y = s_x';
end 

function [linepar, acc] = testLast()
ntheta = 256;
nrho = 256;
gradmagnthresh = 10;
verbose = 0;
nlines = 4;
scale = 4;
img = triangle128;
img = houghtest256;
% img = godthem256;
% img = few256;
[linepar, acc] = houghedgeline(img, scale, gradmagnthresh, nrho, ntheta, nlines, verbose);

end

function [linepar, Accumulator] = houghedgeline(img, scale, gradmagnthresh, nrho, ntheta, nlines, verbose)
curves = extractedge(img, scale, gradmagnthresh, 'same');
magnitude = Lv(img, 'same');

[linepar, Accumulator] = houghline(curves, magnitude, nrho, ntheta, gradmagnthresh, nlines, verbose);
                       
figure
overlaycurves(img, linepar);
axis([1 size(img, 2) 1 size(img, 1)]);                        
title('Result image')

figure
imagesc(Accumulator);
title('Hough Accumulator')


end


function [linepar, Accumulator] = houghline(curves, magn, nrho, ntheta, threshold, nlines, verbose)

Accumulator = zeros(nrho, ntheta);
theta_angles = linspace(-pi/2, pi/2, ntheta);
Divide_plan = sqrt(size(magn, 1).^2 + size(magn, 2).^2);
rhoy = linspace(-Divide_plan, Divide_plan, nrho);
insize = size(curves, 2) ;
ttrypointer = 1; 
while ttrypointer < insize
    polylength = curves(2, ttrypointer);
    disp("polylength : ")
    disp(polylength)
    ttrypointer = ttrypointer + 1;
    for curveidx = 1:polylength 
        x_coord = curves(2, ttrypointer);
        y_coord = curves(1, ttrypointer);
        ttrypointer = ttrypointer + 1;
        magn_xy = abs(magn(round(x_coord), round(y_coord)));
        if magn_xy > threshold
            for theta_index = 1:ntheta
                rhos = x_coord*cos(theta_angles(theta_index)) + y_coord*sin(theta_angles(theta_index));
                j = find(rhoy < rhos, 1, 'last');
%                 Accumulator(j, theta_index) = Accumulator(j, theta_index) + 1;
%                 Accumulator(j, theta_index) = Accumulator(j, theta_index) + magn_xy;
                Accumulator(j, theta_index) = Accumulator(j, theta_index) + log(magn_xy);
            end
        end
    end
end

[pos, value] = locmax8(Accumulator);
[dummy, indexvector] = sort(value);
nmaxima = size(value, 1);
disp("Rho and theta : ")
for index = 1:nlines
    rho_index_accumulator = pos(indexvector(nmaxima - index + 1), 1);
    theta_index_accumulator = pos(indexvector(nmaxima - index + 1), 2);
    rho = rhoy(rho_index_accumulator);
    theta = theta_angles(theta_index_accumulator);
    disp([rho, theta])
    linepar(:,index) = [rho; theta];
    
    x0 = 0;
    y0 = (rho - x0 * cos(theta))./sin(theta);
    dx = Divide_plan.^2;
    dy = (rho - dx * cos(theta))./sin(theta);

    outcurves(1, 4*(index-1) + 1) = 0;          
    outcurves(2, 4*(index-1) + 1) = 3;
    outcurves(2, 4*(index-1) + 2) = x0 - dx;
    outcurves(1, 4*(index-1) + 2) = y0 - dy;
    outcurves(2, 4*(index-1) + 3) = x0;
    outcurves(1, 4*(index-1) + 3) = y0;
    outcurves(2, 4*(index-1) + 4) = x0+dx;
    outcurves(1, 4*(index-1) + 4) = y0+dy;
end

linepar = outcurves;
end


function [linepar, Accumulator] = houghline1(curves, magn, nrho, ntheta, threshold, nlines, verbose)

Accumulator = zeros(nrho, ntheta);
theta_angles = linspace(-pi/2, pi/2, ntheta);
Divide_plan = sqrt(size(magn, 1).^2 + size(magn, 2).^2);
rhoy = linspace(-Divide_plan, Divide_plan, nrho);
insize = size(curves, 2) ;


[M,N] = size(magn);
x = 0.5 * [M N];
x_r = x(1);
y_r = x(2);
d_theta = pi / ntheta;
d_rho = sqrt(M^2 + N^2) / nrho ;
j_zero = nrho / 2;


ttrypointer = 1; 
while ttrypointer < insize
    polylength = curves(2, ttrypointer);
    ttrypointer = ttrypointer + 1;
    for curveidx = 1:polylength 
        x_coord = curves(2, ttrypointer);
        y_coord = curves(1, ttrypointer);
        ttrypointer = ttrypointer + 1;
        magn_xy = abs(magn(round(x_coord), round(y_coord)));
        %
        uv = [x_coord, y_coord] - [x_r, y_r];
        %
        if magn_xy > threshold
            for theta_index = 1:ntheta
                %
                theta = d_theta * theta_index;
                r = uv(1)*cos(theta) + uv(2)*sin(theta);
                j = j_zero + round(r / d_rho, 0) + 1;
                %
%                 rho_val = x_coord*cos(theta_angles(theta_index)) + y_coord*sin(theta_angles(theta_index));
%                 rho_index = find(rhoy < rho_val, 1, 'last');
%                 Accumulator(rho_index, theta_index) = Accumulator(rho_index, theta_index) + 1;
                 Accumulator(j, theta_index) = Accumulator(j, theta_index) + 1;
            end
        end
    end
end

[pos, value] = locmax8(Accumulator);
[dummy, indexvector] = sort(value);
nmaxima = size(value, 1);

for index = 1:nlines
    rho_index_accumulator = pos(indexvector(nmaxima - index + 1), 1);
    theta_index_accumulator = pos(indexvector(nmaxima - index + 1), 2);
    rho = rhoy(rho_index_accumulator);
    theta = theta_angles(theta_index_accumulator);
    linepar(:,index) = [rho; theta];
    
    x0 = 0;
    y0 = (rho - x0 * cos(theta))./sin(theta);
    dx = Divide_plan.^2;
    dy = (rho - dx * cos(theta))./sin(theta);

    outcurves(1, 4*(index-1) + 1) = 0;          
    outcurves(2, 4*(index-1) + 1) = 3;
    outcurves(2, 4*(index-1) + 2) = x0 - dx;
    outcurves(1, 4*(index-1) + 2) = y0 - dy;
    outcurves(2, 4*(index-1) + 3) = x0;
    outcurves(1, 4*(index-1) + 3) = y0;
    outcurves(2, 4*(index-1) + 4) = x0+dx;
    outcurves(1, 4*(index-1) + 4) = y0+dy;
end

linepar = outcurves;
end




function pixels = Lv(inpic, shape)
	if(nargin<2)
		shape = 'same';
	end;

	dymask = fspecial('sobel');
	dxmask = dymask';
	Lx = filter2(dxmask, inpic, shape);
	Ly = filter2(dymask, inpic, shape);
	pixels = sqrt(Lx.^2 + Ly.^2);
end

function edgecurves = extractedge(inpic, scale, threshold, shape)

	Lv_init = Lv(discgaussfft(inpic, scale), shape);
	Lvv = Lvvtilde(discgaussfft(inpic, scale), shape);
	Lvvv = Lvvvtilde(discgaussfft(inpic, scale), shape);
	Lv_mask = (Lv_init > threshold);
	Lvvv_mask = (Lvvv < 0);
	edgecurves = zerocrosscurves(Lvv, Lvvv_mask);
	edgecurves = thresholdcurves(edgecurves, Lv_mask);

end


function pixels = Lvvtilde(inpic, shape)
	dxmask = [0 0 0 0 0;
			  0 0 0 0 0;
			  0 -1/2 0 1/2 0;
			  0 0 0 0 0;
			  0 0 0 0 0];

	dymask = dxmask';

	dxxmask = [0 0 0 0 0;
				0 0 0 0 0;
				0 1 -2 1 0;
				0 0 0 0 0;
				0 0 0 0 0];
	dyymask = dxxmask';

	dxymask = conv2(dxmask, dymask, shape);

	Lx = filter2(dxmask, inpic, shape);
	Ly = filter2(dymask, inpic, shape);
	Lxx = filter2(dxxmask, inpic, shape);
	Lxy = filter2(dxymask, inpic, shape);
	Lyy = filter2(dyymask, inpic, shape);

	pixels = Lx.^2.*Lxx + 2*Lx.*Lxy.*Ly + Ly.^2.*Lyy;

end


function pixels = Lvvvtilde(inpic, shape)

	dxmask = [0 0 0 0 0;
			  0 0 0 0 0;
			  0 -1/2 0 1/2 0;
			  0 0 0 0 0;
			  0 0 0 0 0];
	dymask = dxmask';

	dxxmask = [0 0 0 0 0;
				0 0 0 0 0;
				0 1 -2 1 0;
				0 0 0 0 0;
				0 0 0 0 0];
	dyymask = dxxmask';
	dxymask = conv2(dxmask, dymask, shape);

	dxxxmask = conv2(dxmask, dxxmask, shape);
	dxyymask = conv2(dxmask, dyymask, shape);
	dxxymask = conv2(dxxmask, dymask, shape);
	dyyymask = conv2(dymask, dyymask, shape);

	Lx = filter2(dxmask, inpic, shape);
	Ly = filter2(dymask, inpic, shape);

	Lxxx = filter2(dxxxmask, inpic, shape);
	Lxyy = filter2(dxyymask, inpic, shape);
	Lxxy = filter2(dxxymask, inpic, shape);
	Lyyy = filter2(dyyymask, inpic, shape);

	pixels = Lx.^3.*Lxxx + 3.*Lx.^2.*Ly.*Lxxy + 3.*Lx.*Ly.^2.*Lxyy + Ly.^3.*Lyyy;

end
