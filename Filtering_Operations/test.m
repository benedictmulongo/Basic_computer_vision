%addpath lab_library/;
% Add the directory and its subdirectories 
%addpath(genpath('lab_library/'))

function [] = test(p,q)

    smoothingAndSubsamplingQuestion17_19('gaussian filter')
    
    smoothingAndSubsamplingQuestion17_19('IdealLowPass filter')
    %Question17_gaussFilter_gauss()
    %Question17_gaussFilter_salt()
    %Question17_MedianFilt_salt()
    %Question17_MedianFilt_gauss()
    %Question17_idealLowPass_salt()
    %Question17_idealLowPass_gauss()
    %Question17()
    %Question16()
    %Question15_specific()
    %Question14_15()
    %Question13()
    %rotation12()
    %scaling11()
    %Question10()
    %linearityQuestion9()
    %linearity()
    %question2(p,q)
    %question1()
    %fftwave(p, q, 128)
    %imageCenterFft(p,q)
end 

function [] = smoothingAndSubsamplingQuestion17_19(type)
figure
img = phonecalc256;
smoothimg = img;
N=5;
k = [1.0,4.0,8,16,32];
t = [0.6, 0.3, 0.2 ,1.0, 0.09];
for i=1:N
    if i>1 % generate subsampled versions
        img = rawsubsample(img);
        
        if type(1) == 'g'
            smoothimg = gaussfft(smoothimg, k(i));
        else
            smoothimg = ideal( smoothimg, t(i), 'l');
        end
        %smoothimg = % <call_your_filter_here>(smoothimg, <params>);
        smoothimg = rawsubsample(smoothimg);
    end
    subplot(2, N, i)
    showgrey(img)
    subplot(2, N, i+N)
    showgrey(smoothimg)
    
end

title(type)
end 

function [] = Question17_gaussFilter_gauss()

office = office256;
add = gaussnoise(office, 16);

figure
showgrey(office);
title('Original image ')
figure
showgrey(add);
title('Image gaussian noise')

k = [1.0,4.0,16.0,64];
for i = 1:length(k)
    [filtered, inverted_pic, gauss] = gaussfft(add, k(i));
    figure
    showgrey(filtered);
    title(sprintf('Image GaussNoise sigma = %d', k(i)))
end

end

function [] = Question17_gaussFilter_salt()
office = office256;
sap = sapnoise(office, 0.1, 255);

figure
showgrey(office);
title('Original image ')
figure
showgrey(sap);
title('Image with salt noise')


k = [1.0,4.0,16.0,64];
for i = 1:length(k)
    [filtered, inverted_pic, gauss] = gaussfft(sap, k(i));
    figure
    showgrey(filtered);
    title(sprintf('Image SaltPepper sigma = %d', k(i)))
end

end

function [] = Question17_MedianFilt_salt()
office = office256;
sap = sapnoise(office, 0.1, 255);

figure
showgrey(office);
title('Original image ')
figure
showgrey(sap);
title('Image with salt noise')

k = [2, 4, 6 , 9];
for i = 1:length(k)
    Filtim = medfilt( sap, k(i), k(i));
    figure
    showgrey(Filtim);
    title('Image gauss IdealFilt ') 
    title(sprintf('Image SaltPepper Median Arraysize = %d', k(i)))
end

end

function [] = Question17_MedianFilt_gauss()
office = office256;
add = gaussnoise(office, 16);

figure
showgrey(office);
title('Original image ')
figure
showgrey(add);
title('Image gaussian noise')

k = [2, 4, 6 , 9];
for i = 1:length(k)
    Filtim = medfilt( add, k(i), k(i));
    figure
    showgrey(Filtim);
    title('Image gauss IdealFilt ') 
    title(sprintf('Image gauss Median Arraysize = %d', k(i)))
end

end


function [] = Question17_idealLowPass_gauss()
office = office256;
add = gaussnoise(office, 16);

figure
showgrey(office);
title('Original image ')
figure
showgrey(add);
title('Image gaussian noise')

k = [0.6, 0.3, 0.2 , 0.09];
for i = 1:length(k)
   [ Filtim, MTF] = ideal( add, k(i), 'l');
    figure
    showgrey(Filtim);
    title('Image gauss IdealFilt ') 
    title(sprintf('Image gauss IdealFilt cut = %f', k(i)))
end

end 

function [] = Question17_idealLowPass_salt()
office = office256;
sap = sapnoise(office, 0.1, 255);

figure
showgrey(office);
title('Original image ')
figure
showgrey(sap);
title('Image with salt noise')

k = [0.6, 0.3, 0.2 , 0.09];
for i = 1:length(k)
   [ Filtim, MTF] = ideal( sap, k(i), 'l');
    figure
    showgrey(Filtim);
    title('Image gauss IdealFilt ') 
    title(sprintf('Image SaltPepper IdealFilt cut = %f', k(i)))
end

end 

function [] = Question17()
office = office256;
add = gaussnoise(office, 16);
sap = sapnoise(office, 0.1, 255);

figure
showgrey(office);
title('Original image ')
figure
showgrey(add);
title('Image gaussian noise')
figure
showgrey(sap);
title('Image with salt noise')

end

function [] = Question16()

    img1 = phonecalc128;
    img2 = genevepark128;
    img3 = sunflower128;
    img4 = suburb256;

    pic = deltafcn(128, 128);
    t = [1.0,4.0,16.0,64.0,256.0];
    
    figure
    showgrey(img1);
    title("Original image")
    
    for i = 1:5
        
        [filtered, inverted_pic, gauss] = gaussfft(img1, t(i));
        figure
        showgrey(filtered);
        title(sprintf('F filtered t = %d', t(i)))
    end
    
end 

function [] = Question15_specific()

    mu = [0 0];
    sigma = [1 0; 0 1];
    
    k = [0.1,0.3,1,10];
    for i = 1:4
        t = k(i);
        xdim = 128;
        ydim = 128;
        [x, y] = meshgrid(-xdim/2 : (xdim/2)-1, -ydim/2 : (ydim/2)-1);
        % The gaussian filter
        gaussian_kernel = 1/(2*pi*t)*exp(-(x.^2+y.^2)./(2*t)); 
        gaussian_fft = fft2(gaussian_kernel);
        gaussian_fft = fftshift(gaussian_fft);
        X = [x(:) x(:)];
        ideal_case = mvnpdf(X,mu,t * sigma);
        ideal_case = reshape(ideal_case,ydim,xdim);
        figure
        showgrey(ideal_case);
        title('ideal_case gaussian')
        figure
        showgrey(gaussian_kernel);
        title('Gaussian_kernel Spatial domain')

        %figure
        %showgrey(gaussian_fft);
        %title('Gaussian_fft freq. domain')
    end 
end 


function [] = Question14_15()
    pic = deltafcn(128, 128);
    figure
    showgrey(pic);
    title('Delta func')
    t = [0.1,0.3,1,10,100];
    for i = 1:5
        
        [filtered, inverted_pic, gauss] = gaussfft(pic, t(i));
        figure
        showgrey(filtered);
        title(sprintf('F filtered t = %d', t(i)))
        ideal = t(i).*eye(2);
        %var = variance(gauss)
        var = variance(filtered)
        %title('F filtered')
    end
end


function [] = Question13()

img1 = phonecalc128;
img2 = genevepark128;
img3 = few128;
figure
showgrey(img1);
title('phonecalc128')
a = 10^-10;
pixels = pow2image(img1, a);
figure
showgrey(pixels);
title('phonecalc128 rep. magn')


figure
showgrey(img3);
title('genevepark128')
a = 10^-10;
pixels = pow2image(img3, a);
figure
showgrey(pixels);
title('genevepark128 rep. magn')

figure
showgrey(img3);
title('few128')
pixels = randphaseimage(img3);
figure
showgrey(pixels);
title('few128 new phase')
pixels = randphaseimage(img1);
figure
showgrey(pixels);
title('phonecalc128 new phase')

end 

function [filtered, inverted_pic,gaussian_kernel] = gaussfft(pic,t)
    [xdim, ydim] = size(pic);
    [x, y] = meshgrid(-xdim/2 : (xdim/2)-1, -ydim/2 : (ydim/2)-1);
    % The gaussian filter
    gaussian_kernel = 1/(2*pi*t)*exp(-(x.^2+y.^2)./(2*t)); 
    
    %2. Fourier transform the original image and the Gaussian filter.
    pic_fft = fft2(pic);
    gaussian_fft = fft2(gaussian_kernel);
    
    %3. Multiply the Fourier transforms.
    convolved_image = (pic_fft .* gaussian_fft);
    
    %4. Inversion of Fourier transform.
    inverted_pic = ifft2(convolved_image);
    filtered = fftshift(inverted_pic);
end

function [] = rotation12()
% 128 X 128 test images 
F = [zeros(60, 128); ones(8, 128); zeros(60, 128)] .* [zeros(128, 48) ones(128, 32) zeros(128, 48)];

figure
showgrey(F);
title('F')

figure
F_fft = fft2(F);
%showgrey(G_fft);
showfs( F_fft);
title('Fft of F')

figure
G = rot(F, 30);
showgrey(G);
title('rot(F,30')

figure
G_fft = fft2(G);
%showgrey(G_fft);
showfs( G_fft);
title('Fft(rot(F,30)) ')


%%%
figure
Hhat = rot(fftshift(G_fft), -30);
%showgrey(Hhat);
showgrey(log(1 + abs(Hhat)));
title('rot(rot(F,30),-30)')


end 


function [] = scaling11()
% 128 X 128 test images 
F = [zeros(60, 128); ones(8, 128); zeros(60, 128)] .* [zeros(128, 48) ones(128, 32) zeros(128, 48)];

figure
showgrey(F);
title('F scaling')

Fhat = fft2(F);
figure
showgrey(log(1 + abs(fftshift(Fhat))));
title('Fhat')

figure
showfs( Fhat);
title('Fhat magn')

%%%
F = [ zeros(56, 128); ones(16, 128); zeros(56, 128)];
G = F';
figure 
showgrey(F .* G);
figure 
showfs(fft2(F .* G));
title('F before')

end 

function [] = Question10()
F = [ zeros(56, 128); ones(16, 128); zeros(56, 128)];
G = F';
figure
showgrey(F);
title('F')
figure
showgrey(G);
title('G = F.transpose')
% Element-wise multiplication 
figure 
showgrey(F .* G);
figure 
showfs(fft2(F .* G));

Cfull = conv2(fft2(F),fft2(G))/(128*128);
Cfull = Cfull(1:128,1:128);
figure
%showgrey(Cfull);
showgrey(log(1 + abs(fftshift(Cfull))));
title('conv2(fft2(F),fft2(G))/(128*128)')

end

function [] = linearityQuestion9()
% 128 X 128 test images 
F = [ zeros(56, 128); ones(16, 128); zeros(56, 128)];
G = F';
H = F + 2 * G;

Hhat_1 = fft2(H);
Hhat_2 = fft2(F) + 2*fft2(G);

figure
showgrey(log(1 + abs(fftshift(Hhat_1))));
title('fft2(F + 2 * G)')

figure
showgrey(log(1 + abs(fftshift(Hhat_2))));
title('fft2(F) + 2*fft2(G) ')

end 


function [] = linearity()
% 128 X 128 test images 
F = [ zeros(56, 128); ones(16, 128); zeros(56, 128)];
G = F';
H = F + 2 * G;
figure
showgrey(F);
title('F image')
figure
showgrey(G);
title('G image')
figure
showgrey(H);
title('H image')
% display them with showgrey. Then compute the discrete Fourier transform of the images 
Fhat = fft2(F);
Ghat = fft2(G);
Hhat = fft2(H);
Hhat_t = Fhat + 2*Ghat;
% Show their Fourier spectra 
figure
showgrey(log(1 + abs(fftshift(Fhat))));
title('Spectrum Fhat')
figure
showgrey(log(1 + abs(fftshift(Ghat))));
title('Spectrum Ghat')
figure
showgrey(log(1 + abs(fftshift(Hhat))));
title('Spectrum H_hat')
% Bis
figure
showgrey(log(1 + abs(fftshift(Hhat))));
title('Abs(fftshift(Hhat))')

figure
showgrey(log(1 + abs(fftshift(Hhat_t))));
title('Hhat_t frequency ')

end 

function [] = question2(p,q)

    Fhat = zeros(128,128);
    Fhat(p,q) = 1;
    F = ifft2(Fhat);
    
    figure
    showgrey(Fhat);
    title('F in frequency domain')
    figure
    showgrey(fftshift(Fhat));
    title('F shift at center')
    F_singlePRoj = singelPointProjectionWave(p,q, 128,128);
    figure
    showgrey(real(F_singlePRoj));
    title('F singel Point Projection Wave')
    figure
    [X,Y] = meshgrid(1:128,1:128);
    surf(X,Y,F_singlePRoj')
    title('F singel Point Projection Wave')
    
end 

function [X] = singelPointProjectionWave(p,q, m,n)
X = zeros(m,n);
for i = 1:m
    for j = 1:n
        X(i,j) = 1/(m*n)*( cos((2*pi/n) * (i*p + j*q) ) + i*sin((2*pi/n) * (i*p + j*q) ) ) ;
    end
end

end 

function [] = question1()
s = [[5,9];[9,5]; [17,9];[17,121];[5,1];[125,1] ]
[n,dim] = size(s);
for i = 1:n
    t_i = s(i, :);
    p = t_i(1);
    q = t_i(2);
    figure
    fftwave(p, q, 128)
end

end 


function [] = imageCenterFft(p,q)
    Fhat = zeros(128,128);
    Fhat(p,q) = 1;
    F = ifft2(Fhat);

    Fabsmax = max(abs(F(:)));
    figure
    showgrey(real(F), 64, -Fabsmax, Fabsmax)
    figure
    showgrey(imag(F), 64, -Fabsmax, Fabsmax)
    figure
    showgrey(abs(F), 64, -Fabsmax, Fabsmax)
    figure
    showgrey(angle(F), 64, -pi, pi)

end 

