clear all;
close all;
addpath(genpath(pwd))

%% Load image  and add noise
% load clean image
img = imread('lena.jpg');
img = im2double(img);

% filter size
r = 23;
% add noise
noisy = add_noise(img,r); 

figure();
subplot(1,2,1),imshow(img);title('original image')
subplot(1,2,2),imshow(noisy);title('noisy image')

%% Perform FT on the original image
% implement FT below
img_fq = FT(img);
% center shift (You can read Section: Periodicity and Center shifting)
img_fq = fftshift(img_fq); 
% for displaying the power part by log scale
log_img_fq = log(1+abs(img_fq)); 
% normalize the power to 0-1
log_img_fq = mat2gray(log_img_fq); 

figure();
subplot(2,2,1),imshow(log_img_fq);title('log FT original image')

%% Perform FT on the noisy image
% implement FT below
noisy_fq = FT(noisy);
noisy_fq = fftshift(noisy_fq);
log_noisy_fq = log(1+abs(noisy_fq));
log_noisy_fq = mat2gray(log_noisy_fq);

subplot(2,2,2),imshow(log_noisy_fq);title('log FT noisy image')

%% Creat a filter
% Try different filters for extra credit
% get the size of the input image
[m, n, z] = size(img); 
% create a rectangular filter at center
filter = zeros(m,n);
filter(r:m-r,r:n-r) =1;

subplot(2,2,3),imshow(filter,[]);title('filter')

%% Filter out the noise
fil_img = noisy_fq .* filter; 

log_fil_img = log(1+abs(fil_img));
log_fil_img = mat2gray(log_fil_img);

subplot(2,2,4),imshow(log_fil_img);title('FT image after filter')

%% Inverse Fourier transform
% unshift
fil_img = ifftshift(fil_img); 
% implement IFT below
result = IFT(fil_img); 
result = mat2gray(real(result));
figure();
subplot(1,2,1),imshow(noisy);title('noisy image')
subplot(1,2,2);imshow(result,[]);title('denoised image')

%% Implement your FT/IFT function here
% Hint:
% Images are 3D tensors (height × width × channel).  
% Implement FT and IFT by applying the transform on each channel independently.
%
% Note:
% Do NOT include fftshift or ifftshift inside your FT/IFT functions.
% Calling built-in fft/ifft/fft2/ifft2 functions will not be credited.
% You are expected to implement the Fourier transform manually using the mathematical definition.

function [I_freq] = FT(I)
    [m,n,c] = size(I);
    temp = zeros(m,n,c);
    
    for z = 1:c
        for u = 1:m
            for v = 1:n
                l = 0:m-1;
                j = 0:n-1;
                [L,J] = meshgrid(l,j);
                E = exp(-1i*2*pi*(L*(u-1)/m + J*(v-1)/n));
                temp(u,v,z) = sum(sum(I(:,:,z).*E));
            end
        end
    end
    
    I_freq = temp;
end

function [I] = IFT(I_freq)
    [m,n,c] = size(I_freq);
    temp = zeros(m,n,c);
    
    for z = 1:c
        for u = 1:m
            for v = 1:n
                l = 0:m-1;
                j = 0:n-1;
                [L,J] = meshgrid(l,j);
                E = exp(1i*2*pi*(L*(u-1)/m + J*(v-1)/n));
                temp(u,v,z) = sum(sum(I_freq(:,:,z).*E))/(m*n);
            end
        end
    end
    
    I = temp;
end

