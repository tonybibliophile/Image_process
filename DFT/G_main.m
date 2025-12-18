clear all; close all; clc;

%% 1. 準備影像與雜訊
try
    img = imread('Zebra.jpg');
catch
    img = imread('peppers.png'); 
end
img = im2double(img);
if size(img, 3) == 3
    img_gray = rgb2gray(img); % 分析頻譜時通常看灰階比較清楚
else
    img_gray = img;
end

% 加入高斯雜訊
sigma = 0.05;
noisy_img = img_gray + sigma * randn(size(img_gray));

%% 2. 計算頻譜 (Fourier Transforms)
% A. 原圖的頻譜
F_orig = fftshift(fft2(img_gray));
S_orig = log(1 + abs(F_orig)); % 取 log 方便觀察

% B. 雜訊圖的頻譜
F_noisy = fftshift(fft2(noisy_img));
S_noisy = log(1 + abs(F_noisy));

%% 3. 建立高斯濾波器 (Gaussian Filter)
[M, N] = size(img_gray);
D0 = 80; % 截止頻率

u = -floor(M/2) : floor((M-1)/2);
v = -floor(N/2) : floor((N-1)/2);
[V, U] = meshgrid(v, u);
D_sq = U.^2 + V.^2;

% 高斯公式
H = exp(-D_sq / (2 * D0^2));

%% 4. 濾波後的頻譜
F_filtered = F_noisy .* H;
S_filtered = log(1 + abs(F_filtered));

%% 5. 繪製「講義同款」四格圖
figure('Name', 'Frequency Domain Analysis', 'Position', [100, 100, 800, 800]);

% (a) 原圖頻譜
subplot(2, 2, 1); 
imshow(mat2gray(S_orig)); 
title('(a) log FT original image');
xlabel('原本的訊號集中在中心');

% (b) 雜訊頻譜
subplot(2, 2, 2); 
imshow(mat2gray(S_noisy)); 
title('(b) log FT noisy image');
xlabel('高斯雜訊分布在整張頻譜(灰底)');

% (c) 濾波器 (注意：這次是圓形的！)
subplot(2, 2, 3); 
imshow(H, []); 
title('(c) Gaussian Filter');
xlabel('高斯濾波器 (中間亮，向外柔和變暗)');

% (d) 濾波後頻譜
subplot(2, 2, 4); 
imshow(mat2gray(S_filtered)); 
title('(d) after filter FT image');
xlabel('雜訊被壓暗了，只剩中心訊號');