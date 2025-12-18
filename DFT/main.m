clear all;
close all;
addpath(genpath(pwd))

%% Load image and add noise
% load clean image
% 確保你的資料夾內有 lena.jpg
img = imread('lena.jpg'); 

% 如果是彩圖，某些步驟可能需要處理單通道，但這裡保持原樣
img = im2double(img);

% filter size
r = 20;

% add noise
% 注意：你必須確保目錄下有 add_noise.m 檔案
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

%% Create a filter
% Try different filters for extra credit
% get the size of the input image
[m, n, z] = size(img); 

% create a rectangular filter at center
% 注意：如果是 RGB 彩圖 (z=3)，濾波器需要複製到三層，或者在運算時廣播
filter = zeros(m,n);
filter(r:m-r,r:n-r) = 1;

% 顯示濾波器 (如果是多通道，只顯示第一層以免報錯)
subplot(2,2,3),imshow(filter,[]);title('filter')

%% Filter out the noise
% 如果 noisy_fq 是 3D 而 filter 是 2D，Matlab 新版會自動廣播，舊版需注意
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
function [I_freq] = FT(I)
    % 獲取影像尺寸：高度 M，寬度 N，通道數 C (若是 RGB 則 C=3)
    [M, N, C] = size(I);
    I_freq = zeros(M, N, C);
    
    % --- 準備 DFT 矩陣 (避免多層迴圈，使用矩陣乘法加速) ---
    % 1. 針對「高度 M」的 DFT 矩陣 (用於處理 Column)
    m_idx = 0:M-1;
    k_m = m_idx';
    % W_M 是一個 M x M 的矩陣
    W_M = exp(-1j * 2 * pi * k_m * m_idx / M);
    
    % 2. 針對「寬度 N」的 DFT 矩陣 (用於處理 Row)
    n_idx = 0:N-1;
    k_n = n_idx';
    % W_N 是一個 N x N 的矩陣
    W_N = exp(-1j * 2 * pi * k_n * n_idx / N);
    
    % --- 開始轉換 ---
    % 針對每個顏色通道獨立處理
    for c = 1:C
        channel_data = double(I(:,:,c));
        
        % 數學原理：F = W_M * Image * W_N.'
        % 第一步：W_M * channel_data (對每一行做 DFT)
        % 第二步：結果 * W_N.' (對每一列做 DFT)
        I_freq(:,:,c) = W_M * channel_data * W_N.';
    end
end

function [I] = IFT(I_freq)
    [M, N, C] = size(I_freq);
    I = zeros(M, N, C);
    
    % --- 準備 IDFT 矩陣 ---
    % IDFT 的指數是正的 (exp(j...))
    
    % 1. 高度 M 的 IDFT 矩陣
    m_idx = 0:M-1;
    k_m = m_idx';
    W_M_inv = exp(1j * 2 * pi * k_m * m_idx / M);
    
    % 2. 寬度 N 的 IDFT 矩陣
    n_idx = 0:N-1;
    k_n = n_idx';
    W_N_inv = exp(1j * 2 * pi * k_n * n_idx / N);
    
    % --- 開始反轉換 ---
    for c = 1:C
        channel_freq = I_freq(:,:,c);
        
        % 數學原理：Image = (W_M_inv * F * W_N_inv.') / (M*N)
        temp = W_M_inv * channel_freq * W_N_inv.';
        
        % IDFT 記得要除以總像素數 (Normalization)
        I(:,:,c) = temp / (M * N);
    end
end