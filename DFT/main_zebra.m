clear all;
close all;
addpath(genpath(pwd))

%% Load image and add noise
% load clean image
try
    img = imread('Zebra.jpg'); % 嘗試讀取 Zebra
catch
    img = imread('peppers.png'); % 備用圖
end
% 轉為 double 以便進行數學運算
img = im2double(img);

% parameter for Gaussian Filter (Cutoff Frequency)
D0 = 40; % 截止頻率，相當於以前的 r，控制模糊程度

% add Gaussian noise
% 直接使用 randn 加上常態分佈雜訊
sigma = 0.05; % 雜訊強度
noise = sigma * randn(size(img));
noisy = img + noise;
% 確保數值範圍在 0~1
noisy = max(0, min(noisy, 1));

figure();
subplot(1,2,1),imshow(img);title('original image')
subplot(1,2,2),imshow(noisy);title('noisy image (Gaussian)')

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
noisy_fq = fftshift(noisy_fq); % 注意：這裡做了 shift，直流分量移到了中心
log_noisy_fq = log(1+abs(noisy_fq));
log_noisy_fq = mat2gray(log_noisy_fq);
subplot(2,2,2),imshow(log_noisy_fq);title('log FT noisy image')

%% Create a Gaussian Filter
% get the size of the input image
[M, N, C] = size(img); 

% create coordinate grid (centered)
u = -floor(M/2) : floor((M-1)/2);
v = -floor(N/2) : floor((N-1)/2);
[V, U] = meshgrid(v, u);

% Calculate distance from center (D^2 = u^2 + v^2)
D_sq = U.^2 + V.^2;

% Gaussian Low Pass Filter Formula: H = exp(-D^2 / (2*D0^2))
filter = exp(-D_sq / (2 * D0^2));

% 顯示濾波器 (顯示 2D 遮罩即可)
subplot(2,2,3),imshow(filter,[]);title(['Gaussian Filter (D0=' num2str(D0) ')'])

%% Filter out the noise
% 進行頻域濾波
% Matlab 會自動將 2D 的 filter 廣播乘到 3D 的 noisy_fq (RGB) 上
fil_img = noisy_fq .* filter; 

log_fil_img = log(1+abs(fil_img));
log_fil_img = mat2gray(log_fil_img);
subplot(2,2,4),imshow(log_fil_img);title('FT image after Gaussian filter')

%% Inverse Fourier transform
% unshift (把中心點移回角落，準備做反轉換)
fil_img = ifftshift(fil_img); 

% implement IFT below
result = IFT(fil_img); 

% 取實部並正規化
result = real(result);
result = max(0, min(result, 1)); % Clip to 0-1

figure();
subplot(1,2,1),imshow(noisy);title('noisy image')
subplot(1,2,2);imshow(result,[]);title('denoised image (Gaussian)')

%% Implement your FT/IFT function here
function [I_freq] = FT(I)
    % 獲取影像尺寸：高度 M，寬度 N，通道數 C
    [M, N, C] = size(I);
    I_freq = zeros(M, N, C);
    
    % --- 準備 DFT 矩陣 ---
    % 1. 針對「高度 M」的 DFT 矩陣
    m_idx = 0:M-1;
    k_m = m_idx';
    W_M = exp(-1j * 2 * pi * k_m * m_idx / M);
    
    % 2. 針對「寬度 N」的 DFT 矩陣
    n_idx = 0:N-1;
    k_n = n_idx';
    W_N = exp(-1j * 2 * pi * k_n * n_idx / N);
    
    % --- 開始轉換 ---
    for c = 1:C
        % 數學原理：F = W_M * Image * W_N.'
        % 這裡保留你原本使用的轉置寫法 W_N.' (雖然 DFT 矩陣是對稱的，結果不變)
        I_freq(:,:,c) = W_M * double(I(:,:,c)) * W_N.';
    end
end

function [I] = IFT(I_freq)
    [M, N, C] = size(I_freq);
    I = zeros(M, N, C);
    
    % --- 準備 IDFT 矩陣 ---
    % 1. 高度 M 的 IDFT 矩陣 (指數為正)
    m_idx = 0:M-1;
    k_m = m_idx';
    W_M_inv = exp(1j * 2 * pi * k_m * m_idx / M);
    
    % 2. 寬度 N 的 IDFT 矩陣 (指數為正)
    n_idx = 0:N-1;
    k_n = n_idx';
    W_N_inv = exp(1j * 2 * pi * k_n * n_idx / N);
    
    % --- 開始反轉換 ---
    for c = 1:C
        channel_freq = I_freq(:,:,c);
        
        % 數學原理：Image = (W_M_inv * F * W_N_inv.') / (M*N)
        temp = W_M_inv * channel_freq * W_N_inv.';
        
        % IDFT 記得要除以總像素數
        I(:,:,c) = temp / (M * N);
    end
end