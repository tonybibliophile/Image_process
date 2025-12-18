clear all; close all; clc;

%% 0. 設定環境與讀取圖片
save_dir = 'Color_Zebra_Output';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
fprintf('處理結果將儲存至: %s\n', save_dir);

% 讀取圖片 (確保是彩色)
try
    img = imread('Zebra.jpg');
catch
    warning('找不到 Zebra.jpg，使用內建 peppers.png 代替');
    img = imread('peppers.png'); 
end
img = im2double(img); % 轉為 double (0~1)

% 檢查是否為彩色，若不是則轉為 RGB (為了演示彩色處理流程)
if size(img, 3) == 1
    img = cat(3, img, img, img);
end

imwrite(img, fullfile(save_dir, '0_Original.png'));

%% 1. 產生兩種彩色雜訊 (Color Noise)

% --- A. 彩色高斯雜訊 (Simulating Real Noise) ---
sigma = 0.05; % 雜訊強度
noisy_gauss = add_g_noise(img, sigma); % 調用下方的函式

% --- B. 彩色椒鹽雜訊 (Simulating Impulse Noise) ---
density = 0.05;
noisy_sp = imnoise(img, 'salt & pepper', density);

% 儲存雜訊圖
imwrite(noisy_gauss, fullfile(save_dir, 'A_Noisy_Gaussian.png'));
imwrite(noisy_sp,    fullfile(save_dir, 'B_Noisy_SaltPepper.png'));

%% 2. 建立濾波器 (Filters)
[M, N, C] = size(img);

% 建立頻率座標網格
u = 0:(M-1); v = 0:(N-1);
idx_u = find(u > M/2); u(idx_u) = u(idx_u) - M;
idx_v = find(v > N/2); v(idx_v) = v(idx_v) - N;
[V, U] = meshgrid(v, u);
D = sqrt(U.^2 + V.^2);
D = fftshift(D); % 移到中心

% 設定截止頻率 (Cutoff Frequency)
D0 = 30; % 設在 50，讓理想濾波器的振鈴稍微明顯一點

% --- 理想濾波器 (Ideal LPF - 講義款) ---
H_Ideal = double(D <= D0);

% --- 高斯濾波器 (Gaussian LPF) ---
H_Gauss = exp(-(D.^2) / (2 * D0^2));

%% 3. 執行彩色濾波 (Processing RGB Channels)

% 我們寫了一個 helper function (apply_filter_rgb) 在最下面
% 這樣程式碼比較乾淨，不用重複寫迴圈

% --- 第一組：對付高斯雜訊 ---
res_Gauss_IdealFilter = apply_filter_rgb(noisy_gauss, H_Ideal);
res_Gauss_GaussFilter = apply_filter_rgb(noisy_gauss, H_Gauss);

% --- 第二組：對付椒鹽雜訊 ---
res_SP_IdealFilter    = apply_filter_rgb(noisy_sp, H_Ideal);
res_SP_GaussFilter    = apply_filter_rgb(noisy_sp, H_Gauss);

% 儲存結果
imwrite(res_Gauss_IdealFilter, fullfile(save_dir, 'Result_GaussNoise_IdealFilter.png'));
imwrite(res_Gauss_GaussFilter, fullfile(save_dir, 'Result_GaussNoise_GaussFilter.png'));
imwrite(res_SP_IdealFilter,    fullfile(save_dir, 'Result_SPNoise_IdealFilter.png'));
imwrite(res_SP_GaussFilter,    fullfile(save_dir, 'Result_SPNoise_GaussFilter.png'));

%% 4. 繪圖比較 - 第一戰：高斯雜訊 (真實模擬)
figure('Name', 'Battle 1: Gaussian Noise', 'Position', [50, 100, 1200, 400]);

subplot(1, 3, 1); imshow(noisy_gauss); title('1. 原圖 + 高斯雜訊');
subplot(1, 3, 2); imshow(res_Gauss_IdealFilter); 
title(['2. 理想濾波器 (D0=' num2str(D0) ')']);
xlabel('評語：去噪強，但邊緣有不自然的振鈴(鬼影)');

subplot(1, 3, 3); imshow(res_Gauss_GaussFilter); 
title(['3. 高斯濾波器 (D0=' num2str(D0) ')']);
xlabel('評語：平滑模糊，但視覺自然，無波紋');

%% 5. 繪圖比較 - 第二戰：椒鹽雜訊
figure('Name', 'Battle 2: Salt & Pepper Noise', 'Position', [50, 550, 1200, 400]);

subplot(1, 3, 1); imshow(noisy_sp); title('1. 原圖 + 椒鹽雜訊');
subplot(1, 3, 2); imshow(res_SP_IdealFilter); 
title('2. 理想濾波器');
xlabel('評語：災難。黑白點變成了一圈圈漣漪');

subplot(1, 3, 3); imshow(res_SP_GaussFilter); 
title('3. 高斯濾波器');
xlabel('評語：也不好。黑白點暈開變成灰斑');

fprintf('所有圖片已儲存完畢！\n');

%% === 函式庫 (Functions) ===

% 1. 你的高斯雜訊函式
function output = add_g_noise(input, sigma)
    [m, n, z] = size(input);
    % 產生 RGB 三通道的隨機雜訊
    noise = randn(m, n, z); 
    output = input + (sigma * noise);
    output = max(0, min(output, 1)); % 截斷
end

% 2. 彩色濾波專用函式 (處理 RGB 三層)
function img_filtered = apply_filter_rgb(img_noisy, H_filter)
    [M, N, C] = size(img_noisy);
    img_filtered = zeros(M, N, C);
    
    for c = 1:C
        % 1. 取出單一顏色通道
        channel = img_noisy(:, :, c);
        
        % 2. 轉頻域
        F = fftshift(fft2(channel));
        
        % 3. 乘上濾波器
        F_filtered = F .* H_filter;
        
        % 4. 轉回空域
        res = real(ifft2(ifftshift(F_filtered)));
        
        % 5. 存回通道
        img_filtered(:, :, c) = res;
    end
    
    % 確保數值範圍正確
    img_filtered = max(0, min(img_filtered, 1));
end