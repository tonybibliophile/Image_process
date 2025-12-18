clear all; close all; clc;

%% 0. 準備工作
save_dir = 'Final_Comparison_Output';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
fprintf('圖片將儲存至: %s\n', save_dir);

% 讀取圖片
try
    img = imread('Zebra.jpg');
catch
    img = imread('cameraman.tif');
end
img = im2double(img);
if size(img, 3) == 3, img = rgb2gray(img); end % 轉灰階方便觀察
[M, N] = size(img);

% 寫入原圖
imwrite(img, fullfile(save_dir, '0_Original.png'));

%% 1. 產生三種雜訊 (The 3 Noises)

% --- A. 高斯雜訊 (調用 add_g_noise) ---
sigma = 0.05; % 設定雜訊強度
% [重要] 這裡呼叫下面的 add_g_noise 函式
img_gauss = add_g_noise(img, sigma);

% --- B. 椒鹽雜訊 (Salt & Pepper) ---
density = 0.05;
img_sp = imnoise(img, 'salt & pepper', density);

% --- C. 週期性雜訊 (Periodic) ---
freq_u = 60; freq_v = 60; % 設定頻率
[X, Y] = meshgrid(1:N, 1:M);

% [修正處] 這裡變數名稱統一為 periodic_interference
periodic_interference = 0.3 * cos(2*pi*(freq_u*X/N + freq_v*Y/M));
img_periodic = img + periodic_interference;
img_periodic = max(0, min(img_periodic, 1));

%% 2. 準備兩種濾波器 (The 2 Filters)

% 頻率座標與距離 D
u = 0:(M-1); v = 0:(N-1);
idx_u = find(u > M/2); u(idx_u) = u(idx_u) - M;
idx_v = find(v > N/2); v(idx_v) = v(idx_v) - N;
[V, U] = meshgrid(v, u);
D = sqrt(U.^2 + V.^2);
D = fftshift(D); % 移到中心

% 設定截止頻率 (你可以改這裡來觀察不同程度的糊/振鈴)
D0 = 40; 

% --- Filter 1: 高斯濾波器 (Gaussian LPF) ---
H_Gauss = exp(-(D.^2) / (2 * D0^2));

% --- Filter 2: 理想濾波器 (Ideal LPF - 講義款) ---
H_Ideal = double(D <= D0);

% 顯示濾波器形狀
figure('Name', 'Filters', 'Position', [100, 100, 600, 300]);
subplot(1,2,1); imshow(H_Gauss, []); title('Gaussian Filter (軟)');
subplot(1,2,2); imshow(H_Ideal, []); title('Ideal Filter (硬)');

%% 3. 執行濾波 (處理 3 種雜訊 x 2 種濾波器)

% 準備雜訊圖陣列，方便跑迴圈
noise_images = {img_gauss, img_sp, img_periodic};
noise_names = {'Gaussian', 'SaltPepper', 'Periodic'};
filenames = {'A', 'B', 'C'};

% 儲存結果的容器
results_gauss_filter = cell(1, 3);
results_ideal_filter = cell(1, 3);

for i = 1:3
    current_noisy = noise_images{i};
    
    % FFT
    F = fftshift(fft2(current_noisy));
    
    % 1. 套用高斯濾波
    F_G = F .* H_Gauss;
    res_G = real(ifft2(ifftshift(F_G)));
    results_gauss_filter{i} = res_G;
    
    % 2. 套用理想濾波
    F_I = F .* H_Ideal;
    res_I = real(ifft2(ifftshift(F_I)));
    results_ideal_filter{i} = res_I;
    
    % 存檔
    imwrite(current_noisy, fullfile(save_dir, [filenames{i} '_Noise.png']));
    imwrite(res_G, fullfile(save_dir, [filenames{i} '_Result_GaussianFilter.png']));
    imwrite(res_I, fullfile(save_dir, [filenames{i} '_Result_IdealFilter.png']));
end

%% 4. 繪圖大比較 (3x3 Grid)
figure('Name', 'Filter Comparison Matrix', 'Position', [50, 50, 1000, 900]);

% --- Row 1: 高斯雜訊 ---
subplot(3, 3, 1); imshow(noise_images{1}); title('A. 高斯雜訊 (Gaussian)');
subplot(3, 3, 2); imshow(results_gauss_filter{1}); title('高斯濾波 (平滑自然)');
subplot(3, 3, 3); imshow(results_ideal_filter{1}); title('理想濾波 (邊緣有振鈴)');

% --- Row 2: 椒鹽雜訊 ---
subplot(3, 3, 4); imshow(noise_images{2}); title('B. 椒鹽雜訊 (Salt & Pepper)');
subplot(3, 3, 5); imshow(results_gauss_filter{2}); title('高斯濾波 (糊掉的髒點)');
subplot(3, 3, 6); imshow(results_ideal_filter{2}); title('理想濾波 (嚴重的漣漪)');

% --- Row 3: 週期性雜訊 ---
subplot(3, 3, 7); imshow(noise_images{3}); title('C. 週期性雜訊 (Periodic)');
subplot(3, 3, 8); imshow(results_gauss_filter{3}); title('高斯濾波 (有效但模糊)');
subplot(3, 3, 9); imshow(results_ideal_filter{3}); title('理想濾波 (振鈴明顯)');

fprintf('處理完成！圖片已儲存至 %s\n', save_dir);

%% === 函式區 ===
% 這是你要調用的 add_g_noise
function output = add_g_noise(input, sigma)
    % 模擬真實的高斯加性雜訊
    [m, n, z] = size(input);
    noise = randn(m, n, z); 
    output = input + (sigma * noise);
    output = max(0, min(output, 1)); % 截斷數值
end