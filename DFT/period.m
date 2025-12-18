clear all; close all; clc;

%% 1. 模擬摩爾紋/印刷網點 (Simulate Halftone/Moiré)
% 讀取圖片
try
    img = imread('Zebra.jpg'); % 辣椒圖色彩豐富，效果較好
catch
    img = imread('Zebra.jpg');
end
img = im2double(img);
if size(img, 3) == 1, img = cat(3, img, img, img); end

[M, N, C] = size(img);

% --- 產生高頻網點干擾 ---
% 摩爾紋通常是「乘法性」的干擾，模擬印刷的網格
% 設定網點頻率 (數字越小，網點越密，越像印刷品)
freq = 4; 
[X, Y] = meshgrid(1:N, 1:M);

% 建立網點紋路 (Checkerboard pattern)
% 這模擬了報紙印刷的半色調 (Halftone) 效果
pattern = 0.5 * (sin(2*pi*X/freq) .* sin(2*pi*Y/freq)) + 1;

% 將紋路「乘」上去，而不是「加」上去
% 這樣更像真實的摩爾紋/網點物理現象
img_moire = img .* repmat(pattern, [1, 1, C]);

% 稍微調整亮度讓圖片好看一點
img_moire = (img_moire - min(img_moire(:))) / (max(img_moire(:)) - min(img_moire(:)));

imwrite(img_moire, 'Moire_Input.png');

%% 2. 建立頻率域濾波器
% 建立頻率座標網格
u = 0:(M-1); v = 0:(N-1);
idx_u = find(u > M/2); u(idx_u) = u(idx_u) - M;
idx_v = find(v > N/2); v(idx_v) = v(idx_v) - N;
[V, U] = meshgrid(v, u);
D = sqrt(U.^2 + V.^2);
D = fftshift(D); 

% --- 設定截止頻率 D0 ---
% 必須設得比網點頻率低。
% 這裡 freq=4，代表週期是4個像素。頻譜位置大約在 N/4 處。
% 假設圖片 512x512，干擾點約在 128 附近。
% 我們把 D0 設在 80，試圖切斷它。
D0 = 80; 

% A. 理想濾波器 (Ideal LPF)
H_Ideal = double(D <= D0);

% B. 高斯濾波器 (Gaussian LPF)
H_Gauss = exp(-(D.^2) / (2 * D0^2));

%% 3. 執行濾波
res_Ideal = apply_filter_rgb(img_moire, H_Ideal);
res_Gauss = apply_filter_rgb(img_moire, H_Gauss);

%% 4. 顯示結果比較
figure('Name', 'Moire Pattern Removal', 'Position', [100, 100, 1400, 500]);

% 原圖 + 摩爾紋
subplot(1, 3, 1); 
imshow(img_moire); 
title('1. 帶有摩爾紋/網點的圖像');
xlabel('模擬掃描印刷品的效果');

% 理想濾波器
subplot(1, 3, 2); 
imshow(res_Ideal); 
title(['2. 理想濾波器 (D0=' num2str(D0) ')']);
xlabel('評語：網點去掉了，但物體邊緣有「光暈/振鈴」');

% 高斯濾波器
subplot(1, 3, 3); 
imshow(res_Gauss); 
title(['3. 高斯濾波器 (D0=' num2str(D0) ')']);
xlabel('評語：像柔焦鏡頭，去除網點且視覺柔和');

%% 5. 頻譜分析 (關鍵：看摩爾紋在頻譜長什麼樣)
gray_moire = img_moire(:,:,2); % 取綠色通道分析
F = fftshift(fft2(gray_moire));
S = log(1 + abs(F));

figure('Name', 'Spectrum Analysis', 'Position', [100, 600, 500, 400]);
imshow(S, []); colormap jet; colorbar;
title('摩爾紋的頻譜特徵');
xlabel('請注意中心以外的 4 個對稱亮點，那就是摩爾紋');

%% === 函式庫 ===
function img_filtered = apply_filter_rgb(img_noisy, H_filter)
    [M, N, C] = size(img_noisy);
    img_filtered = zeros(M, N, C);
    for c = 1:C
        F = fftshift(fft2(img_noisy(:, :, c)));
        F_filtered = F .* H_filter;
        res = real(ifft2(ifftshift(F_filtered)));
        img_filtered(:, :, c) = res;
    end
    % 拉伸對比度以便觀察 (Normalize)
    img_filtered = (img_filtered - min(img_filtered(:))) / (max(img_filtered(:)) - min(img_filtered(:)));
end