clear all;
close all;
clc;

%% 1. 準備工作
disp('=== 傅立葉轉換大亂鬥：四種方法效能評比 ===');

% 讀取圖片
try
    img = imread('lena.jpg');
catch
    img = rand(256, 256, 3);
end
img = im2double(img);

% [重要] 設定測試尺寸
% 建議設定 64 或 128。
% 絕對不要超過 128，否則方法 1 (暴力迴圈) 會跑到當機！
N_size = 256; 
img = imresize(img, [N_size, N_size]);
[M, N, C] = size(img);
fprintf('測試圖片尺寸: %d x %d x %d\n', M, N, C);
fprintf('注意：尺寸若大於 128，暴力法(Method 1)將耗時極久！\n\n');

%% 2. 比賽開始

% --- 方法 1: 暴力迴圈 DFT (Brute Force Loops) ---
disp('1. 正在執行 [暴力迴圈法] (O(N^4))... 請耐心等待...');
tic;
F1 = FT_Loops(img);
t1 = toc;
fprintf('>> 耗時: %.4f 秒\n', t1);

% --- 方法 2: 矩陣運算 DFT (Matrix DFT) ---
disp('2. 正在執行 [DFT 矩陣法] (O(N^3))...');
tic;
F2 = FT_Matrix(img);
t2 = toc;
fprintf('>> 耗時: %.4f 秒\n', t2);

% --- 方法 3: 遞迴 FFT (Recursive FFT) ---
disp('3. 正在執行 [遞迴 FFT] (O(N^2 logN) + Overhead)...');
tic;
F3 = FT_Recursive(img);
t3 = toc;
fprintf('>> 耗時: %.4f 秒\n', t3);

% --- 方法 4: 迭代 FFT (Iterative FFT) ---
disp('4. 正在執行 [迭代 FFT] (O(N^2 logN) Optimized)...');
tic;
F4 = FT_Iterative(img);
t4 = toc;
fprintf('>> 耗時: %.4f 秒\n', t4);

%% 3. 結果驗證 (確保大家算的都一樣)
% 以矩陣法為基準
err1 = norm(F1(:) - F2(:)) / numel(img);
err3 = norm(F3(:) - F2(:)) / numel(img);
err4 = norm(F4(:) - F2(:)) / numel(img);

fprintf('\n=== 準確度檢查 (誤差應接近 0) ===\n');
fprintf('Method 1 誤差: %.10f\n', err1);
fprintf('Method 3 誤差: %.10f\n', err3);
fprintf('Method 4 誤差: %.10f\n', err4);

%% 4. 繪製圖表
figure('Name', 'Four Methods Comparison', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 500]);

% 時間長條圖
times = [t1, t2, t3, t4];
b = bar(times);
b.FaceColor = 'flat';
b.CData(1,:) = [0.8 0.2 0.2]; % 紅 (慢)
b.CData(2,:) = [0.2 0.6 0.8]; % 藍 (中)
b.CData(3,:) = [0.9 0.6 0.2]; % 橘 (遞迴)
b.CData(4,:) = [0.2 0.8 0.4]; % 綠 (快)

set(gca, 'XTickLabel', {'Loop (N^4)', 'Matrix (N^3)', 'Rec. FFT', 'Iter. FFT'});
ylabel('Time (seconds)');
title(sprintf('Performance Comparison (N=%d)', N_size));
grid on;

% 標示數值
for i = 1:4
    text(i, times(i), num2str(times(i), '%.4f s'), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

%% === 演算法實作區 ===

% 1. 暴力迴圈法 (最慢)
function F = FT_Loops(I)
    [M, N, C] = size(I);
    F = zeros(M, N, C);
    % 四層迴圈 + 通道迴圈 = 5層
    for c = 1:C
        for u = 0:M-1
            for v = 0:N-1
                sum_val = 0;
                for x = 0:M-1
                    for y = 0:N-1
                        % 公式直接實作
                        angle = -2*pi * (u*x/M + v*y/N);
                        sum_val = sum_val + I(x+1, y+1, c) * exp(1j * angle);
                    end
                end
                F(u+1, v+1, c) = sum_val;
            end
        end
    end
end

% 2. DFT 矩陣法 (MATLAB 優化)
function F = FT_Matrix(I)
    [M, N, C] = size(I);
    F = zeros(M, N, C);
    m = 0:M-1; W_M = exp(-1j * 2 * pi * (m' * m) / M);
    n = 0:N-1; W_N = exp(-1j * 2 * pi * (n' * n) / N);
    for c = 1:C
        F(:,:,c) = W_M * I(:,:,c) * W_N.';
    end
end

% 3. 遞迴 FFT (邏輯簡單，開銷大)
function F = FT_Recursive(I)
    [M, N, C] = size(I);
    F = zeros(M, N, C);
    for c = 1:C
        temp = zeros(M, N);
        for m = 1:M, temp(m, :) = fft_rec_1d(I(m, :, c).').'; end
        for n = 1:N, F(:, n, c) = fft_rec_1d(temp(:, n)); end
    end
end

function X = fft_rec_1d(x)
    N = length(x);
    if N == 1, X = x; return; end
    X_even = fft_rec_1d(x(1:2:end));
    X_odd  = fft_rec_1d(x(2:2:end));
    k = (0 : N/2 - 1)'; W = exp(-1j * 2 * pi * k / N);
    term = W .* X_odd;
    X = [X_even + term; X_even - term];
end

% 4. 迭代 FFT (移除遞迴，速度最穩)
function F = FT_Iterative(I)
    [M, N, C] = size(I);
    F = zeros(M, N, C);
    for c = 1:C
        temp = zeros(M, N);
        for m = 1:M, temp(m, :) = fft_iter_1d(I(m, :, c).').'; end
        for n = 1:N, F(:, n, c) = fft_iter_1d(temp(:, n)); end
    end
end

function X = fft_iter_1d(x)
    N = length(x);
    % 1. 位元反轉排序 (Bit-Reversal)
    % 這一步是將遞迴的輸入順序轉為迭代的順序
    p = bitrevorder(1:N);
    X = x(p);
    
    % 2. 迭代運算 (Butterfly)
    stage = 1;
    while stage < N
        step_size = 2 * stage;
        % 預計算旋轉因子
        k = (0 : stage-1)';
        W = exp(-1j * pi * k / stage);
        
        % 向量化處理每一層的蝴蝶
        for start_idx = 1 : step_size : N
            upper_idx = start_idx : (start_idx + stage - 1);
            lower_idx = upper_idx + stage;
            
            upper = X(upper_idx);
            lower = X(lower_idx);
            
            term = W .* lower;
            X(upper_idx) = upper + term;
            X(lower_idx) = upper - term;
        end
        stage = step_size;
    end
end