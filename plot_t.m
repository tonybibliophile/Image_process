clear all; close all; clc;

%% 1. 實驗設定
% 測試的尺寸清單
N_list = [32, 64, 128, 256, 512, 1024, 2048]; 
num_tests = length(N_list);

% 預先分配空間存時間
t_loop      = nan(1, num_tests);
t_matrix    = nan(1, num_tests);
t_recursive = nan(1, num_tests);
t_iterative = nan(1, num_tests);
t_builtin   = nan(1, num_tests);

fprintf('=== 開始效能評測 (N 從 %d 到 %d) ===\n', N_list(1), N_list(end));

% --- [新增] 初始化進度條 ---
h_bar = waitbar(0, '準備開始...', 'Name', 'DFT Algorithm Benchmark');

%% 2. 迴圈測試
for i = 1:num_tests
    N = N_list(i);
    
    % --- [新增] 更新進度條文字 ---
    msg = sprintf('正在測試 N = %d (%d/%d)... 請稍候', N, i, num_tests);
    waitbar((i-1)/num_tests, h_bar, msg);
    
    fprintf('正在測試 N = %4d ... ', N);
    
    % 產生隨機單通道圖片
    img = rand(N, N);
    
    % --- 1. 暴力迴圈法 (只測到 N=64，不然會跑不完) ---
    if N <= 64
        tic; FT_Loops(img); t_loop(i) = toc;
    end
    
    % --- 2. 矩陣法 (DFT Matrix) ---
    tic; FT_Matrix(img); t_matrix(i) = toc;
    
    % --- 3. 遞迴 FFT (Recursive) ---
    tic; FT_Recursive(img); t_recursive(i) = toc;
    
    % --- 4. 迭代 FFT (Iterative - Turbo版) ---
    tic; FT_Iterative(img); t_iterative(i) = toc;
    
    % --- 5. 內建 FFT (Reference) ---
    tic; fft2(img); t_builtin(i) = toc;
    
    fprintf('完成 (矩陣法耗時: %.4fs)\n', t_matrix(i));
end

% --- [新增] 關閉進度條 ---
waitbar(1, h_bar, '測試完成！');
pause(0.5);
close(h_bar);

%% 3. 繪製趨勢圖
figure('Name', 'DFT Algorithm Performance Analysis', 'Position', [100, 100, 900, 600]);

% 使用 semilogy (Y軸 Log Scale) 
semilogy(N_list, t_loop,      'o-', 'LineWidth', 2, 'DisplayName', 'Brute Force Loop O(N^4)'); hold on;
semilogy(N_list, t_matrix,    's-', 'LineWidth', 2, 'DisplayName', 'Matrix DFT O(N^3)');
semilogy(N_list, t_recursive, '^-', 'LineWidth', 2, 'DisplayName', 'Recursive FFT O(N^2 logN)');
semilogy(N_list, t_iterative, 'd-', 'LineWidth', 2, 'DisplayName', 'Iterative FFT O(N^2 logN)');
semilogy(N_list, t_builtin,   'x-', 'LineWidth', 3, 'Color', 'k', 'DisplayName', 'Built-in fft2 (Compiled C)');

% 圖表美化
grid on;
xlabel('Image Size (N)');
ylabel('Execution Time (seconds) [Log Scale]');
title('Execution Time vs. Image Size (N)');
legend('Location', 'northwest', 'FontSize', 10);
xticks(N_list);
xticklabels(string(N_list));
xlim([N_list(1), N_list(end)]);

% 在圖上標註關鍵資訊
if ~isnan(t_loop(2))
    text(N_list(2), t_loop(2), '  \leftarrow Too slow!', 'Color', 'r', 'FontSize', 10);
end
text(N_list(end), t_matrix(end), '  \leftarrow Exponential Growth', 'FontSize', 10, 'VerticalAlignment', 'bottom');

fprintf('\n=== 測試完成，請查看跳出的圖表 ===\n');


%% === 演算法實作區 (Core Functions) ===

% 1. 暴力迴圈 (O(N^4))
function F = FT_Loops(I)
    [M, N] = size(I); F = zeros(M, N);
    for u = 0:M-1
        for v = 0:N-1
            sum_val = 0;
            for x = 0:M-1
                for y = 0:N-1
                    angle = -2*pi * (u*x/M + v*y/N);
                    sum_val = sum_val + I(x+1, y+1) * exp(1j * angle);
                end
            end
            F(u+1, v+1) = sum_val;
        end
    end
end

% 2. 矩陣法 (O(N^3))
function F = FT_Matrix(I)
    [M, N] = size(I);
    m = 0:M-1; W = exp(-1j * 2 * pi * (m' * m) / M);
    F = W * I * W.';
end

% 3. 遞迴 FFT (O(N^2 logN))
function F = FT_Recursive(I)
    [M, N] = size(I); F = zeros(M, N);
    for m = 1:M, F(m, :) = fft_rec_1d(I(m, :).').'; end
    for n = 1:N, F(:, n) = fft_rec_1d(F(:, n)); end
end
function X = fft_rec_1d(x)
    N = length(x); if N == 1, X = x; return; end
    X_even = fft_rec_1d(x(1:2:end)); X_odd  = fft_rec_1d(x(2:2:end));
    k = (0 : N/2 - 1)'; W = exp(-1j * 2 * pi * k / N);
    term = W .* X_odd; X = [X_even + term; X_even - term];
end

% 4. 迭代 FFT (Turbo Optimized)
function F = FT_Iterative(I)
    [M, N] = size(I); F = zeros(M, N);
    for m = 1:M, F(m, :) = fft_iter_1d(I(m, :).').'; end
    for n = 1:N, F(:, n) = fft_iter_1d(F(:, n)); end
end
function X = fft_iter_1d(x)
    N = length(x);
    % 純數值位元反轉 (最快)
    n_idx = 0:N-1; k_bits = log2(N); p_idx = zeros(1, N);
    for b = 1:k_bits
        p_idx = bitset(p_idx, k_bits - b + 1, bitget(n_idx, b));
    end
    X = x(p_idx + 1);
    
    stage = 1;
    while stage < N
        step_size = 2 * stage;
        k = (0 : stage-1)'; W = exp(-1j * pi * k / stage);
        for start_idx = 1 : step_size : N
            upper_idx = start_idx : (start_idx + stage - 1);
            lower_idx = upper_idx + stage;
            upper = X(upper_idx); lower = X(lower_idx);
            term = W .* lower;
            X(upper_idx) = upper + term; X(lower_idx) = upper - term;
        end
        stage = step_size;
    end
end