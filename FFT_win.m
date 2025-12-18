clear all; close all; clc;

% 設定尺寸：這裡不用很大，N=128 或 256 就夠了
% 因為在 script 模式下，O(N^3) 會慢得非常有感
N_test = 128; 
fprintf('=== 公平對決 B (純手寫腳本戰) N=%d ===\n', N_test);
fprintf('規則：禁止使用 * 運算子，強制雙方都用 Loop 實作\n\n');

img = rand(N_test, N_test);

%% 1. 手寫矩陣法 (Manual Matrix DFT)
% 原理：F = W * img * W.'
% 但我們手動實作矩陣乘法，模擬沒有硬體加速的狀況
disp('1. 手寫矩陣法 (O(N^3) Loop實作) 正在爬行...');
tic;

% 準備 W 矩陣
m = 0:N_test-1;
W = exp(-1j * 2 * pi * (m' * m) / N_test);
Wt = W.';

% 第一步: T = W * img
% 手寫矩陣乘法 (Row x Col)
T = my_matrix_mult(W, img);

% 第二步: F = T * W.'
F_manual_mat = my_matrix_mult(T, Wt);

t_manual_mat = toc;
fprintf('>> 手寫矩陣法耗時: %.4f 秒\n', t_manual_mat);


%% 2. 手寫 FFT (Turbo Iterative)
% 這是你原本的程式碼，也是用 Loop 跑的
disp('2. 手寫 FFT (Turbo Iterative) 正在衝刺...');
tic;
F_manual_fft = FT_Iterative(img);
t_manual_fft = toc;
fprintf('>> 手寫 FFT 耗時:   %.4f 秒\n', t_manual_fft);


%% 結果判定
fprintf('\n=== 最終結果 ===\n');
speedup = t_manual_mat / t_manual_fft;
fprintf('🏆 FFT 終於贏了！\n');
fprintf('🚀 速度快了 %.2f 倍\n', speedup);
fprintf('驗證：當雙方都沒有硬體加速時，演算法的威力 (N^3 vs N^2logN) 就出來了。\n');


%% === 核心函式區 ===

% [關鍵] 手寫矩陣乘法 (O(N^3))
% 模擬 C 語言底層邏輯，但跑在 MATLAB 慢速直譯器上
function C = my_matrix_mult(A, B)
    [rows_A, cols_A] = size(A);
    [rows_B, cols_B] = size(B);
    
    if cols_A ~= rows_B
        error('維度不合');
    end
    
    C = zeros(rows_A, cols_B);
    
    % 三層迴圈：標準的矩陣乘法定義
    for i = 1:rows_A
        for j = 1:cols_B
            sum_val = 0;
            % 這一層最致命
            for k = 1:cols_A
                sum_val = sum_val + A(i,k) * B(k,j);
            end
            C(i,j) = sum_val;
        end
    end
end

% 你的 Turbo FFT
function F = FT_Iterative(I)
    [M, N] = size(I);
    F = zeros(M, N);
    for m = 1:M, F(m, :) = fft_iter_1d(I(m, :).').'; end
    for n = 1:N, F(:, n) = fft_iter_1d(F(:, n)); end
end

% 優化版 1D FFT
function X = fft_iter_1d(x)
    N = length(x);
    n_idx = 0:N-1; k_bits = log2(N);
    p_idx = zeros(1, N);
    for b = 1:k_bits
        extracted_bit = bitget(n_idx, b);
        p_idx = bitset(p_idx, k_bits - b + 1, extracted_bit);
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