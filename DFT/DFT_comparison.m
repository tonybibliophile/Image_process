clear all;
close all;
clc;

%% 1. 設定測試環境
N_size = 128; % 建議設 64 或 80，設 128 純迴圈會跑很久
fprintf('=== DFT 實作效能評比 (N = %d) ===\n', N_size);

% 產生測試圖片
img = rand(N_size, N_size, 3); 
img = im2double(img);

%% 2. 執行測試 (順序沒差，重點是最後怎麼排)

% --- 選手 A: 矩陣運算法 (Matrix) ---
fprintf('1. 正在執行 [矩陣運算法] (Matrix)... ');
tic;
F_matrix = FT_Matrix(img);
t_matrix = toc;
fprintf('完成! 耗時: %.4f 秒\n', t_matrix);

% --- 選手 B: 半向量化法 (User Code) ---
fprintf('2. 正在執行 [半向量化法] (User Loop+Sum)... ');
h = waitbar(0, '正在執行 User Code...');
tic;
F_user = FT_UserCode(img, h);
t_user = toc;
close(h);
fprintf('完成! 耗時: %.4f 秒\n', t_user);

% --- 選手 C: 純 For 迴圈法 (Pure Loops) ---
if N_size > 128
    fprintf('3. [純 For 迴圈法]: N=%d 太大了，略過以防當機。\n', N_size);
    t_loop = NaN; 
else
    fprintf('3. 正在執行 [純 For 迴圈法]... ');
    tic;
    F_loop = FT_PureLoop(img);
    t_loop = toc;
    fprintf('完成! 耗時: %.4f 秒\n', t_loop);
end

%% 3. 結果繪圖比較 (由快到慢排序)
figure('Name', 'Performance: Fastest to Slowest', 'Position', [200, 200, 700, 450]);

% --- [關鍵修改] 這裡決定顯示順序 ---
if isnan(t_loop)
    % 如果純迴圈沒跑，只比前兩名
    times = [t_matrix, t_user];
    labels = {'1. Matrix (最快)', '2. User Code (中等)'};
    bar_colors = [0.2 0.8 0.4;  0.2 0.6 0.8]; % 綠 -> 藍
else
    % 如果全跑了，依序排列: Matrix -> User -> Loop
    times = [t_matrix, t_user, t_loop];
    labels = {'1. Matrix (最快)', '2. User Code (中等)', '3. Pure Loop (最慢)'};
    bar_colors = [0.2 0.8 0.4;  0.2 0.6 0.8;  0.8 0.2 0.2]; % 綠 -> 藍 -> 紅
end

% 繪製長條圖
b = bar(times);
b.FaceColor = 'flat';

% 套用顏色
for i = 1:length(times)
    b.CData(i,:) = bar_colors(i,:);
end

set(gca, 'XTickLabel', labels, 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Time (seconds)');
title(sprintf('DFT Performance Ranking (N=%d) - Lower is Better', N_size));
grid on;

% 標示數值
for i = 1:length(times)
    text(i, times(i), num2str(times(i), '%.4f s'), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
end

% 顯示倍率比較
fprintf('\n=== 速度大比拚 ===\n');
fprintf('冠軍 (Matrix) 比 亞軍 (User Code) 快了 %.2f 倍\n', t_user / t_matrix);
if ~isnan(t_loop)
    fprintf('冠軍 (Matrix) 比 季軍 (Pure Loop) 快了 %.2f 倍\n', t_loop / t_matrix);
end


%% === 函式實作區 (保持不變) ===

% 1. 矩陣法 (O(N^3) - 最快)
function [I_freq] = FT_Matrix(I)
    [M, N, C] = size(I);
    I_freq = zeros(M, N, C);
    m_idx = 0:M-1; k_m = m_idx'; W_M = exp(-1j * 2 * pi * k_m * m_idx / M);
    n_idx = 0:N-1; k_n = n_idx'; W_N = exp(-1j * 2 * pi * k_n * n_idx / N);
    for c = 1:C
        I_freq(:,:,c) = W_M * double(I(:,:,c)) * W_N.';
    end
end

% 2. 半向量化法 (O(N^4) 但有 sum 加速 - 中間)
function [I_freq] = FT_UserCode(I, h_waitbar)
    [m,n,c] = size(I);
    temp = zeros(m,n,c);
    total_steps = c * m; step_count = 0;
    
    for z = 1:c
        for u = 1:m
            if mod(u, 10) == 0 && ~isempty(h_waitbar)
                step_count = step_count + 10;
                waitbar(step_count/total_steps, h_waitbar);
            end
            for v = 1:n
                l = 0:m-1; j = 0:n-1;
                [L,J] = meshgrid(l,j); 
                E = exp(-1i*2*pi*(L*(u-1)/m + J*(v-1)/n));
                temp(u,v,z) = sum(sum(I(:,:,z).*E));
            end
        end
    end
    I_freq = temp;
end

% 3. 純迴圈法 (O(N^4) 純純的慢 - 最慢)
function [I_freq] = FT_PureLoop(I)
    [M, N, C] = size(I);
    I_freq = zeros(M, N, C);
    for z = 1:C
        for u = 0:M-1
            for v = 0:N-1
                sum_val = 0;
                for x = 0:M-1
                    for y = 0:N-1
                        angle = -2*pi * (u*x/M + v*y/N);
                        sum_val = sum_val + I(x+1, y+1, z) * exp(1j * angle);
                    end
                end
                I_freq(u+1, v+1, z) = sum_val;
            end
        end
    end
end