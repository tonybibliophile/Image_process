clear all;
close all;
clc;

%% 1. 準備工作
disp('=== 傅立葉轉換效能測試 (附進度條與 ETA) ===');

% 讀取圖片
try
    img = imread('lena.jpg');
catch
    img = rand(256, 256, 3);
end
img = im2double(img);

% [設定] 測試尺寸
% 建議先設 64 跑一次看看感覺。
% 若想挑戰 M1 極限，可設 128 (Method 1 會跑比較久，這時進度條就很有用)
% 設 256 你會看到進度條告訴你要跑好幾十分鐘...
N_size = 256; 
img = imresize(img, [N_size, N_size]);
[M, N, C] = size(img);
fprintf('測試圖片尺寸: %d x %d x %d\n\n', M, N, C);

%% 2. 比賽開始

% --- 方法 1: 暴力迴圈 DFT (附 ETA 進度條) ---
fprintf('1. 正在執行 [暴力迴圈法]... (請看跳出的進度視窗)\n');

% 呼叫下方附帶進度條的函式
tic;
F1 = FT_Loops_With_Bar(img); 
t1 = toc;

fprintf('>> Method 1 實際耗時: %.4f 秒\n', t1);


% --- 方法 2: 矩陣運算 DFT ---
disp('2. 正在執行 [DFT 矩陣法]...');
tic;
F2 = FT_Matrix(img);
t2 = toc;
fprintf('>> 耗時: %.4f 秒\n', t2);


% --- 方法 3: 遞迴 FFT ---
disp('3. 正在執行 [遞迴 FFT]...');
tic;
F3 = FT_Recursive(img);
t3 = toc;
fprintf('>> 耗時: %.4f 秒\n', t3);


% --- 方法 4: 迭代 FFT ---
disp('4. 正在執行 [迭代 FFT]...');
tic;
F4 = FT_Iterative(img);
t4 = toc;
fprintf('>> 耗時: %.4f 秒\n', t4);


%% 3. 繪圖結果
figure('Name', 'Execution Time', 'NumberTitle', 'off', 'Position', [200, 200, 600, 400]);
times = [t1, t2, t3, t4];
bar(times);
set(gca, 'XTickLabel', {'Loop (N^4)', 'Matrix (N^3)', 'Rec. FFT', 'Iter. FFT'});
ylabel('Time (s)');
title(['Performance (N=' num2str(N_size) ')']);
grid on;
for i=1:4, text(i, times(i), num2str(times(i),'%.3fs'),'Vert','bottom','Horiz','center'); end


%% === 核心函式區 ===

% [重點] 附帶進度條的暴力迴圈法
function F = FT_Loops_With_Bar(I)
    [M, N, C] = size(I);
    F = zeros(M, N, C);
    
    % 初始化進度條
    h = waitbar(0, '準備開始...', 'Name', 'Method 1 Progress');
    total_steps = C * M; % 總共要跑的「大迴圈」次數
    start_time = tic;    % 紀錄開始時間
    
    step_count = 0;
    
    for c = 1:C
        for u = 0:M-1
            % --- 內層運算 (最花時間的地方) ---
            for v = 0:N-1
                sum_val = 0;
                for x = 0:M-1
                    for y = 0:N-1
                        angle = -2*pi * (u*x/M + v*y/N);
                        sum_val = sum_val + I(x+1, y+1, c) * exp(1j * angle);
                    end
                end
                F(u+1, v+1, c) = sum_val;
            end
            % ------------------------------
            
            % 更新進度條 (每做完一個 Row 更新一次，避免更新太頻繁拖慢速度)
            step_count = step_count + 1;
            progress = step_count / total_steps;
            
            % 計算 ETA (預估剩餘時間)
            elapsed_time = toc(start_time);
            estimated_total_time = elapsed_time / progress;
            remaining_time = estimated_total_time - elapsed_time;
            
            % 更新顯示文字
            msg = sprintf('進度: %.1f%% | 已用: %.0fs | 預估剩餘 (ETA): %.0f 秒', ...
                          progress*100, elapsed_time, remaining_time);
            waitbar(progress, h, msg);
        end
    end
    
    % 關掉進度條
    close(h);
end

% 其他快速方法 (無需進度條)
function F = FT_Matrix(I)
    [M, N, C] = size(I); F = zeros(M, N, C);
    m = 0:M-1; W_M = exp(-1j * 2 * pi * (m' * m) / M);
    n = 0:N-1; W_N = exp(-1j * 2 * pi * (n' * n) / N);
    for c = 1:C, F(:,:,c) = W_M * I(:,:,c) * W_N.'; end
end

function F = FT_Recursive(I)
    [M, N, C] = size(I); F = zeros(M, N, C);
    for c = 1:C
        temp = zeros(M, N);
        for m = 1:M, temp(m, :) = fft_rec_1d(I(m, :, c).').'; end
        for n = 1:N, F(:, n, c) = fft_rec_1d(temp(:, n)); end
    end
end
function X = fft_rec_1d(x)
    N = length(x); if N == 1, X = x; return; end
    X_even = fft_rec_1d(x(1:2:end)); X_odd  = fft_rec_1d(x(2:2:end));
    k = (0 : N/2 - 1)'; W = exp(-1j * 2 * pi * k / N);
    term = W .* X_odd; X = [X_even + term; X_even - term];
end

function F = FT_Iterative(I)
    [M, N, C] = size(I); F = zeros(M, N, C);
    for c = 1:C
        temp = zeros(M, N);
        for m = 1:M, temp(m, :) = fft_iter_1d(I(m, :, c).').'; end
        for n = 1:N, F(:, n, c) = fft_iter_1d(temp(:, n)); end
    end
end
function X = fft_iter_1d(x)
    N = length(x); p = bitrevorder(1:N); X = x(p);
    stage = 1;
    while stage < N
        step_size = 2 * stage; k = (0 : stage-1)'; W = exp(-1j * pi * k / stage);
        for start_idx = 1 : step_size : N
            upper_idx = start_idx : (start_idx + stage - 1);
            lower_idx = upper_idx + stage;
            upper = X(upper_idx); lower = X(lower_idx);
            term = W .* lower; X(upper_idx) = upper + term; X(lower_idx) = upper - term;
        end
        stage = step_size;
    end
end