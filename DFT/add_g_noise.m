function [output] = add_noise(input, sigma)
    % input: 輸入影像 (建議先轉為 double，範圍 0~1)
    % sigma: 標準差 (Standard Deviation)，控制雜訊強度
    %        建議值：0.01 (輕微) ~ 0.1 (強烈)
    %        注意：不要傳入原本的 r=20，那樣雜訊會大到看不見圖片！

    %% 1. 產生高斯雜訊 (Gaussian Noise Generation)
    % randn 會產生平均值為 0，標準差為 1 的常態分佈亂數
    [m, n, z] = size(input);
    
    % 如果想要每次雜訊都一樣以便除錯，可以把下面這行取消註解
     rng(1); 
    
    noise = randn(m, n, z); 
    
    %% 2. 將雜訊加入影像 (Additive Noise)
    % 公式：Image_new = Image_old + (強度 * Noise)
    output = input + (sigma * noise);
    
    %% 3. 數值截斷 (Clipping)
    % 加完雜訊後，數值可能會超過 1 或小於 0，必須修剪回來
    output = max(0, min(output, 1));
end