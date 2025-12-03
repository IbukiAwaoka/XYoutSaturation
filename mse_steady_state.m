% S2/S1を変化させて定常状態のMSEを記録するプログラム
clc;
clear all;

tic;    % ストップウォッチタイマー開始

% 固定パラメータ
t = 5000;
N = 200;
ens = 100;
a = 10;  % メモリ幅
mu = 0.1;
xi = 0.0;
S1 = 2.0;

% S2/S1の範囲を設定
ratio_values = 0:0.1:3.0;  % 0から3.0まで0.1刻み
num_ratios = length(ratio_values);

% 結果を保存する配列
steady_state_mse = zeros(num_ratios, 1);

% 出力フォルダの作成
outputDir = 'output';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% 各S2/S1比率についてシミュレーション実行
for r_idx = 1:num_ratios
    ratio = ratio_values(r_idx);
    S2 = ratio * S1;
    
    fprintf('S2/S1 = %.2f (S2 = %.3f) を計算中...\n', ratio, S2);
    
    n = N * t + 1; % 更新回数
    MSE = zeros(a * t + 1, 1); % MSEの初期化
    
    for i = 1:ens
        % 初期化
        g = randn(N, 1);   % 未知システム
        w = zeros(N, 1);   % 適応フィルタ
        u = randn(N, 1) / sqrt(N); % N行1列の入力ベクトル
        
        for j = 1:n
            % 入力ベクトルの更新
            u = circshift(u, 1);
            u(1) = randn / sqrt(N);
            
            % 未知システムの出力
            x = g.' * u;
            
            % 非線形飽和型関数を適用
            if x > S1
                x_h = S1;
            elseif x < -S1
                x_h = -S1;
            else
                x_h = x;
            end
            
            % 適応フィルタの出力
            y = w.' * u;
            
            % 非線形飽和型関数を適用
            if y > S2
                y_h = S2;
            elseif y < -S2
                y_h = -S2;
            else
                y_h = y;
            end
            
            % 誤差の生成
            e = x_h - y_h + sqrt(xi) * randn;
            
            % MSEの計算
            if rem(j, N / a) == 1
                MSE((j - 1) / (N / a) + 1) = MSE((j - 1) / (N / a) + 1) + e^2 / ens;
            end
            
            % 適応フィルタの更新
            w = w + mu * e * u;
        end
    end
    
    % 定常状態のMSE（最後の値）を記録
    steady_state_mse(r_idx) = MSE(end);
    
    fprintf('  完了: 定常MSE = %.6f\n', steady_state_mse(r_idx));
end

% 結果の保存
data = [ratio_values; steady_state_mse'];
header = 'Steady State MSE';
fname = char([header, '(N=', num2str(N), ',T=', num2str(t), ',S1=', num2str(S1), ',mu=', num2str(mu), ',ens=', num2str(ens), ',xi=', num2str(xi), ').txt']);
fname = fullfile(outputDir, fname);

Fid = fopen(fname, 'w');
S = '#S2/S1     MSE';
fprintf(Fid, '%s\n', S);
fprintf(Fid, '%.2f   %1.6f\n', data);
fclose(Fid);

fprintf('\n結果を保存しました: %s\n', fname);

toc;    % タイマー終了
