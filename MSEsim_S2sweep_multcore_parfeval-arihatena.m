clc;
clear all;

tic;

% パラメータ設定
t = 5000;      % 十分に長い時間
N = 200;       % 次元
ens = 1000;     % 試行回数
a = 10;        % メモリ幅
mu = 0.1;      % ステップサイズ
xi = 0.0;      % 背景雑音
S1 = 2;        % 入力側飽和幅（固定）

% 出力フォルダの作成
outputDir = 'Ioutput';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% S2/S1の範囲設定
S2_S1_ratio = 0:0.1:3.0;  % 0から3.0まで0.1刻み
num_ratios = length(S2_S1_ratio);

% 結果保存用の配列
MSE_mean = zeros(num_ratios, 1);
MSE_std = zeros(num_ratios, 1);

fprintf('=== S2スイープ計算開始（並列計算版 - Simulation with Errorbar） ===\n');
fprintf('S1 = %.2f (固定)\n', S1);
fprintf('t = %d, N = %d, ens = %d\n', t, N, ens);
fprintf('S2/S1 範囲: %.2f ~ %.2f\n', S2_S1_ratio(1), S2_S1_ratio(end));
fprintf('計算点数: %d\n', num_ratios);
fprintf('==========================================\n\n');

% 並列プールの起動
if isempty(gcp('nocreate'))
    parpool;
end

% parfevalを使った並列実行
fprintf('並列計算開始\n');
f(num_ratios) = parallel.FevalFuture;
for idx = 1:num_ratios
    ratio = S2_S1_ratio(idx);
    S2 = ratio * S1;
    f(idx) = parfeval(@computeSteadyStateMSE, 2, S2, S1, t, N, ens, a, mu, xi);
end

% 結果を順次取得しながら進捗表示
for idx = 1:num_ratios
    [completedIdx, mse_mean_val, mse_std_val] = fetchNext(f);
    MSE_mean(completedIdx) = mse_mean_val;
    MSE_std(completedIdx) = mse_std_val;
    
    fprintf('完了: %d/%d (%.2f%%) - S2/S1=%.2f, MSE=%.6e±%.6e\n', idx, num_ratios, ...
        (idx/num_ratios)*100, S2_S1_ratio(completedIdx), mse_mean_val, mse_std_val);
end

% 結果をファイルに保存
fname = char(['parfeval_SteadyStateMSE_errorbar(N=', num2str(N), ', S1=', num2str(S1), ...
    ', T=', num2str(t), ', mu=', num2str(mu), ', ens=', num2str(ens), ', xi=', num2str(xi), ').txt']);
fname = fullfile(outputDir, fname);

Fid = fopen(fname, 'w');
fprintf(Fid, '#S2/S1 MSE_mean MSE_std\n');
for idx = 1:num_ratios
    fprintf(Fid, '%.2f\t%.6e\t%.6e\n', S2_S1_ratio(idx), MSE_mean(idx), MSE_std(idx));
end
fclose(Fid);

fprintf('\n=== 計算完了 ===\n');
fprintf('結果ファイル: %s\n', fname);
toc;

% プロット
figure('Position', [100, 100, 800, 600]);
errorbar(S2_S1_ratio, MSE_mean, MSE_std, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
grid on;
xlabel('S2/S1', 'FontSize', 14);
ylabel('Steady State MSE', 'FontSize', 14);
title(sprintf('MSE vs S2/S1 at t=%d (S1=%.2f fixed)', t, S1), 'FontSize', 16);
set(gca, 'YScale', 'log');

plot_fname = fullfile(outputDir, ['SteadyStateMSE_errorbar_S1', num2str(S1), '.png']);
saveas(gcf, plot_fname);
disp(['プロット保存: ', plot_fname]);

%% 並列計算用の関数
function [mse_mean, mse_std] = computeSteadyStateMSE(S2, S1, t, N, ens, a, mu, xi)
    % 各試行でのMSEを記録
    MSE_trials = zeros(ens, 1);
    
    n = N * t + 1; % 更新回数
    
    for trial = 1:ens
        % 初期化
        g = randn(N, 1);   % 未知システム
        w = zeros(N, 1);   % 適応フィルタ
        u = randn(N, 1) / sqrt(N);
        
        MSE_temp = zeros(a * t + 1, 1);
        
        for j = 1:n
            % 入力ベクトルの更新
            u = circshift(u, 1);
            u(1) = randn / sqrt(N);
            
            % 未知システムの出力
            x = g.' * u;
            if x > S1
                x_h = S1;
            elseif x < -S1
                x_h = -S1;
            else
                x_h = x;
            end
            
            % 適応フィルタの出力
            y = w.' * u;
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
                MSE_temp((j - 1) / (N / a) + 1) = e^2;
            end
            
            % 適応フィルタの更新
            w = w + mu * e * u;
        end
        
        % 定常状態のMSE（最後の値）を記録
        MSE_trials(trial) = MSE_temp(end);
    end
    
    % 平均と標準偏差を計算
    mse_mean = mean(MSE_trials);
    mse_std = std(MSE_trials);
end
