clc;
clear all;

% outputフォルダの作成
if ~exist('output', 'dir')
    mkdir('output');
end

% 並列プールの開始
if isempty(gcp('nocreate'))
    poolobj = parpool('local');
    fprintf('並列プール起動: %d コアを使用中\n', poolobj.NumWorkers);
else
    poolobj = gcp;
    fprintf('既存の並列プールを使用: %d コアを使用中\n', poolobj.NumWorkers);
end

t = 500;
N = 200;
ens = 100;

mu = 0.1; % ステップサイズ
xi = 0; % 背景雑音

% S1とS2の範囲を0～3まで0.1刻みで設定
S1_values = 0:0.1:3;
S2_values = 0:0.1:3;

num_S1_values = numel(S1_values);
num_S2_values = numel(S2_values);

% MSEの値を保存する行列 (行:S2, 列:S1)
MSE_matrix = zeros(num_S2_values, num_S1_values);

n = N*t+1; % 更新回数

tic; % ストップウォッチタイマー開始

total_iterations = num_S1_values * num_S2_values;

% S1とS2の二重ループを並列化
parfor idx_S1 = 1:num_S1_values
    S1 = S1_values(idx_S1);
    
    MSE_column = zeros(num_S2_values, 1); % 各S1に対するS2のMSE列
    
    for idx_S2 = 1:num_S2_values
        S2 = S2_values(idx_S2);
        
        MSE = zeros(ens, 1); % MSE初期化
        
        for i = 1:ens
            g = randn(N, 1); % 未知システム
            w = zeros(N, 1); % 適応フィルタの初期化
            x = randn(N, 1) / sqrt(N); % N行1列の入力ベクトル
            
            for j = 1:n
                x = circshift(x, 1); % 入力ベクトルをすべてずらす
                x(1) = randn / sqrt(N); % 値を上書きする
                
                d = g.' * x; % 未知システムの出力
                y = w.' * x; % 適応フィルタの出力
                
                % S1による飽和処理
                if d < -S1
                    d = -S1;
                elseif d > S1
                    d = S1;
                end
                
                % S2による飽和処理
                if y < -S2
                    y = -S2;
                elseif y > S2
                    y = S2;
                end
                
                e = d - y + sqrt(xi) * randn; % 誤差の生成
                
                if j >= n-9
                    MSE(i) = MSE(i) + e^2; % 最後の10回分の誤差の二乗を足す
                end
                
                w = w + mu * e * x; % 適応フィルタの更新
            end
            MSE(i) = MSE(i)/10; % 最後の10回分の平均を計算
        end
        
        MSE_column(idx_S2) = median(MSE); % MSEの中央値を保存
    end
    
    MSE_matrix(:, idx_S1) = MSE_column; % 列ごとに結果を保存
    
    fprintf('S1=%d/%d 完了\n', idx_S1, num_S1_values);
end

toc; % タイマー終了

% ヒートマップの描画
figure('Position', [100, 100, 800, 600]);
imagesc(S1_values, S2_values, MSE_matrix);
set(gca, 'YDir', 'normal'); % Y軸を正常な向きに
colorbar;
colormap('jet'); % カラーマップの設定（他に'hot', 'parula'なども選択可能）

xlabel('S1', 'FontSize', 14);
ylabel('S2', 'FontSize', 14);
title(['準定常MSEヒートマップ (N=', num2str(N), ', T=', num2str(t), ', \mu=', num2str(mu), ', ens=', num2str(ens), ')'], 'FontSize', 14);

% グリッドの追加
grid on;
set(gca, 'FontSize', 12);

% ヒートマップの保存
heatmap_fname = fullfile('output', ['MSE_heatmap(N=', num2str(N), ', T=', num2str(t), ', mu=', num2str(mu), ', ens=', num2str(ens), ', xi=', num2str(xi), ').png']);
saveas(gcf, heatmap_fname);

% データの保存
data_fname = fullfile('output', ['MSE_heatmap_data(N=', num2str(N), ', T=', num2str(t), ', mu=', num2str(mu), ', ens=', num2str(ens), ', xi=', num2str(xi), ').mat']);
save(data_fname, 'S1_values', 'S2_values', 'MSE_matrix');

disp(['ヒートマップを保存しました: ', heatmap_fname]);
disp(['データを保存しました: ', data_fname]);