clc; clear; close all;

% 並列プールの開始
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool;
end

% Ioutputフォルダの作成（親ディレクトリ）
if ~exist('Ioutput', 'dir')
    mkdir('Ioutput');
end

% パラメータ設定
t = 5000;    % シミュレーション時間
N = 200;     % 次元数
ens = 100;  % 試行回数
mu = 0.1;    % ステップサイズ
xi = 0.0;      % 背景雑音
n = N * t + 1; % 更新回数

% 飽和パラメータの範囲
S1_vals = linspace(0, 3, 31); % 初期値，終了値，分割数
S2_vals = linspace(0, 3, 31); % 初期値，終了値，分割数
MSE_values = zeros(length(S1_vals), length(S2_vals)); % MSE を格納する行列

total_start_time = tic; % 全体の計算時間測定開始

for idx1 = 1:length(S1_vals)
    S1 = S1_vals(idx1);
    for idx2 = 1:length(S2_vals)
        S2 = S2_vals(idx2);

        loop_start_time = tic; % ループごとの時間測定開始

        MSE_array = zeros(ens, 1); % 各試行のMSEを格納

        parfor i = 1:ens
            g = randn(N, 1); % 未知システム
            w = zeros(N, 1); % 適応フィルタの初期化
            x = randn(N, 1) / sqrt(N); % 入力ベクトル

            for j = 1:n
                x = circshift(x, 1); % 入力ベクトルをすべてずらす
                x(1) = randn / sqrt(N); % 値を上書きする
                d = g.' * x; % 未知システムの出力
                y = w.' * x; % 適応フィルタの出力

                % 飽和非線形性の適用
                d = max(min(d, S1), -S1);
                y = max(min(y, S2), -S2);

                e = d - y + sqrt(xi) * randn; % 誤差の計算
                w = w + mu * e * x; % 適応フィルタの更新

                % 最終時刻の MSE
                if j == n
                    MSE_array(i) = e^2;
                end
            end
        end

        % 平均 MSE を計算して保存
        MSE_values(idx1, idx2) = sum(MSE_array) / ens;

        % 経過時間の取得
        elapsed_time = toc(loop_start_time);
        fprintf('S1 = %.2f, S2 = %.2f, MSE = %.5f\n, Time = %.2f sec\n', S1, S2, MSE_values(idx1, idx2), elapsed_time);
    end
end

% データの保存
filename = sprintf('mse_data_xi=%.2f.txt', xi);
fname = fullfile('Ioutput', filename);
fid = fopen(fname, 'w');
fprintf(fid, '# S1 S2 MSE\n');

for idx1 = 1:length(S1_vals)
    for idx2 = 1:length(S2_vals)
        fprintf(fid, '%.2f %.2f %.5f\n', S1_vals(idx1), S2_vals(idx2), MSE_values(idx1, idx2));
    end
    fprintf(fid, "\n"); % S1 のブロックが終わるごとに空行を挿入
end

fclose(fid);

total_elapsed_time = toc(total_start_time);
fprintf('データ保存完了: %s\n', filename);
fprintf('全体の計算時間: %.2f sec\n', total_elapsed_time);

% 並列プールのクローズ
delete(gcp('nocreate'));
