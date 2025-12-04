% MSEtheoryheatmap_1int.m
% 理論値版のS1×S2ヒートマップ生成スクリプト
% 
% 機能:
%   - S1とS2を0～3まで0.1刻みでスイープ
%   - parfeval並列計算により理論MSE値を計算（1次元積分使用）
%   - ヒートマップとして可視化・保存
%   - 統計データをテキストファイルとして保存
%
% 出力ファイル:
%   - output/MSEtheory_heatmap_t{tEnd}_mu{mu}.png
%   - output/MSEtheory_heatmap_data_t{tEnd}_mu{mu}.mat
%   - output/MSEtheory_heatmap_stats_t{tEnd}_mu{mu}.txt

clc;
clear all;

tic;

% パラメータ設定
tEnd = 300;
rho = 1;
xi = 0;
mu = 1;
sgm_g = sqrt(1);
PRS = 1e-10;
RNG = 8;

% 出力フォルダの作成
outputDir = 'output';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% S1とS2の範囲設定（0～3まで0.1刻み）
S1_values = 0:0.1:3;
S2_values = 0:0.1:3;
num_S1 = length(S1_values);
num_S2 = length(S2_values);
total_points = num_S1 * num_S2;

fprintf('=== 理論MSEヒートマップ計算開始（並列計算版 - parfeval） ===\n');
fprintf('t = %d\n', tEnd);
fprintf('μ = %g\n', mu);
fprintf('S1 範囲: %.1f ~ %.1f (刻み: 0.1)\n', S1_values(1), S1_values(end));
fprintf('S2 範囲: %.1f ~ %.1f (刻み: 0.1)\n', S2_values(1), S2_values(end));
fprintf('計算点数: %d x %d = %d\n', num_S1, num_S2, total_points);
fprintf('==========================================\n\n');

% MSE結果を保存する行列（行:S2, 列:S1）
MSE_matrix = zeros(num_S2, num_S1);

% 並列プールの起動
if isempty(gcp('nocreate'))
    parpool;
end

% parfevalを使った並列実行
fprintf('並列計算開始\n');
f(total_points) = parallel.FevalFuture;

% S1とS2の全組み合わせについてparfeval起動
idx = 0;
for idx_S1 = 1:num_S1
    for idx_S2 = 1:num_S2
        idx = idx + 1;
        S1 = S1_values(idx_S1);
        S2 = S2_values(idx_S2);
        f(idx) = parfeval(@computeMSEforS1S2, 3, S1, S2, tEnd, rho, xi, mu, sgm_g, PRS, RNG, idx_S1, idx_S2);
    end
end

% 結果を順次取得しながら進捗表示
for idx = 1:total_points
    [completedIdx, MSE_val, idx_S1_result, idx_S2_result] = fetchNext(f);
    MSE_matrix(idx_S2_result, idx_S1_result) = MSE_val;
    
    fprintf('完了: %d/%d (%.2f%%) - S1=%.1f, S2=%.1f, MSE=%.6e\n', idx, total_points, ...
        (idx/total_points)*100, S1_values(idx_S1_result), S2_values(idx_S2_result), MSE_val);
end

fprintf('\n=== 計算完了 ===\n');

% ヒートマップの描画と保存
fprintf('ヒートマップ描画中...\n');
figure('Position', [100, 100, 800, 600]);
imagesc(S1_values, S2_values, MSE_matrix);
set(gca, 'YDir', 'normal'); % Y軸を正常な向きに
colorbar;
colormap('jet');

xlabel('S1', 'FontSize', 14);
ylabel('S2', 'FontSize', 14);
title(sprintf('理論MSEヒートマップ (t=%d, μ=%g)', tEnd, mu), 'FontSize', 14);

grid on;
set(gca, 'FontSize', 12);

% ヒートマップファイル名
plot_fname = fullfile(outputDir, sprintf('MSEtheory_heatmap_t%d_mu%g.png', tEnd, mu));
saveas(gcf, plot_fname);
fprintf('ヒートマップ保存: %s\n', plot_fname);

% MATデータファイルの保存
mat_fname = fullfile(outputDir, sprintf('MSEtheory_heatmap_data_t%d_mu%g.mat', tEnd, mu));
save(mat_fname, 'S1_values', 'S2_values', 'MSE_matrix');
fprintf('データファイル保存: %s\n', mat_fname);

% 統計データテキストファイルの保存
fprintf('統計データ保存中...\n');
stats_fname = fullfile(outputDir, sprintf('MSEtheory_heatmap_stats_t%d_mu%g.txt', tEnd, mu));
fid = fopen(stats_fname, 'w');
fprintf(fid, '#S1 S2 MSE_mean MSE_median MSE_variance\n');

% MSE_matrix全体の統計量を計算
MSE_median_all = median(MSE_matrix(:));
MSE_variance_all = var(MSE_matrix(:));

% 各(S1, S2)ペアについて書き込み
for idx_S1 = 1:num_S1
    for idx_S2 = 1:num_S2
        S1 = S1_values(idx_S1);
        S2 = S2_values(idx_S2);
        MSE_val = MSE_matrix(idx_S2, idx_S1);
        
        % 理論値では各点で1つの値のみなので、平均値=その値自体
        % 中央値と分散は全体のMSE_matrixから計算
        fprintf(fid, '%.1f\t%.1f\t%.10e\t%.10e\t%.10e\n', ...
            S1, S2, MSE_val, MSE_median_all, MSE_variance_all);
    end
end

fclose(fid);
fprintf('統計データ保存: %s\n', stats_fname);

fprintf('\n=== 全処理完了 ===\n');
toc;

%% 並列計算用の関数（各ワーカーで実行される）
function [MSE_final, idx_S1, idx_S2] = computeMSEforS1S2(S1, S2, tEnd, rho, xi, mu, sgm_g, PRS, RNG, idx_S1, idx_S2)
    % 定数で表されるサンプル平均（S1に依存）
    fd2 = S1^2 - S1 * sqrt(2 * rho^2 * sgm_g^2 / pi) * ...
        exp(-(S1^2)/(2 * (rho^2) * sgm_g^2)) + ...
        (rho^2 * sgm_g^2 - S1^2) * erf(S1/sqrt(2 * rho^2 * sgm_g^2));
    
    dfd = rho^2 * sgm_g^2 * erf(S1/sqrt(2 * rho^2 * sgm_g^2));
    
    % ODEソルバーのオプション
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'InitialStep', 1e-3, 'MaxStep', 1);
    
    % ODEを解く
    [t, y] = ode45(@(t, y) odefun_local(t, y, S1, S2, fd2, dfd, rho, xi, mu, sgm_g, PRS, RNG), ...
        [0 tEnd], [1e-9, 1e-9], options);
    
    % 最終値を取得
    Q_final = y(end, 1);
    r_final = y(end, 2);
    
    % MSEを計算
    dycov = rho^2 * [sgm_g^2 r_final; r_final Q_final];
    invcov = inv(dycov);
    
    MSE_final = calculateMSE_local(Q_final, r_final, dycov, invcov, S2, rho, xi, fd2, S1, PRS, RNG);
end

%% ローカル関数群（parfor内で使用可能にするため）
function dydt = odefun_local(t, y, S1, S2, fd2, dfd, rho, xi, mu, sgm_g, PRS, RNG)
    dydt = zeros(2, 1);
    
    Q = max(y(1), 1e-12);
    r = max(y(2), 1e-12);
    
    erf_arg1 = S2/sqrt(2 * rho^2 * Q);
    erf_arg2 = S1/sqrt(2 * rho^2 * sgm_g^2);
    
    if isreal(erf_arg1) && isfinite(erf_arg1)
        fy2 = S2^2 - S2 * sqrt(2 * rho^2 * Q/pi) * ...
            exp(-S2^2/(2 * rho^2 * Q)) + (rho^2 * Q - S2^2) * erf(erf_arg1);
    else
        fy2 = 0;
    end
    
    if isreal(erf_arg1) && isfinite(erf_arg1)
        dfy = rho^2 * r * erf(erf_arg1);
    else
        dfy = 0;
    end
    
    if isreal(erf_arg2) && isfinite(erf_arg2)
        yfd = rho^2 * r * erf(erf_arg2);
    else
        yfd = 0;
    end
    
    if isreal(erf_arg1) && isfinite(erf_arg1)
        yfy = rho^2 * Q * erf(erf_arg1);
    else
        yfy = 0;
    end
    
    dycov = rho^2 * [sgm_g^2 r; r Q];
    invcov = inv(dycov);
    
    fdfy = calculateIntegral_local(Q, dycov, invcov, S1, S2, PRS, RNG, rho);
    
    dydt(1) = mu^2 * rho^2 * (fd2 + fy2 - 2 * fdfy + xi) + 2 * mu * (yfd - yfy);
    dydt(2) = mu * (dfd - dfy);
end

function result = calculateMSE_local(Q, r, dycov, invcov, S2, rho, xi, fd2, S1, PRS, RNG)
    erf_arg1 = S2/sqrt(2 * rho^2 * Q);
    if isreal(erf_arg1) && isfinite(erf_arg1)
        fy2 = S2^2 - S2 * sqrt(2 * rho^2 * Q/pi) * ...
            exp(-S2^2/(2 * rho^2 * Q)) + (rho^2 * Q - S2^2) * erf(erf_arg1);
    else
        fy2 = 0;
    end
    
    fdfy = calculateIntegral_local(Q, dycov, invcov, S1, S2, PRS, RNG, rho);
    result = fd2 + fy2 - 2 * fdfy + xi;
end

function result = calculateIntegral_local(Q, dycov, invcov, S1, S2, PRS, RNG, rho)
    yr = RNG * sqrt(rho^2) * max(sqrt(Q), 1e-12);
    
    result1 = integral(@(y) integratedXLower_local(y, S1, dycov, invcov, S2), -yr, -S2, ...
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result2 = integral(@(y) integratedXLower_local(y, S1, dycov, invcov, S2), -S2, S2, ...
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result3 = integral(@(y) integratedXLower_local(y, S1, dycov, invcov, S2), S2, yr, ...
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
        
    result4 = integral(@(y) integratedXMiddle_local(y, S1, S1, dycov, invcov, S2), -yr, -S2, ...
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result5 = integral(@(y) integratedXMiddle_local(y, S1, S1, dycov, invcov, S2), -S2, S2, ...
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result6 = integral(@(y) integratedXMiddle_local(y, S1, S1, dycov, invcov, S2), S2, yr, ...
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
        
    result7 = integral(@(y) integratedXUpper_local(y, S1, dycov, invcov, S2), -yr, -S2, ...
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result8 = integral(@(y) integratedXUpper_local(y, S1, dycov, invcov, S2), -S2, S2, ...
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result9 = integral(@(y) integratedXUpper_local(y, S1, dycov, invcov, S2), S2, yr, ...
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
        
    result = result1 + result2 + result3 + result4 + result5 + ...
             result6 + result7 + result8 + result9;
end

function result = integratedXLower_local(y, x1, dycov, invcov, S2)
    a = invcov(1,1);
    b = invcov(1,2);
    c = invcov(2,2);
    detCov = det(dycov);
    
    if y <= -S2
        coef = -S2;
    elseif y >= S2
        coef = S2;
    else
        coef = y;
    end
    
    erf_arg = sqrt(a) * (-x1 + b*y/a) / sqrt(2);
    
    if isreal(erf_arg) && all(isfinite(erf_arg))
        result = (-x1) * coef / (2*pi*sqrt(detCov)) * ...
            sqrt(pi./a./2) .* exp(-(c - b.^2./a).*y.^2./2) .* ...
            (1 + erf(erf_arg));
    else
        result = 0;
    end
end

function result = integratedXMiddle_local(y, x1, x2, dycov, invcov, S2)
    a = invcov(1,1);
    b = invcov(1,2);
    c = invcov(2,2);
    detCov = det(dycov);
    
    if y <= -S2
        coef = -S2;
    elseif y >= S2
        coef = S2;
    else
        coef = y;
    end
    
    erf_arg1 = sqrt(a./2) .* (x1 + b.*y./a);
    erf_arg2 = sqrt(a./2) .* (-x1 + b.*y./a);
    
    if isreal(erf_arg1) && all(isfinite(erf_arg1)) && isreal(erf_arg2) && all(isfinite(erf_arg2))
        result = coef ./ (2.*pi.*sqrt(detCov)) .* exp(-(c - b.^2./a).*y.^2./2) .* ...
            ((-1./a) .* (exp(-a.*(x1 + b.*y./a).^2./2) - exp(-a.*(-x1 + b.*y./a).^2./2)) - ...
            (b.*y./a) .* sqrt(pi./2./a) .* (erf(erf_arg1) - erf(erf_arg2)));
    else
        result = 0;
    end
end

function result = integratedXUpper_local(y, x2, dycov, invcov, S2)
    S1 = x2;
    a = invcov(1,1);
    b = invcov(1,2);
    c = invcov(2,2);
    detCov = det(dycov);
    
    if y <= -S2
        coef = -S2;
    elseif y >= S2
        coef = S2;
    else
        coef = y;
    end
    
    erf_arg = sqrt(a) * (S1 + b*y/a) / sqrt(2);
    
    if isreal(erf_arg) && all(isfinite(erf_arg))
        result = S1 * coef / (2*pi*sqrt(detCov)) * ...
            sqrt(pi./a./2) .* exp(-(c - b.^2./a).*y.^2./2) .* ...
            (1 - erf(erf_arg));
    else
        result = 0;
    end
end
