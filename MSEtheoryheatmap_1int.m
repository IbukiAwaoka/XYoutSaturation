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

% S1とS2の範囲設定
S1_values = 0:0.1:3;
S2_values = 0:0.1:3;

num_S1_values = length(S1_values);
num_S2_values = length(S2_values);

% MSEの値を保存する行列 (行:S2, 列:S1)
MSE_matrix = zeros(num_S2_values, num_S1_values);

fprintf('=== 理論MSEヒートマップ計算開始（並列計算版 - parfeval） ===\n');
fprintf('t = %d, mu = %.2f\n', tEnd, mu);
fprintf('S1 範囲: %.2f ~ %.2f (刻み幅: %.1f)\n', S1_values(1), S1_values(end), S1_values(2)-S1_values(1));
fprintf('S2 範囲: %.2f ~ %.2f (刻み幅: %.1f)\n', S2_values(1), S2_values(end), S2_values(2)-S2_values(1));
fprintf('計算点数: %d x %d = %d\n', num_S1_values, num_S2_values, num_S1_values * num_S2_values);
fprintf('==========================================\n\n');

% 並列プールの起動
if isempty(gcp('nocreate'))
    poolobj = parpool;
    fprintf('並列プール起動: %d コアを使用中\n', poolobj.NumWorkers);
else
    poolobj = gcp;
    fprintf('既存の並列プールを使用: %d コアを使用中\n', poolobj.NumWorkers);
end

% parfevalを使った並列実行
fprintf('\n並列計算開始\n');
total_points = num_S1_values * num_S2_values;
f(total_points) = parallel.FevalFuture;

idx = 1;
for idx_S1 = 1:num_S1_values
    S1 = S1_values(idx_S1);
    for idx_S2 = 1:num_S2_values
        S2 = S2_values(idx_S2);
        f(idx) = parfeval(@computeMSEforS1S2, 1, S1, S2, tEnd, rho, xi, mu, sgm_g, PRS, RNG);
        idx = idx + 1;
    end
end

% 結果を順次取得しながら進捗表示
for idx = 1:total_points
    [completedIdx, MSE_val] = fetchNext(f);
    
    % インデックスを(idx_S1, idx_S2)に変換
    idx_S1 = floor((completedIdx - 1) / num_S2_values) + 1;
    idx_S2 = mod(completedIdx - 1, num_S2_values) + 1;
    
    MSE_matrix(idx_S2, idx_S1) = MSE_val;
    
    if mod(idx, 50) == 0 || idx == total_points
        fprintf('完了: %d/%d (%.2f%%) - S1=%.2f, S2=%.2f, MSE=%.6e\n', ...
            idx, total_points, (idx/total_points)*100, S1_values(idx_S1), S2_values(idx_S2), MSE_val);
    end
end

fprintf('\n=== 計算完了 ===\n');
toc;

% 統計情報の計算
MSE_median_all = median(MSE_matrix(:));
MSE_variance_all = var(MSE_matrix(:));

fprintf('\n統計情報:\n');
fprintf('MSE 最小値: %.6e\n', min(MSE_matrix(:)));
fprintf('MSE 最大値: %.6e\n', max(MSE_matrix(:)));
fprintf('MSE 中央値: %.6e\n', MSE_median_all);
fprintf('MSE 分散: %.6e\n', MSE_variance_all);

% 複素数のチェック
if ~isreal(MSE_matrix)
    fprintf('警告: MSE_matrixに複素数が含まれています。\n');
    % 複素数が含まれる位置を特定
    [row, col] = find(imag(MSE_matrix) ~= 0);
    fprintf('複素数が検出された位置: %d箇所\n', length(row));
    for k = 1:min(5, length(row))
        fprintf('  S1=%.2f, S2=%.2f: MSE=%.6e+%.6ei\n', ...
            S1_values(col(k)), S2_values(row(k)), ...
            real(MSE_matrix(row(k), col(k))), imag(MSE_matrix(row(k), col(k))));
    end
    MSE_matrix = real(MSE_matrix);
end

% txtファイルの保存（S1, S2, MSE形式）
txt_fname = fullfile(outputDir, ['MSEtheory_heatmap_t', num2str(tEnd), '_mu', num2str(mu), '.txt']);
Fid_txt = fopen(txt_fname, 'w');
fprintf(Fid_txt, '#S1\tS2\tMSE\n');

for idx_S1 = 1:num_S1_values
    for idx_S2 = 1:num_S2_values
        fprintf(Fid_txt, '%.2f\t%.2f\t%.10e\n', ...
            S1_values(idx_S1), S2_values(idx_S2), MSE_matrix(idx_S2, idx_S1));
    end
end
fclose(Fid_txt);
fprintf('txtデータ保存: %s\n', txt_fname);

% 統計データをテキストファイルに保存（正規化なし）
stats_fname = fullfile(outputDir, ['MSEtheory_heatmap_stats_t', num2str(tEnd), '_mu', num2str(mu), '.txt']);
Fid = fopen(stats_fname, 'w');
fprintf(Fid, '#S1\tS2\tMSE\tMSE_median(all)\tMSE_variance(all)\n');

for idx_S1 = 1:num_S1_values
    for idx_S2 = 1:num_S2_values
        fprintf(Fid, '%.2f\t%.2f\t%.10e\t%.10e\t%.10e\n', ...
            S1_values(idx_S1), S2_values(idx_S2), MSE_matrix(idx_S2, idx_S1), ...
            MSE_median_all, MSE_variance_all);
    end
end
fclose(Fid);
fprintf('統計データ保存: %s\n', stats_fname);

% MATデータの保存
data_fname = fullfile(outputDir, ['MSEtheory_heatmap_data_t', num2str(tEnd), '_mu', num2str(mu), '.mat']);
save(data_fname, 'S1_values', 'S2_values', 'MSE_matrix', 'tEnd', 'mu', 'rho', 'xi', 'sgm_g');
fprintf('MATデータ保存: %s\n', data_fname);

% ヒートマップの描画
figure('Position', [100, 100, 800, 600]);
imagesc(S1_values, S2_values, MSE_matrix);
set(gca, 'YDir', 'normal'); % Y軸を正常な向きに
colorbar;
colormap('jet');

xlabel('S1', 'FontSize', 14);
ylabel('S2', 'FontSize', 14);
title(['理論MSEヒートマップ (t=', num2str(tEnd), ', \mu=', num2str(mu), ')'], 'FontSize', 14);

grid on;
set(gca, 'FontSize', 12);

% ヒートマップの保存
heatmap_fname = fullfile(outputDir, ['MSEtheory_heatmap_t', num2str(tEnd), '_mu', num2str(mu), '.png']);
saveas(gcf, heatmap_fname);
fprintf('ヒートマップ保存: %s\n', heatmap_fname);

fprintf('\n=== 全処理完了 ===\n');

%% 並列計算用の関数（各ワーカーで実行される）
function MSE_final = computeMSEforS1S2(S1, S2, tEnd, rho, xi, mu, sgm_g, PRS, RNG)
    % S1またはS2が0の場合の特別処理
    if S1 < 1e-10 || S2 < 1e-10
        MSE_final = 1.0; % デフォルト値
        return;
    end
    
    % 定数で表されるサンプル平均（S1に依存）
    fd2 = S1^2 - S1 * sqrt(2 * rho^2 * sgm_g^2 / pi) * ...
        exp(-(S1^2)/(2 * (rho^2) * sgm_g^2)) + ...
        (rho^2 * sgm_g^2 - S1^2) * erf(S1/sqrt(2 * rho^2 * sgm_g^2));
    
    dfd = rho^2 * sgm_g^2 * erf(S1/sqrt(2 * rho^2 * sgm_g^2));
    
    % ODEソルバーのオプション（精度を上げる）
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'InitialStep', 1e-4, 'MaxStep', 0.5);
    
    try
        % ODEを解く
        [t, y] = ode45(@(t, y) odefun_local(t, y, S1, S2, fd2, dfd, rho, xi, mu, sgm_g, PRS, RNG), ...
            [0 tEnd], [1e-12, 1e-12], options);
        
        % 最終値を取得
        Q_final = max(real(y(end, 1)), 1e-12);
        r_final = real(y(end, 2));
        
        % 共分散行列が正定値かチェック
        dycov = rho^2 * [sgm_g^2 r_final; r_final Q_final];
        if det(dycov) <= 0
            warning('共分散行列が正定値ではありません: S1=%.2f, S2=%.2f', S1, S2);
            MSE_final = NaN;
            return;
        end
        
        invcov = inv(dycov);
        
        MSE_final = real(calculateMSE_local(Q_final, r_final, dycov, invcov, S2, rho, xi, fd2, S1, PRS, RNG));
    catch ME
        warning('計算エラー: S1=%.2f, S2=%.2f - %s', S1, S2, ME.message);
        MSE_final = NaN;
    end
end

%% ローカル関数群（parfor内で使用可能にするため）
function dydt = odefun_local(t, y, S1, S2, fd2, dfd, rho, xi, mu, sgm_g, PRS, RNG)
    dydt = zeros(2, 1);
    
    Q = max(real(y(1)), 1e-12);
    r = real(y(2));
    
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
    
    % 共分散行列が正定値かチェック
    if det(dycov) <= 0
        dydt(1) = 0;
        dydt(2) = 0;
        return;
    end
    
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