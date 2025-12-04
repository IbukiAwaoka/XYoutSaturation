clc;
clear all;

tEnd = 300;  % 時刻を300に設定

tic; % ストップウォッチタイマー開始

global fd2 dfd PRS RNG mu S1 S2 xi rho sgm_g

rho = 1;
S1 = 2;  % S1を固定
xi = 0;
mu = 1;
sgm_g = sqrt(1);

PRS = 1e-10; % 数値積分の精度 RelTol, AbsTol
RNG = 8; % 数値積分の範囲を決めるパラメータ

% 出力フォルダの作成
outputDir = 'output';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% S2/S1の範囲設定
S2_S1_ratio = linspace(0, 3, 31); % 0から3まで31点（0.1刻み）
num_ratios = length(S2_S1_ratio);

% 結果保存用の配列
MSE_results = zeros(num_ratios, 1);

% 出力ファイルの作成
fname = char(['MSE_S2sweep_t', num2str(tEnd), '_mu', num2str(mu), ...
    '_S1', num2str(S1), '.txt']);
fname = fullfile(outputDir, fname);
Fid = fopen(fname, 'w');
Header = '#S2/S1 S2 MSE(t=300)';
fprintf(Fid, '%s\n', Header);

fprintf('=== S2スイープ計算開始 ===\n');
fprintf('S1 = %.2f (固定)\n', S1);
fprintf('t = %d\n', tEnd);
fprintf('S2/S1 範囲: %.2f ~ %.2f\n', S2_S1_ratio(1), S2_S1_ratio(end));
fprintf('==========================\n\n');

% 各S2値について計算
for idx = 1:num_ratios
    S2 = S2_S1_ratio(idx) * S1;  % S2を計算
    
    fprintf('Progress: %d/%d (%.1f%%) | S2/S1 = %.3f, S2 = %.3f\n', ...
        idx, num_ratios, (idx/num_ratios)*100, S2_S1_ratio(idx), S2);
    
    % 定数で表されるサンプル平均（S1に依存）
    fd2 = S1^2 - S1 * sqrt(2 * rho^2 * sgm_g^2 / pi) * ...
        exp(-(S1^2)/(2 * (rho^2) * sgm_g^2)) + ... 
        (rho^2 * sgm_g^2 - S1^2) * erf(S1/sqrt(2 * rho^2 * sgm_g^2));
    
    dfd = rho^2 * sgm_g^2 * erf(S1/sqrt(2 * rho^2 * sgm_g^2));
    
    % ODEソルバーのオプション
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'InitialStep', 1e-3, 'MaxStep', 1);
    
    % ODEを解く
    [t, y] = ode45(@(t, y) odefun(t, y), [0 tEnd], [1e-9, 1e-9], options);
    
    % t=tEndでのQ, rを取得（最終値）
    Q_final = y(end, 1);
    r_final = y(end, 2);
    
    % MSEを計算
    dycov = rho^2 * [sgm_g^2 r_final; r_final Q_final];
    invcov = inv(dycov);
    
    MSE_final = calculateMSE(Q_final, r_final, dycov, invcov);
    MSE_results(idx) = MSE_final;
    
    % 結果をファイルに書き込み
    fprintf(Fid, '%.6f\t%.6f\t%.10e\n', S2_S1_ratio(idx), S2, MSE_final);
    
    fprintf('  -> Q(t=300) = %.6e, r(t=300) = %.6e, MSE = %.6e\n\n', ...
        Q_final, r_final, MSE_final);
end

fclose(Fid);

fprintf('\n=== 計算完了 ===\n');
disp(['結果ファイル: ', fname]);
toc;

% 結果のプロット
figure('Position', [100, 100, 800, 600]);
plot(S2_S1_ratio, MSE_results, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
grid on;
xlabel('S2/S1', 'FontSize', 14);
ylabel('MSE (t=300)', 'FontSize', 14);
title(sprintf('MSE vs S2/S1 at t=%d (S1=%.2f fixed)', tEnd, S1), 'FontSize', 16);

% プロットを保存
plot_fname = fullfile(outputDir, ['MSE_S2sweep_t', num2str(tEnd), '_S1', num2str(S1), '.png']);
saveas(gcf, plot_fname);
disp(['プロット保存: ', plot_fname]);

%% 補助関数群

function dydt = odefun(t, y)
    global fd2 dfd mu S1 S2 xi rho sgm_g
    
    dydt = zeros(2, 1);
    
    Q = max(y(1), 1e-12); % Qが負値にならないように制限
    r = max(y(2), 1e-12); % rが負値にならないように制限
    
    erf_arg1 = S2/sqrt(2 * rho^2 * Q);
    erf_arg2 = S1/sqrt(2 * rho^2 * sgm_g^2);
    
    % erf_argが実数で有限な値かをチェック
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
    
    fdfy = calculateIntegral(Q, dycov, invcov);
    
    dydt(1) = mu^2 * rho^2 * (fd2 + fy2 - 2 * fdfy + xi) + 2 * mu * (yfd - yfy);
    dydt(2) = mu * (dfd - dfy);
end

function result = calculateMSE(Q, r, dycov, invcov)
    global S2 rho xi fd2
    
    erf_arg1 = S2/sqrt(2 * rho^2 * Q);
    if isreal(erf_arg1) && isfinite(erf_arg1)
        fy2 = S2^2 - S2 * sqrt(2 * rho^2 * Q/pi) * ...
            exp(-S2^2/(2 * rho^2 * Q)) + (rho^2 * Q - S2^2) * erf(erf_arg1);
    else
        fy2 = 0;
    end
    
    fdfy = calculateIntegral(Q, dycov, invcov);
    result = fd2 + fy2 - 2 * fdfy + xi;
end

function result = calculateIntegral(Q, dycov, invcov)
    global S1 S2 PRS RNG rho
    
    % yに関する積分範囲
    yr = RNG * sqrt(rho^2) * max(sqrt(Q), 1e-12);
    
    % yに関して積分済みの関数を定義 (x -> -Inf to -S1 の範囲)
    result1 = integral(@(y) integratedXLower(y, S1, dycov, invcov), -yr, -S2, ... 
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result2 = integral(@(y) integratedXLower(y, S1, dycov, invcov), -S2, S2, ...
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result3 = integral(@(y) integratedXLower(y, S1, dycov, invcov), S2, yr, ...
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
        
    % -S1 to S1 の範囲
    result4 = integral(@(y) integratedXMiddle(y, S1, S1, dycov, invcov), -yr, -S2, ...
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result5 = integral(@(y) integratedXMiddle(y, S1, S1, dycov, invcov), -S2, S2, ... 
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result6 = integral(@(y) integratedXMiddle(y, S1, S1, dycov, invcov), S2, yr, ...
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
        
    % S1 to Inf の範囲
    result7 = integral(@(y) integratedXUpper(y, S1, dycov, invcov), -yr, -S2, ... 
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result8 = integral(@(y) integratedXUpper(y, S1, dycov, invcov), -S2, S2, ...
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result9 = integral(@(y) integratedXUpper(y, S1, dycov, invcov), S2, yr, ...
        'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
        
    result = result1 + result2 + result3 + result4 + result5 + ... 
             result6 + result7 + result8 + result9;
end

function result = integratedXLower(y, x1, dycov, invcov)
    % -Inf から -x1 までの x に関する積分結果
    global S2
    
    % 2次形式の係数
    a = invcov(1,1);
    b = invcov(1,2);
    c = invcov(2,2);
    
    detCov = det(dycov);
    
    % yの値に応じて係数を決定
    if y <= -S2
        coef = -S2;
    elseif y >= S2
        coef = S2;
    else
        coef = y;
    end
    
    % -S1 に対する積分結果
    erf_arg = sqrt(a) * (-x1 + b*y/a) / sqrt(2);
    
    if isreal(erf_arg) && all(isfinite(erf_arg))
        result = (-x1) * coef / (2*pi*sqrt(detCov)) * ... 
            sqrt(pi./a./2) .* exp(-(c - b.^2./a).*y.^2./2) .* ...
            (1 + erf(erf_arg));
    else
        result = 0;
    end
end

function result = integratedXMiddle(y, x1, x2, dycov, invcov)
    % -x1 から x1 までの x に関する積分結果
    global S2
    
    % 2次形式の係数
    a = invcov(1,1);
    b = invcov(1,2);
    c = invcov(2,2);
    
    detCov = det(dycov);
    
    % yの値に応じて係数を決定
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

function result = integratedXUpper(y, x2, dycov, invcov)
    % x2 (S1) から Inf までの x に関する積分結果
    global S2 S1
    
    % 2次形式の係数
    a = invcov(1,1);
    b = invcov(1,2);
    c = invcov(2,2);
    
    detCov = det(dycov);
    
    % yの値に応じて係数を決定
    if y <= -S2
        coef = -S2;
    elseif y >= S2
        coef = S2;
    else
        coef = y;
    end
    
    % S1 に対する積分結果
    erf_arg = sqrt(a) * (S1 + b*y/a) / sqrt(2);
    
    if isreal(erf_arg) && all(isfinite(erf_arg))
        result = S1 * coef / (2*pi*sqrt(detCov)) * ...
            sqrt(pi./a./2) .* exp(-(c - b.^2./a).*y.^2./2) .* ...
            (1 - erf(erf_arg));
    else
        result = 0;
    end
end