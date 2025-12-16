clc;
clear all;

tEnd = 300;
tic;

% パラメータ設定
rho = 1;
S1 = 0.5;
xi = 0;
mu = 1;
sgm_g = sqrt(1);
PRS = 1e-10;
RNG = 8;

% 出力フォルダの作成
outputDir = 'Ioutput';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% S2/S1の範囲設定
S2_S1_ratio = linspace(0, 3, 901);
num_ratios = length(S2_S1_ratio);

% 結果保存用の配列
MSE_results = zeros(num_ratios, 1);
Q_results = zeros(num_ratios, 1);
r_results = zeros(num_ratios, 1);

fprintf('=== S2スイープ計算開始（並列計算版 - parfeval） ===\n');
fprintf('S1 = %.2f (固定)\n', S1);
fprintf('t = %d\n', tEnd);
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
    S2 = S2_S1_ratio(idx) * S1;
    f(idx) = parfeval(@computeMSEforS2, 3, S2, S1, tEnd, rho, xi, mu, sgm_g, PRS, RNG);
end

% 結果を順次取得しながら進捗表示
for idx = 1:num_ratios
    [completedIdx, MSE_final, Q_final, r_final] = fetchNext(f);
    MSE_results(completedIdx) = MSE_final;
    Q_results(completedIdx) = Q_final;
    r_results(completedIdx) = r_final;
    
    fprintf('完了: %d/%d (%.2f%%) - S2/S1=%.3f, MSE=%.6e\n', idx, num_ratios, ...
        (idx/num_ratios)*100, S2_S1_ratio(completedIdx), MSE_final);
end

% 結果をファイルに保存
fname = char(['parfeval_MSE_S2sweep_t=', num2str(tEnd), '_mu=', num2str(mu), ...
    '_S1=', num2str(S1), '.txt']);
fname = fullfile(outputDir, fname);
Fid = fopen(fname, 'w');
Header = '#S2/S1 S2 MSE(t=300) Q(t=300) r(t=300)';
fprintf(Fid, '%s\n', Header);

for idx = 1:num_ratios
    fprintf(Fid, '%.6f\t%.6f\t%.10e\t%.10e\t%.10e\n', ...
        S2_S1_ratio(idx), S2_S1_ratio(idx)*S1, MSE_results(idx), Q_results(idx), r_results(idx));
end
fclose(Fid);

fprintf('\n=== 計算完了 ===\n');
disp(['結果ファイル: ', fname]);
toc;

% プロット
figure('Position', [100, 100, 800, 600]);
plot(S2_S1_ratio, MSE_results, 'b-', 'LineWidth', 2);
grid on;
xlabel('S2/S1', 'FontSize', 14);
ylabel('MSE (t=300)', 'FontSize', 14);
title(sprintf('MSE vs S2/S1 at t=%d (S1=%.2f fixed)', tEnd, S1), 'FontSize', 16);

plot_fname = fullfile(outputDir, ['MSE_S2sweep_parfeval_t', num2str(tEnd), '_S1', num2str(S1), '.png']);
saveas(gcf, plot_fname);
disp(['プロット保存: ', plot_fname]);

%% 並列計算用の関数（各ワーカーで実行される）
function [MSE_final, Q_final, r_final] = computeMSEforS2(S2, S1, tEnd, rho, xi, mu, sgm_g, PRS, RNG)
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
