clc;
clear all;

% 並列プールの開始
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool;
end

% Ioutputフォルダの作成（親ディレクトリ）
if ~exist('Ioutput', 'dir')
    mkdir('Ioutput');
end

tEnd = 300;

tic;    % ストップウォッチタイマー開始

global PRS RNG mu S1 S2 xi rho sgm_g

rho = 1;
S1loop = 0:0.1:3;
S2loop = 0:0.1:3; 
xi = 0.1;
mu = 0.1;
sgm_g = sqrt(1);

fname = char(['dyhouwaSTVMSETheoryheatmap',  ',xi=', num2str(xi),'.txt']);
fname = fullfile('Ioutput', fname);  % Ioutputフォルダ内のパスを生成

% 結果を格納する行列
results = zeros(length(S1loop), length(S2loop));

parfor i = 1:length(S1loop)
    S1 = S1loop(i);
    local_results = zeros(1, length(S2loop));

    for j = 1:length(S2loop)
        S2 = S2loop(j);
        fprintf('Progress: S1 = %.2f, S2 = %.2f\n', S1, S2);
        
        % ODEソルバーで時間発展を計算
        options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);  % 精度の調整
        [t, y] = ode45(@(t, y) odefun(t, y, mu, S1, S2, xi, rho, sgm_g), [0 tEnd], [1e-9, 1e-9], options);
        
        % 最終時刻のQ, rを取得
        Q_final = y(end, 1);
        r_final = y(end, 2);
        
        % 最終時刻での共分散行列とその逆行列を計算
        dycov_final = rho^2 * [sgm_g^2, r_final; r_final, Q_final];
        invcov_final = inv(dycov_final);
        
        % 最終時刻でのMSEを計算
        MSE_final = calculateMSE(Q_final, r_final, dycov_final, invcov_final, S1, S2, rho, xi, sgm_g);
        
        local_results(j) = MSE_final;
    end
    results(i, :) = local_results;
end

% 結果をファイルに書き込み
Fid = fopen(fname, 'w');
Header = '#S1         S2       STVMSE';
fprintf(Fid, '%s\n', Header);

for i = 1:length(S1loop)
    S1 = S1loop(i);
    for j = 1:length(S2loop)
        S2 = S2loop(j);
        fprintf(Fid, '%g   %g   %g\n', S1, S2, results(i, j));
    end
    fprintf(Fid, '\n');
end

fclose(Fid);
toc;

% 並列プールのクローズ
delete(gcp('nocreate'));

function dydt = odefun(t, y, mu, S1, S2, xi, rho, sgm_g)
    dydt = zeros(2, 1);
    Q = max(y(1), 1e-12);  % Qが負値にならないように制限
    r = max(y(2), 1e-12);  % rが負値にならないように制限

erf_arg1 = S2 / sqrt(2 * rho^2 * Q);
erf_arg2 = S1 / sqrt(2 * rho^2 * sgm_g^2);

% erf_argが実数で有限な値かをチェック
if isreal(erf_arg1) && isfinite(erf_arg1)
    fy2 = S2^2 - S2 * sqrt(2 * rho^2 * Q / pi) * exp(-S2^2 / (2 * rho^2 * Q)) + (rho^2 * Q - S2^2) * erf(erf_arg1);
    dfy = rho^2 * r * erf(erf_arg1);
    yfy = rho^2 * Q * erf(erf_arg1);
else
    fy2 = 0;  % erfの代わりに0
    dfy =  0;
    yfy =  0;
end


if isreal(erf_arg2) && isfinite(erf_arg2)
    yfd = rho^2 * r * erf(erf_arg2);
    fd2 = S1^2 - S1 * sqrt(2 * rho^2 * sgm_g^2 / pi) * exp(-(S1^2) / (2 * (rho^2) * sgm_g^2)) + (rho^2 * sgm_g^2 - S1^2) * erf(erf_arg2);
    dfd = rho^2 * sgm_g^2 * erf(erf_arg2);
else
    yfd =  0;  % erfの代わりに0
    fd2 = 0;
    dfd = 0;
end

    dycov = rho.^2 .* [sgm_g.^2 r; r Q];
    
    invcov = inv(dycov);

    fdfy = calculateIntegral(Q, dycov, invcov, S1, S2, rho);
    
    dydt(1) = mu^2 * rho^2 * (fd2+fy2-2*fdfy+xi) + 2*mu*(yfd-yfy);
    dydt(2) = mu * (dfd - dfy);
end

function result = calculateMSE(Q, r, dycov, invcov, S1, S2, rho, xi, sgm_g)
    global  S2 rho xi S1 sgm_g
    
    erf_arg1 = S2 / sqrt(2 * rho^2 * Q);
    erf_arg2 = S1 / sqrt(2 * rho^2 * sgm_g^2);

if isreal(erf_arg1) && isfinite(erf_arg1)
    fy2 = S2^2 - S2 * sqrt(2 * rho^2 * Q / pi) * exp(-S2^2 / (2 * rho^2 * Q)) + (rho^2 * Q - S2^2) * erf(erf_arg1);
else
    % erf_argが無効な場合は、代わりに0
    fy2 = 0;  % erfの代わりに0
end

if isreal(erf_arg2) && isfinite(erf_arg2)
    fd2 = S1^2 - S1 * sqrt(2 * rho^2 * sgm_g^2 / pi) * exp(-(S1^2) / (2 * (rho^2) * sgm_g^2)) + (rho^2 * sgm_g^2 - S1^2) * erf(erf_arg2);
else
    fd2 = 0;
end

    fdfy = calculateIntegral(Q, dycov, invcov, S1, S2, rho);
    
    result = fd2 + fy2 - 2 * fdfy + xi;

    if isnan(result) || result < 0
        result = 0;
    end
end

function result = calculateIntegral(Q, dycov, invcov, S1, S2, rho)
    PRS = 1e-6;
    RNG = 7;
    
    % y に関する積分範囲
    yr = RNG * sqrt(rho^2) * max(sqrt(Q), 1e-9);
    
    % xに関して積分済みの関数を定義（-Inf to -S1 の範囲）
    result1 = integral(@(y) integratedXLower(y, -S1, dycov, invcov, S2, S1), -yr, -S2, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result2 = integral(@(y) integratedXLower(y, -S1, dycov, invcov, S2, S1), -S2, S2, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result3 = integral(@(y) integratedXLower(y, -S1, dycov, invcov, S2, S1), S2, yr, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    
    % -S1 to S1 の範囲
    result4 = integral(@(y) integratedXMiddle(y, -S1, S1, dycov, invcov, S2), -yr, -S2, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result5 = integral(@(y) integratedXMiddle(y, -S1, S1, dycov, invcov, S2), -S2, S2, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result6 = integral(@(y) integratedXMiddle(y, -S1, S1, dycov, invcov, S2), S2, yr, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    
    % S1 to Inf の範囲
    result7 = integral(@(y) integratedXUpper(y, S1, dycov, invcov, S2, S1), -yr, -S2, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result8 = integral(@(y) integratedXUpper(y, S1, dycov, invcov, S2, S1), -S2, S2, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result9 = integral(@(y) integratedXUpper(y, S1, dycov, invcov, S2, S1), S2, yr, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    
    result = result1 + result2 + result3 + result4 + result5 + result6 + result7 + result8 + result9;
end

function result = integratedXLower(y, x1, dycov, invcov, S2, S1)
    % -Inf から x1 までの x に関する積分結果
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
    
    % -S1 に対する積分結果
erf_arg = sqrt(a) * (x1 + b * y / a) / sqrt(2);
if isreal(erf_arg) && all(isfinite(erf_arg))  % 実数かつ有限な値であることを確認
    result = (-S1) .* coef ./ (2.*pi.*sqrt(detCov)) .* ...
             sqrt(pi./a./2) .* exp(-(c - b.^2./a).*y.^2./2) .* ...
             (1 + erf(erf_arg));
else
    % erfの引数が無効な場合
    result = 0;
end

end

function result = integratedXMiddle(y, x1, x2, dycov, invcov, S2)
    % x1 から x2 までの x に関する積分結果
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
    
erf_arg1 = sqrt(a./2) * (x2 + b * y ./ a);
    erf_arg2 = sqrt(a./2) * (x1 + b * y ./ a);
    
    % erfの引数が実数かつ有限な値であることを確認
    if isreal(erf_arg1) && all(isfinite(erf_arg1)) && isreal(erf_arg2) && all(isfinite(erf_arg2))
        result = coef ./ (2.*pi.*sqrt(detCov)) .* exp(-(c - b.^2./a).*y.^2./2) .* ...
            ((-1./a) .* (exp(-a*(x2 + b.*y./a).^2./2) - exp(-a.*(x1 + b.*y./a).^2./2)) - ...
            (b.*y./a) .* sqrt(pi./2./a) .* (erf(erf_arg1) - erf(erf_arg2)));
    else
        result = 0;  % erfの引数が無効な場合、0を返す
    end

end

function result = integratedXUpper(y, x2, dycov, invcov, S2, S1)
    % x2 から Inf までの x に関する積分結果
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
    erf_arg = sqrt(a) * (x2 + b * y / a) / sqrt(2);
    
    % erfの引数が実数かつ有限な値であることを確認
    if isreal(erf_arg) && all(isfinite(erf_arg))
        result = S1 .* coef ./ (2.*pi.*sqrt(detCov)) .* ...
                 sqrt(pi./a./2) .* exp(-(c - b.^2./a).*y.^2./2) .* ...
                 (1 - erf(erf_arg));
    else
        result = 0;  % erfの引数が無効な場合、0を返す
    end
end