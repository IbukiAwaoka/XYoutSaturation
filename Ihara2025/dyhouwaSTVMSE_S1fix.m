clc;
clear all;

% outputフォルダの作成
if ~exist('Ioutput', 'dir')
    mkdir('Ioutput');
end

tEnd = 300;

tic;    % ストップウォッチタイマー開始

global fd2 dfd PRS RNG mu S1 S2 xi rho sgm_g

rho = 1;
S1 = 1;
S2loop = 0:S1/100:S1*3; %初期値，刻み幅，終了値
xi = 0.1;
mu = 0.1;
sgm_g = sqrt(1);

PRS = 1e-6;  % 数値積分の精度RelTol,AbsTol
RNG = 7;  % 数値積分の範囲を決めるパラメータ

fname = char(['dyhouwa1intS2loop', ',tEnd=', num2str(tEnd), ',mu=', num2str(mu), ',S1=', num2str(S1), ',S2=', num2str(S2), ',xi=', num2str(xi),'.txt']);
Fid = fopen(fname, 'w');
Header = '#S2/S1          STVMSE';
fprintf(Fid, '%s\n', Header);

% 定数で表されるサンプル平均
fd2 = S1^2 - S1 * sqrt(2 * rho^2 * sgm_g^2 / pi) * exp(-(S1^2) / (2 * (rho^2) * sgm_g^2)) + (rho^2 * sgm_g^2 - S1^2) * erf(S1 / sqrt(2 * rho^2 * sgm_g^2));
dfd = rho^2 * sgm_g^2 * erf(S1 / sqrt(2 * rho^2 * sgm_g^2));

for i = 1:length(S2loop)
    S2 = S2loop(i);
    fprintf('Progress for S2 = %.2f: %.2f%%\n', S2, (i / length(S2loop)) * 100);
    
    % ODEソルバーで時間発展を計算
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);  % 精度の調整
    [t, y] = ode45(@(t, y) odefun(t, y), [0 tEnd], [1e-9, 1e-9], options);
    
    % 最終時刻のQ, rを取得
    Q_final = y(end, 1);
    r_final = y(end, 2);
    
    % 最終時刻での共分散行列とその逆行列を計算
    dycov_final = rho^2 * [sgm_g^2, r_final; r_final, Q_final];
    invcov_final = inv(dycov_final);
    
    % 最終時刻でのMSEを計算
    MSE_final = calculateMSE(Q_final, r_final, dycov_final, invcov_final);
    
    % 結果をファイルに書き込み
    fprintf(Fid, '%g   %g\n', S2/S1, MSE_final);
end

fclose(Fid);
toc;

function dydt = odefun(t, y)
    global fd2 dfd mu S1 S2 xi rho sgm_g 
    dydt = zeros(2, 1);
    Q = max(y(1), 1e-12);  % Qが負値にならないように制限
    r = max(y(2), 1e-12);  % rが負値にならないように制限

erf_arg1 = S2 / sqrt(2 * rho^2 * Q);
erf_arg2 = S1 / sqrt(2 * rho^2 * sgm_g^2);

% erf_argが実数で有限な値かをチェック
if isreal(erf_arg1) && isfinite(erf_arg1)
    fy2 = S2^2 - S2 * sqrt(2 * rho^2 * Q / pi) * exp(-S2^2 / (2 * rho^2 * Q)) + (rho^2 * Q - S2^2) * erf(erf_arg1);
else
    % erf_argが無効な場合は、代わりに0
    fy2 = 0;  % erfの代わりに0
end

if isreal(erf_arg1) && isfinite(erf_arg1)
    dfy = rho^2 * r * erf(erf_arg1);
else
    dfy =  0;  % erfの代わりに0
end

if isreal(erf_arg2) && isfinite(erf_arg2)
    yfd = rho^2 * r * erf(erf_arg2);
else
    yfd =  0;  % erfの代わりに0
end

if isreal(erf_arg1) && isfinite(erf_arg1)
    yfy = rho^2 * Q * erf(erf_arg1);
else
    yfy =  0;  % erfの代わりに0
end


    dycov = rho.^2 .* [sgm_g.^2 r; r Q];
    

    invcov = inv(dycov);

    
    fdfy = calculateIntegral(Q, dycov, invcov);
    
    dydt(1) = mu^2 * rho^2 * (fd2+fy2-2*fdfy+xi) + 2*mu*(yfd-yfy);
    dydt(2) = mu * (dfd - dfy);
end

function result = calculateMSE(Q, r, dycov, invcov)
    global  S2 rho xi fd2
    
    erf_arg1 = S2 / sqrt(2 * rho^2 * Q);
if isreal(erf_arg1) && isfinite(erf_arg1)
    fy2 = S2^2 - S2 * sqrt(2 * rho^2 * Q / pi) * exp(-S2^2 / (2 * rho^2 * Q)) + (rho^2 * Q - S2^2) * erf(erf_arg1);
else
    % erf_argが無効な場合は、代わりに0
    fy2 = 0;  % erfの代わりに0
end
    fdfy = calculateIntegral(Q, dycov, invcov);
    
    result = fd2 + fy2 - 2 * fdfy + xi;
end

function result = calculateIntegral(Q, dycov, invcov)
    global S1 S2 PRS RNG rho
    
    % y に関する積分範囲
    yr = RNG * sqrt(rho^2) * max(sqrt(Q), 1e-9);
    
    % xに関して積分済みの関数を定義（-Inf to -S1 の範囲）
    result1 = integral(@(y) integratedXLower(y, -S1, dycov, invcov), -yr, -S2, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result2 = integral(@(y) integratedXLower(y, -S1, dycov, invcov), -S2, S2, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result3 = integral(@(y) integratedXLower(y, -S1, dycov, invcov), S2, yr, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    
    % -S1 to S1 の範囲
    result4 = integral(@(y) integratedXMiddle(y, -S1, S1, dycov, invcov), -yr, -S2, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result5 = integral(@(y) integratedXMiddle(y, -S1, S1, dycov, invcov), -S2, S2, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result6 = integral(@(y) integratedXMiddle(y, -S1, S1, dycov, invcov), S2, yr, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    
    % S1 to Inf の範囲
    result7 = integral(@(y) integratedXUpper(y, S1, dycov, invcov), -yr, -S2, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result8 = integral(@(y) integratedXUpper(y, S1, dycov, invcov), -S2, S2, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    result9 = integral(@(y) integratedXUpper(y, S1, dycov, invcov), S2, yr, 'RelTol', PRS, 'AbsTol', PRS, 'ArrayValued', true);
    
    result = result1 + result2 + result3 + result4 + result5 + result6 + result7 + result8 + result9;
end

function result = integratedXLower(y, x1, dycov, invcov)
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

function result = integratedXMiddle(y, x1, x2, dycov, invcov)
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

function result = integratedXUpper(y, x2, dycov, invcov)
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
