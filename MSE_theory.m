clc;
clear all;

tEnd = 100;

tic; % ストップウォッチタイマー開始

global fd2 dfd PRS RNG mu S1 S2 xi rho sgm_g

rho = 1;
S1 = 2;
S2 = 1;
xi = 0;
mu = 0.1;
sgm_g = sqrt(1);

PRS = 1e-10; % 数値積分の精度 RelTol, AbsTol
RNG = 7; % 数値積分の範囲を決めるパラメータ

% 出力フォルダの作成
outputDir = 'output';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

fname = char(['dyhouwa2MSE', ', tEnd=', num2str(tEnd), ', mu=', num2str(mu), ...
    ', S1=', num2str(S1), ', S2=', num2str(S2), '.txt']);
fname = fullfile(outputDir, fname);  % outputフォルダ内のパスを生成
Fid = fopen(fname, 'w');
Header = '#time Q r MSE';
fprintf(Fid, '%s\n', Header);

% 定数で表されるサンプル平均
fd2 = S1^2 - S1 * sqrt(2 * rho^2 * sgm_g^2 / pi) * ...
    exp(-(S1^2)/(2 * (rho^2) * sgm_g^2)) + ...
    (rho^2 * sgm_g^2 - S1^2) * erf(S1/sqrt(2 * rho^2 * sgm_g^2));

dfd = rho^2 * sgm_g^2 * erf(S1/sqrt(2 * rho^2 * sgm_g^2));

options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10); % 精度の調整
[t, y] = ode45(@(t, y) odefun(t, y), [0 tEnd], [1e-9, 1e-9], options);

Q = y(:, 1);
r = y(:, 2);

[sizet, ~] = size(t);
MSE = zeros(sizet, 1);

% 進行状況の表示のための変数
totalSteps = sizet; % ループの総ステップ数

for k = 1:sizet

    % 現在の進捗率を計算
    progress = (k / totalSteps) * 100;
    
    % 進行状況を表示
    fprintf('Progress: %.2f%%\n', progress);

    dycov = rho^2 * [sgm_g^2 r(k); r(k) Q(k)];
    invcov = inv(dycov);
    time = t(k);

    fun1 = @(d,y) (-S1).*(-S2)./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*...
        invcov(2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun2 = @(d,y) (-S1).*y./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*...
        invcov(2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun3 = @(d,y) (-S1).*S2./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*...
        invcov(2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun4 = @(d,y) d.*(-S2)./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*...
        invcov(2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun5 = @(d,y) d.*y./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*invcov...
        (2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun6 = @(d,y) d.*S2./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*invcov...
        (2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun7 = @(d,y) S1.*(-S2)./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*...
        invcov(2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun8 = @(d,y) S1.*y./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*invcov...
        (2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun9 = @(d,y) S1.*S2./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*...
        invcov(2,1).*d.*y+invcov(2,2).*y.^2)/2);


    fy2 = S2^2 - S2 * sqrt(2 * rho^2 * Q(k) / pi) * ...
        exp(-S2^2/(2 * rho^2 * Q(k))) + (rho^2 * Q(k) - S2^2) * ...
        erf(S2/sqrt(2 * rho^2 * Q(k)));


    % 積分範囲に上限を設ける
    dr = RNG * sqrt(rho^2 * sgm_g^2);
    yr = RNG * sqrt(rho^2) * max(sqrt(Q(k)), 1);

    fdfy = integral2(fun1, -dr, -S1, -yr, -S2, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun2, -dr, -S1, -S2, S2, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun3, -dr, -S1, S2, yr, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun4, -S1, S1, -yr, -S2, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun5, -S1, S1, -S2, S2, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun6, -S1, S1, S2, yr, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun7, S1, dr, -yr, -S2, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun8, S1, dr, -S2, S2, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun9, S1, dr, S2, yr, 'RelTol', PRS, 'AbsTol', PRS);

    MSE(k) = fd2 + fy2 - 2 * fdfy + xi;

    fprintf(Fid, '%g\t%g\t%g\t%g\n', time, Q(k), r(k), MSE(k));
end

fclose(Fid); % 記録用ファイルのクローズ

disp(['ファイル出力完了: ', fname]);
toc;

% 微分方程式の関数
function dydt = odefun(t, y)
    global fd2 dfd PRS RNG mu S1 S2 xi rho sgm_g
    
    dydt = zeros(2, 1); % 行列ベクトルとして連立方程式を表示
    
    Q = y(1);
    r = y(2);
    
    fy2 = S2^2 - S2 * sqrt(2 * rho^2 * Q/pi) * ...
        exp(-S2^2/(2 * rho^2 * Q)) + (rho^2 * Q - S2^2) * erf(S2/sqrt(2 * rho^2 * Q));
        
    dfy = rho^2 * r * erf(S2/sqrt(2 * rho^2 * Q));
    
    yfd = rho^2 * r * erf(S1/sqrt(2 * rho^2 * sgm_g^2));
    
    yfy = rho^2 * Q * erf(S2/sqrt(2 * rho^2 * Q));
    
    dycov = rho^2 * [sgm_g^2 r; r Q];
    invcov = inv(dycov);
    
    fun1 = @(d,y) (-S1).*(-S2)./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*...
        invcov(2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun2 = @(d,y) (-S1).*y./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*...
        invcov(2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun3 = @(d,y) (-S1).*S2./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*...
        invcov(2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun4 = @(d,y) d.*(-S2)./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*...
        invcov(2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun5 = @(d,y) d.*y./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*invcov...
        (2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun6 = @(d,y) d.*S2./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*invcov...
        (2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun7 = @(d,y) S1.*(-S2)./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*...
        invcov(2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun8 = @(d,y) S1.*y./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*invcov...
        (2,1).*d.*y+invcov(2,2).*y.^2)/2);

    fun9 = @(d,y) S1.*S2./(2.*pi.*sqrt(det(dycov))).*exp(-(invcov(1,1).*d.^2+2*...
        invcov(2,1).*d.*y+invcov(2,2).*y.^2)/2);

    dr = RNG * sqrt(rho^2 * sgm_g^2);
    yr = RNG * sqrt(rho^2) * max(sqrt(Q), 1);
    
    fdfy = integral2(fun1, -dr, -S1, -yr, -S2, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun2, -dr, -S1, -S2, S2, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun3, -dr, -S1, S2, yr, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun4, -S1, S1, -yr, -S2, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun5, -S1, S1, -S2, S2, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun6, -S1, S1, S2, yr, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun7, S1, dr, -yr, -S2, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun8, S1, dr, -S2, S2, 'RelTol', PRS, 'AbsTol', PRS) + ...
           integral2(fun9, S1, dr, S2, yr, 'RelTol', PRS, 'AbsTol', PRS);

    dydt(1) = mu^2 * rho^2 * (fd2 + fy2 - 2 * fdfy + xi) + 2 * mu * (yfd - yfy);
    dydt(2) = mu * (dfd - dfy);
end