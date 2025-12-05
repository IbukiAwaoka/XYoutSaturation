clc;
clear all;

t = 5000;
N = 200;
ens = 100;
mu = 0.1; % ステップサイズ
xi = 0.1; % 背景雑音
S2 = 2;
S1_values = [S2*0.3, S2*0.7, S2, S2*1.2, S2*1.6, S2*2, S2*2.4, S2*2.8]; % 任意のSの値を指定0.3 , 0.7 , 1.0 , 1.2 , 1.6 , 2.0 , 2.4 , 2.8
num_S1_values = numel(S1_values);
MSE_values = zeros(num_S1_values, 1); % MSEの値を保存するベクトル
MSE_std_p = zeros(num_S1_values, 1); % MSEの標準偏差を保存するベクトル
MSE_std_n = zeros(num_S1_values, 1);

tic; % ストップウォッチタイマー開始
for k = 1:num_S1_values
    S1 = S1_values(k);
    n = N * t + 1; % 更新回数
    MSE = zeros(ens, 1); % MSE 初期化
    for i = 1:ens
        g = randn(N, 1); % 未知システム
        w = zeros(N, 1); % 適応フィルタの初期化
        x = randn(N, 1) / sqrt(N); % N行1列の入力ベクトル
        for j = 1:n
            x = circshift(x, 1); % 入力ベクトルをすべてずらす
            x(1) = randn / sqrt(N); % 値を上書きする
            d = g.' * x; % 未知システムの出力
            y = w.' * x; % 適応フィルタの出力
            if d < -S1
                d = -S1;
            elseif d > S1
                d = S1;
            else
                d = d;
            end
 
            if y < -S2
                y = -S2;
            elseif y > S2
                y = S2;
            else
                y = y;                                                      
            end
            e = d - y + sqrt(xi) * randn; % 誤差の生成

            if j >= n - 9
                MSE(i) = MSE(i) + e^2; % 最後の10回分の誤差の二乗を足す
            end

            w = w + mu * e * x; % 適応フィルタの更新
        end
        MSE(i) = MSE(i) / 10; % 最後の10回分の平均を計算
    end
    MSE_values(k) = median(MSE);  % MSEの中央値
    MSE_descend = sort(MSE, 'descend'); % MSEを降順にソート
    MSE_std_p(k) = MSE_descend(round(0.1587 * ens)); % 標準偏差の近似(上位16％)
    MSE_std_n(k) = MSE_descend(round((1-0.1587) * ens));% 標準偏差の近似(下位16％)
    disp([num2str(k * 100 / num_S1_values), '%']); % 進行状況の確認
end

% 結果の保存
data = [S1_values/S2; MSE_values'; MSE_std_p'; MSE_std_n'];
header = 'dyhouwaerrorbar';
fname = char([header, '(N=', num2str(N), ',S2=', num2str(S2), ',T=', num2str(t), ',mu=', num2str(mu), ',ens=', num2str(ens), ',xi=', num2str(xi), '_errorbars).txt']);
Fid = fopen(fname, 'w');
fprintf(Fid, 'S1/S2\tMSE\tMSE_std_p\tMSE_std_n\n');
fprintf(Fid, '%.2f\t%.5f\t%.5f\t%.5f\n', data);
fclose(Fid);
toc; % タイマー終了
