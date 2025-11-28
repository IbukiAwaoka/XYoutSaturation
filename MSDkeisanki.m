clc;
clear all;

t = input('時刻を入力してください（標準値50）');
N = input('次元を入力してください（標準値200）');
ens = input('試行回数を入力してください');
a = 2;

tic;    % ストップウォッチタイマー開始

mu = 0.1;    % ステップサイズ
xi = 1;    % 背景雑音
n = N * t + 1;    % 更新回数
S1 = 1;    % 入力側飽和特性
S2 = 1;    % 出力側飽和特性
MSD = zeros(a * t + 1, 1);    % MSD初期化

for i = 1:ens

    g = randn(N, 1);   % 未知システム
    w = zeros(N, 1);   % 適応フィルタの初期化
    
    
    % 信号の時系列を表現する形で初期化します。
    u_history = randn(N, 1) / sqrt(N);
    % ***************************************

    for j = 1:n
        
        
        % 新しい入力信号を生成し、古い信号をシフトさせます。
        new_u_sample = randn / sqrt(N);
        u_history = circshift(u_history, 1);
        u_history(1) = new_u_sample;
        % ******************************************************

        x = g.' * u_history;    % 未知システムの出力

        % *** 修正3: 非線形飽和型関数をxに直接適用し、結果をx_hに格納 ***
        x_h = x;
        if x > S1
            x_h = S1;
        elseif x < -S1
            x_h = -S1;
        end
        % ****************************************************************

        y = w.' * u_history;    % 適応フィルタの出力

        y_h = y;
        if y > S2
            y_h = S2;
        elseif y < -S2
            y_h = -S2;
        end

        % *** 修正4: 誤差信号eの計算に非線形出力を利用 ***
        e = x_h - y_h + sqrt(xi) * randn;    % 誤差の生成
        % ****************************************************

        sd = g - w;    % 偏差の生成
        msd = ((sd).' * (sd));    % MSDの計算
        
        if rem(j, N / a) == 1    % MSDの計算
            MSD((j - 1) / (N / a) + 1) = MSD(((j - 1) / (N / a) + 1)) + (msd / N) / ens;
        end
        
        % *** 修正5: 適応フィルタの更新にu_historyを使用 ***
        w = w + mu * e * u_history;    % 適応フィルタの更新
        % ******************************************************
    end
    disp([num2str(i * 100 / ens), '%'])    % 進行状況の確認
end

data = [(0:1/a:t); MSD'];
header = 'SimulationMSD';
fname = char([header, '(N=', num2str(N), ',T=', num2str(t), ',S1=', num2str(S1), ',S2=', num2str(S2), ',mu=', num2str(mu), ',ens=', num2str(ens), ',xi=', num2str(xi), ',a=', num2str(a), ').txt']);
Fid = fopen(fname, 'w');
S = '#time     MSD';    % 記録ファイルの一行目
fprintf(Fid, '%s\n', S);
fprintf(Fid, '%f   %1.5f\n', data);
fclose(Fid);
toc;    % タイマー終了