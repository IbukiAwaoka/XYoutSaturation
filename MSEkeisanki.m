%未知システムの出力側に飽和型非線形要素があった場合のシミュレーション
clc;
clear all;

% ユーザー入力
fprintf('時刻を入力してください（標準値100）: ');
t = input('');
fprintf('次元を入力してください（標準値200）: ');
N = input('');
fprintf('試行回数を入力してください（標準値1000）: ');
ens = input('');
fprintf('メモリ幅の決定: ');
a = input('');

tic;    % ストップウォッチタイマー開始

% 定数の設定
mu = 0.3;    % ステップサイズ
xi = 0.0;      % 背景雑音
n = N * t + 1; % 更新回数
S1 = 1;     % 入力側飽和特性の幅
S2 = 1;     % 出力側飽和特性の幅
MSE = zeros(a * t + 1, 1); % MSEの初期化

for i = 1:ens

    % 初期化
    g = randn(N, 1);   % 未知システム
    w = zeros(N, 1);   % 適応フィルタ
    u = randn(N, 1) / sqrt(N); % N行1列の入力ベクトル

    for j = 1:n

        % 入力ベクトルの更新
        u = circshift(u, 1); % 入力ベクトルを全てずらす
        u(1) = randn / sqrt(N); % 新しいランダム値を追加
       
        % 未知システムの出力
        x = g.' * u;
        
        x_h=x;
         % 非線形飽和型関数を適用
        if x > S1
            x_h = S1;
        elseif x < -S1
            x_h = -S1;
        else
            x_h = x;
        end

        % 適応フィルタの出力
        y = w.' * u;

        y_h=y;
         % 非線形飽和型関数を適用
        if y > S2
            y_h = S2;
        elseif y < -S2
            y_h = -S2;
        else
            y_h = y;
        end

        % 誤差の生成
        e = x_h - y_h + sqrt(xi) * randn;
        % MSEの計算
        if rem(j, N / a) == 1
            MSE((j - 1) / (N / a) + 1) = MSE((j - 1) / (N / a) + 1) + e^2 / ens;
        end

        % 適応フィルタの更新
        w = w + mu * e * u;
    end

    % 進行状況の表示
    disp([num2str(i * 100 / ens), '%'])
end

% 結果の保存
% 出力フォルダの作成
outputDir = 'output';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

data = [(0:1 / a:t); MSE'];
header = 'MSE Simulation';
fname = char([header, '(N=', num2str(N), ',T=', num2str(t), ',S1=', num2str(S1), ',S2=', num2str(S2), ',mu=', num2str(mu), ',ens=', num2str(ens), ',xi=', num2str(xi), ',a=', num2str(a), ').txt']);
fname = fullfile(outputDir, fname);  % outputフォルダ内のパスを生成
Fid = fopen(fname, 'w');
S = '#time     MSE';    % 記録ファイルの一行目
fprintf(Fid, '%s\n', S);
fprintf(Fid, '%f   %1.5f\n', data);
fclose(Fid);

toc;    % タイマー終了
