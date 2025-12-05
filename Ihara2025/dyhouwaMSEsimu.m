clc ;
clear all ;

% outputフォルダの作成
if ~exist('Ioutput', 'dir')
    mkdir('Ioutput');
end

t = 5000; %input(' 時刻を入力してください(標準値 30)') ;
N = 200; %input(' 次元を入力してください(標準値 200)') ;
ens = 1000; %input(' 試行回数を入力してください ') ;
a = 2; %input(' データを記録する間隔を調整するパラメータを入力してください ') ;


tic;  %ストップウォッチタイマー開始

mu = 0.1 ;   % ステップサイズ
xi = 0 ;   % 背景雑音
n = N*t+1;   % 更新回数
S1 = 1;   % 未知システム飽和値
S2 = 1;   % 適応フィルタ飽和値
MSE = zeros(a*t+1 , 1);   % MSE 初期化

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

            if(rem(j,N/a)==1)   % MSEの計算
            MSE((j-1)/(N/a)+1) = MSE((j-1)/(N/a)+1) + e^2/ens;
            end
            
            w = w + mu * e * x; % 適応フィルタの更新
                 
        end

    disp([num2str(i*100/ens),'%'])   % 進行状況の確認

end

data = [(0:1/a:t);MSE'];
header = 'dyhouwaMSE Simulation';
fname = char([header,'(N=',num2str(N),',T=',num2str(t),',S1=',num2str(S1),',S2=',num2str(S2),',mu=',num2str(mu),',ens=',num2str(ens),',xi=',num2str(xi),').txt']);
Fid = fopen(fname,'w');
F = '#time    MSE';   % 記録ファイルの一行目


fprintf(Fid,'%s\n',F);
fprintf(Fid,'%f   %1.5f\n',data);
fclose(Fid);
toc;   % タイマー終了