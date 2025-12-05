clc ;
clear all ;

t = 100 ;
N = 200 ;
ens = 1000 ;
a = 1 ;
xi = 0 ;   % 背景雑音
n = N*t+1;   % 更新回数
S1 = 1;   % 未知システム飽和値
S2 = 1;   % 適応フィルタ飽和値
mu = 1 ;   % ステップサイズ
MSD = zeros(a*t+1 , 1);   % MSD 初期化


tic;  %ストップウォッチタイマー開始

for i = 1 : ens

    g = randn(N,1);   % 未知システム
    w = zeros(N,1);   % 適応フィルタの初期化
    x = randn(N,1)/sqrt(N);   % N 行1 列の入力ベクトル

    for j = 1 : n

        x = circshift(x,1);   % 入力ベクトルをすべてずらす
        x(1) = randn/sqrt(N);   % 値を上書きする

        d = g.'*x;   % 未知システムの出力
        y = w.'*x;   % 適応フィルタの出力

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


        gw = g-w;    %　偏差の生成
        e = d - y + sqrt(xi) * randn;

        if(rem(j,N/a)==1)   % MSDの計算
            MSD((j-1)/(N/a)+1) = MSD((j-1)/(N/a)+1) + ((gw.'*gw)/N)/ens;
        end

        w = w + mu*e*x;   % 適応フィルタの更新

    end

    disp([num2str(i*100/ens),'%'])   % 進行状況の確認

end

data = [(0:1/a:t);MSD'];
header = 'MSD Simulation';
fname = char([header,'(N=',num2str(N),',T=',num2str(t),',S1=',num2str(S1),',S2=',num2str(S2),',mu=',num2str(mu),',ens=',num2str(ens),',xi=',num2str(xi),').txt']);
Fid = fopen(fname,'w');
S = '#time    MSD';   % 記録ファイルの一行目
fprintf(Fid,'%s\n',S);
fprintf(Fid,'%f   %1.5f\n',data);
fclose(Fid);
toc;   % タイマー終了
