clc;
clear;

%% 参数设置

N_sc=128;      %系统子载波数（不包括直流载波）、number of subcarrierA
N_fft=128;            % FFT 长度
N_cp=8;             % 循环前缀长度、Cyclic prefix
N_symbo=N_fft+N_cp;        % 1个完整OFDM符号长度
%N_c=53;             % 包含直流载波的总的子载波数、number of carriers
M=4;               %4PSK调制  改成16qam
SNR=0:5:30;         %仿真信噪比
N_frm=1;           %仿真时 仿真帧数
%N_frm=53;            % 每种信噪比下的仿真帧数、frame
%Nd=19040;               % 每帧包含的bits数
Nd=N_symbo*20;               % 每帧包含的bits数
sample_rate = 1.92e6; %采样率
P_t_inter=8;      %时域导频间隔
P_f_inter=4;      %频域导频间隔
data_station=[];    %导频位置
L=7;                %卷积码约束长度
tblen=6*L;          %Viterbi译码器回溯深度
stage = 3;          % m序列的阶数
ptap1 = [1 3];      % m序列的寄存器连接方式
regi1 = [1 1 1];    % m序列的寄存器初始值
N = Nd*N_frm;       %时域bit数

%% 基带数据数据产生
P_data=randi([0 1],1,N_sc*Nd*N_frm);


%% 信道编码（卷积码、或交织器）
%卷积码：前向纠错非线性码
%交织：使突发错误最大限度的分散化
trellis = poly2trellis(7,[133 171]);       %(2,1,7)卷积编码
code_data=convenc(P_data,trellis);


%% qpsk调制
data_temp1= reshape(code_data,log2(M),[])';             %以每组2比特进行分组，M=4
data_temp2= bi2de(data_temp1);                             %二进制转化为十进制
modu_data=pskmod(data_temp2,M,pi/M);              % 4PSK调制
% figure(1);
scatterplot(modu_data),grid;                  %星座图(也可以取实部用plot函数)

%% 扩频
%————————————————————————————————————————————————————————%
%扩频通信信号所占有的频带宽度远大于所传信息必需的最小带宽
%根据香农定理，扩频通信就是用宽带传输技术来换取信噪比上的好处，这就是扩频通信的基本思想和理论依据。
%扩频就是将一系列正交的码字与基带调制信号内积
%扩频后数字频率变成了原来的m倍。码片数量 = 2（符号数）* m（扩频系数）
%————————————————————————————————————————————————————————%

code = mseq(stage,ptap1,regi1,N_sc);     % 扩频码的生成
code = code * 2 - 1;         %将1、0变换为1、-1
modu_data=reshape(modu_data,N_sc,length(modu_data)/N_sc);
spread_data = spread(modu_data,code);        % 扩频
spread_data=reshape(spread_data,[],1);

%% 插入导频(C)

%P_f_inter:频域导频间隔，P_t_inter:时域导频间隔

 P_f=3+3*1i;                       %Pilot frequency
P_f_station=[1:P_f_inter:N_fft];%频域导频位置
Length_P_f=length(P_f_station);
Length_P_t=floor(length(spread_data)/(N_fft*P_t_inter-Length_P_f));%时域导频的个数
remain_P=mod(length(spread_data),(N_fft*P_t_inter-Length_P_f));%多出来的数据
Length_P_t1=ceil(length(spread_data)/(N_fft*P_t_inter-Length_P_f));%矩阵时域的最大的长度
P_t_station=[1:P_t_inter:(Length_P_t-1)*P_t_inter+1]; %时域导频位置
%构建数据矩阵
if (remain_P==0)
    data=zeros(N_fft,Length_P_t*P_t_inter);
else if remain_P<N_fft-Length_P_f+1
    data=zeros(N_fft,Length_P_t*P_t_inter+1);
    else
        r=ceil((remain_P-N_fft+Length_P_f)/N_fft);
        data=zeros(N_fft,Length_P_t*P_t_inter+1+r);
    end
end

%将导频放入矩阵
for i=1:Length_P_f
    for j=1:Length_P_t
        data(sub2ind(size(data),P_f_station(i),P_t_station(j))) = P_f;
    end
end
[data_row,data_col] = size(data);
data = reshape(data,[],1);
ind = find(~real(data));
add_data = zeros(data_row*data_col-length(spread_data)-Length_P_t*Length_P_f,1);
sp_data = [spread_data;add_data];
data(ind(1:end),1)=sp_data;

%% 串并转换
data=reshape(data,data_row,data_col);

%% IFFT
ifft_data=ifft(data)*sqrt(N_fft); 

%% 插入保护间隔、循环前缀
Tx_cd=[ifft_data(N_fft-N_cp+1:end,:);ifft_data];%把ifft的末尾N_cp个数补充到最前面

%% 并串转换
Tx_data=reshape(Tx_cd,[],1);%由于传输需要

%% 信道（通过多经瑞利信道+AWGN信道）
    
% 'PathDelays',[0 30 150 310 370 710 1090 1730 2510]*1e-9, ...
%     'AveragePathGains',[0.0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9], ...

%  'PathDelays',[0 1 2 3 4 5]/sample_rate, ...
%     'AveragePathGains',[0.0 -1.5  -3.6  -7.0 -12.0 -16.9], ...

%多径信道
 RayleighMultiPath =comm.RayleighChannel(...
    'SampleRate',sample_rate, ...
    'PathDelays',[0 30 150 310 370 710 1090 1730 2510]*1e-9, ...
    'AveragePathGains',[0.0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9], ...
    'NormalizePathGains',true, ...
    'MaximumDopplerShift',111, ...
     'DopplerSpectrum',doppler('Jakes'), ...
    'RandomStream','mt19937ar with seed', ...
    'Seed',22, ...
    'PathGainsOutputPort',true);
    [n1,h_pathGains]=RayleighMultiPath(Tx_data);%经过多径信道和多普勒频移

 
 %准确信道(多径)
%PathDelays=[0 1 2 3 4 5]/sample_rate;
 PathDelays=[0 30 150 310 370 710 1090 1730 2510]*1e-9;
 %PathDelays=[0];
 k=0:1:N_fft-1;
 k=repmat(k,data_col,1);
 k=k';
 H_real=0;
 for i = 1:length(PathDelays)
     h_i=h_pathGains(:,i);
     h_i=reshape(h_i,N_fft+N_cp,[]);
     h_i=h_i(N_cp+1:end,:);
     t_i=PathDelays(1,i);
     H_real = H_real+h_i.*exp((-2*pi*1i)*k*t_i*sample_rate/N_fft);
 end
 
 %高斯信道
 S = RandStream('mt19937ar','Seed',5489);
 Ber=zeros(1,length(SNR));
 Ber2=zeros(1,length(SNR));
for jj=1:length(SNR)
    rx_channel=awgn(n1,SNR(jj),'measured',S);%添加高斯白噪声

%% 串并转换
    %舍去信道延时，在矩阵后面补零
      rx_channel_dp = rx_channel(8:end);
      rx_channel_b0 = [rx_channel_dp;zeros(7,1)];
    Rx_data1=reshape(rx_channel_b0,N_fft+N_cp,[]);
    
%% 去掉保护间隔、循环前缀
    Rx_data2=Rx_data1(N_cp+1:end,:);
%% FFT
    fft_data=fft(Rx_data2)/sqrt(N_fft);
    
 %% 信道估计与插值（均衡）

      data3=fft_data(1:N_fft,:); 
      Rx_pilot=data3(P_f_station(1:end),P_t_station(1:end)); %接收到的导频
      pilot_seq = data(P_f_station(1:end),P_t_station(1:end));%发送的导频
      h=Rx_pilot./pilot_seq; 
      
      loc_f = 1:N_fft;
      loc_t = 1:data_col;
      %先插频域
     H_f=interp1( P_f_station(1:end),h,loc_f,'spline');%分段线性插值：插值点处函数值由连接其最邻近的两侧点的线性函数预测。对超出已知点集的插值点用指定插值方法计算函数值
     %再插时域
     H_t=interp1( P_t_station(1:end)',H_f',loc_t','spline');%分段线性插值：插值点处函数值由连接其最邻近的两侧点的线性函数预测。对超出已知点集的插值点用指定插值方法计算函数值
     H = H_t';
%     ERRA=abs((h-H_real(P_f_station(1:end),:))./H_real(P_f_station(1:end),:));
%     MSEA(jj)=sum(sum(ERRA.^2))/(length(P_f_station)*data_col);
    
%% 信道校正
      data_s = reshape(data3,[],1);
      H_s = reshape(H,[],1);
      data_aftereq=data_s(ind(1:end))./H_s(ind(1:end));
%% 并串转换
     data_aftereq=data_aftereq(1:length(spread_data));
     data_aftereq=reshape(data_aftereq,N_sc,length(data_aftereq)/N_sc);
    
%% 解扩
    demspread_data = despread(data_aftereq,code);       % 数据解扩
    
%% QPSK解调
    demodulation_data=pskdemod(demspread_data,M,pi/M);    
    De_data1 = reshape(demodulation_data,[],1);
    De_data2 = de2bi(De_data1);
    De_Bit = reshape(De_data2',1,[]);

%% （解交织）
%% 信道译码（维特比译码）
    trellis = poly2trellis(7,[133 171]);
    rx_c_de = vitdec(De_Bit,trellis,tblen,'trunc','hard');   %硬判决

%% 计算误码率
    [err,Ber2(jj)] = biterr(De_Bit(1:length(code_data)),code_data);%译码前的误码率
    [err, Ber(jj)] = biterr(rx_c_de(1:length(P_data)),P_data);%译码后的误码率

end
 figure(2);
 semilogy(SNR,Ber2,'b-s');
 hold on;
 semilogy(SNR,Ber,'r-o');
 hold on;
 legend('4PSK调制、卷积码译码前（有扩频）','4PSK调制、卷积码译码后（有扩频）');
 hold on;
 xlabel('SNR');
 ylabel('BER');
 title('AWGN信道下误比特率曲线');

 figure(3)
 subplot(2,1,1);
 x=0:1:30;
 stem(x,P_data(1:31));
 ylabel('amplitude');
 title('发送数据（以前30个数据为例)');
 legend('4PSK调制、卷积译码、有扩频');

 subplot(2,1,2);
 x=0:1:30;
 stem(x,rx_c_de(1:31));
 ylabel('amplitude');
 title('接收数据(以前30个数据为例)');
 legend('4PSK调制、卷积译码、有扩频');
 



