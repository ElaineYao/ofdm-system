clc;
clear;

%% ��������

N_sc=128;      %ϵͳ���ز�����������ֱ���ز�����number of subcarrierA
N_fft=128;            % FFT ����
N_cp=8;             % ѭ��ǰ׺���ȡ�Cyclic prefix
N_symbo=N_fft+N_cp;        % 1������OFDM���ų���
%N_c=53;             % ����ֱ���ز����ܵ����ز�����number of carriers
M=4;               %4PSK����  �ĳ�16qam
SNR=0:5:30;         %���������
N_frm=1;           %����ʱ ����֡��
%N_frm=53;            % ÿ��������µķ���֡����frame
%Nd=19040;               % ÿ֡������bits��
Nd=N_symbo*2;               % ÿ֡������bits��
sample_rate = 1.92e6; %������
P_t_inter=32;      %ʱ��Ƶ���
P_f_inter=16;      %Ƶ��Ƶ���
data_station=[];    %��Ƶλ��
L=7;                %�����Լ������
tblen=6*L;          %Viterbi�������������
stage = 3;          % m���еĽ���
ptap1 = [1 3];      % m���еļĴ������ӷ�ʽ
regi1 = [1 1 1];    % m���еļĴ�����ʼֵ
N = Nd*N_frm;       %ʱ��bit��

%% �����������ݲ���
P_data=randi([0 1],1,N_sc*Nd*N_frm);


%% �ŵ����루����롢��֯����
%����룺ǰ������������
%��֯��ʹͻ����������޶ȵķ�ɢ��
trellis = poly2trellis(7,[133 171]);       %(2,1,7)�������
code_data=convenc(P_data,trellis);


%% qpsk����
data_temp1= reshape(code_data,log2(M),[])';             %��ÿ��2���ؽ��з��飬M=4
data_temp2= bi2de(data_temp1);                             %������ת��Ϊʮ����
modu_data=pskmod(data_temp2,M,pi/M);              % 4PSK����
% figure(1);
scatterplot(modu_data),grid;                  %����ͼ(Ҳ����ȡʵ����plot����)

%% ��Ƶ
%����������������������������������������������������������������������������������������������������������������%
%��Ƶͨ���ź���ռ�е�Ƶ�����Զ����������Ϣ�������С����
%������ũ������Ƶͨ�ž����ÿ�����似������ȡ������ϵĺô����������Ƶͨ�ŵĻ���˼����������ݡ�
%��Ƶ���ǽ�һϵ����������������������ź��ڻ�
%��Ƶ������Ƶ�ʱ����ԭ����m������Ƭ���� = 2����������* m����Ƶϵ����
%����������������������������������������������������������������������������������������������������������������%

code = mseq(stage,ptap1,regi1,N_sc);     % ��Ƶ�������
code = code * 2 - 1;         %��1��0�任Ϊ1��-1
modu_data=reshape(modu_data,N_sc,length(modu_data)/N_sc);
spread_data = spread(modu_data,code);        % ��Ƶ
spread_data=reshape(spread_data,[],1);

%% ���뵼Ƶ
P_f=3+3*1i;                       %Pilot frequency
P_f_station=[1:P_f_inter:N_fft];%��Ƶλ�ã���Ƶλ�ú���Ҫ��why?��
pilot_num=length(P_f_station);%��Ƶ����

for img=1:N_fft                        %����λ��
    if mod(img,P_f_inter)~=1          %mod(a,b)���������a����b������
        data_station=[data_station,img];
    end
end
data_row=length(data_station);
data_col=ceil(length(spread_data)/data_row);

pilot_seq=ones(pilot_num,data_col)*P_f;%����Ƶ�������
data=zeros(N_fft,data_col);%Ԥ����������
data(P_f_station(1:end),:)=pilot_seq;%��pilot_seq����ȡ

if data_row*data_col>length(spread_data)
    data2=[spread_data;zeros(data_row*data_col-length(spread_data),1)];%�����ݾ����룬��0������Ƶ~
end;data



%% ����ת��
data_seq=reshape(data2,data_row,data_col);
data(data_station(1:end),:)=data_seq;%����Ƶ�����ݺϲ�

%% IFFT
ifft_data=ifft(data)*sqrt(N_fft); 

%% ���뱣�������ѭ��ǰ׺
Tx_cd=[ifft_data(N_fft-N_cp+1:end,:);ifft_data];%��ifft��ĩβN_cp�������䵽��ǰ��

%% ����ת��
Tx_data=reshape(Tx_cd,[],1);%���ڴ�����Ҫ

%% �ŵ���ͨ���ྭ�����ŵ�+AWGN�ŵ���
% %�����ŵ�
% RayleighSinglePath = comm.RayleighChannel(...
%     'SampleRate',sample_rate, ...                  
%     'MaximumDopplerShift',1, ...
%     'DopplerSpectrum',doppler('Jakes'),...
%     'PathGainsOutputPort',true);
%     %'Visualization','Impulse and frequency responses');
%     
%     [n1,h_pathGains]=RayleighSinglePath(Tx_data);%���������ŵ�
% %     Rx_data11=reshape(n1,N_fft+N_cp,[]);
% %     Rx_data21=Rx_data11(N_cp+1:end,:);
% %     fft_data1=fft(Rx_data21)/sqrt(N_fft);
% %     %׼ȷ�ŵ�
% %     H_real=fft_data1./data;
    
% 'PathDelays',[0 30 150 310 370 710 1090 1730 2510]*1e-9, ...
%     'AveragePathGains',[0.0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9], ...

%�ྶ�ŵ�
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
   % Tx_data = randi([0,1],1,30);
    [n1,h_pathGains]=RayleighMultiPath(Tx_data);%�����ྶ�ŵ��Ͷ�����Ƶ��
   %[n1,h_pathGains]=RayleighMultiPath(Tx_data');
%     Rx_data11=reshape(n1,N_fft+N_cp,[]);
%     Rx_data21=Rx_data11(N_cp+1:end,:);
%     fft_data1=fft(Rx_data21)/sqrt(N_fft);
%     %��׼ȷ�ŵ�
%     H_real=fft_data1./data;
 
 
 %׼ȷ�ŵ�(�ྶ)
%PathDelays=[0 4]/sample_rate;
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
 
 %��˹�ŵ�
 S = RandStream('mt19937ar','Seed',5489);
 Ber=zeros(1,length(SNR));
 Ber2=zeros(1,length(SNR));
for jj=1:length(SNR)
    rx_channel=awgn(n1,SNR(jj),'measured',S);%��Ӹ�˹������

%% ����ת��
    Rx_data1=reshape(rx_channel,N_fft+N_cp,[]);
    
%% ȥ�����������ѭ��ǰ׺
    Rx_data2=Rx_data1(N_cp+1:end,:);
%% FFT
    fft_data=fft(Rx_data2)/sqrt(N_fft);
    
%% �ŵ��������ֵ�����⣩

    data3=fft_data(1:N_fft,:); 
    Rx_pilot=data3(P_f_station(1:end),:); %���յ��ĵ�Ƶ
    h=Rx_pilot./pilot_seq; 
    H=interp1( P_f_station(1:end)',h,data_station(1:end)','linear','extrap');%�ֶ����Բ�ֵ����ֵ�㴦����ֵ�����������ڽ������������Ժ���Ԥ�⡣�Գ�����֪�㼯�Ĳ�ֵ����ָ����ֵ�������㺯��ֵ
    
    ERRA=abs((h-H_real(P_f_station(1:end),:))./H_real(P_f_station(1:end),:));
    MSEA(jj)=sum(sum(ERRA.^2))/(length(P_f_station)*data_col);
    
%% �ŵ�У��
    data_aftereq=data3(data_station(1:end),:)./H;
%% ����ת��
    data_aftereq=reshape(data_aftereq,[],1);
    data_aftereq=data_aftereq(1:length(spread_data));
    data_aftereq=reshape(data_aftereq,N_sc,length(data_aftereq)/N_sc);
    
%% ����
    demspread_data = despread(data_aftereq,code);       % ���ݽ���
    
%% QPSK���
    demodulation_data=pskdemod(demspread_data,M,pi/M);    
    De_data1 = reshape(demodulation_data,[],1);
    De_data2 = de2bi(De_data1);
    De_Bit = reshape(De_data2',1,[]);

%% ���⽻֯��
%% �ŵ����루ά�ر����룩
    trellis = poly2trellis(7,[133 171]);
    rx_c_de = vitdec(De_Bit,trellis,tblen,'trunc','hard');   %Ӳ�о�

%% ����������
    [err,Ber2(jj)] = biterr(De_Bit(1:length(code_data)),code_data);%����ǰ��������
    [err, Ber(jj)] = biterr(rx_c_de(1:length(P_data)),P_data);%������������

end
 figure(2);
 semilogy(SNR,Ber2,'b-s');
 hold on;
 semilogy(SNR,Ber,'r-o');
 hold on;
 legend('4PSK���ơ����������ǰ������Ƶ��','4PSK���ơ���������������Ƶ��');
 hold on;
 xlabel('SNR');
 ylabel('BER');
 title('AWGN�ŵ��������������');

 figure(3)
 subplot(2,1,1);
 x=0:1:30;
 stem(x,P_data(1:31));
 ylabel('amplitude');
 title('�������ݣ���ǰ30������Ϊ��)');
 legend('4PSK���ơ�������롢����Ƶ');

 subplot(2,1,2);
 x=0:1:30;
 stem(x,rx_c_de(1:31));
 ylabel('amplitude');
 title('��������(��ǰ30������Ϊ��)');
 legend('4PSK���ơ�������롢����Ƶ');
 
 