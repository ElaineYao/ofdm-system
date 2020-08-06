clc;
clear;

%% Parameters
s=rng;
N_sc=128;      % Number of subcarrierA
N_fft=128;            % Length of FFT 
N_cp=8;             % Cyclic prefix
N_symbo=N_fft+N_cp;        % An OFDM symbol
%N_c=53;             % Number of carriers(with DC carriers included)
M=4;               % Modulation - 4PSK 
SNR=0:5:30;         % SNR 
N_frm=1;           % Number of frames
%N_frm=53;            % Number of frames at per SNR
%Nd=19040;               % Number of bits in each frame
Nd=N_symbo*20;               % Number of bits in each frame
sample_rate = 1.92e6; % Sample rate
P_t_inter=8;      % Interval between pilot frequencies in time domain
P_f_inter=4;      % Interval between pilot frequencies in time domain
data_station=[];    % Position of pilot frequencies
L=7;                % Constraint length of Convolutional Code
tblen=6*L;          % Depth of Viterbi decoder
stage = 3;          % Order of m-sequence
ptap1 = [1 3];      % Shift-register connections for m-sequences
regi1 = [1 1 1];    % Initial value for registers in m-sequnces
N = Nd*N_frm;       % Number of bits

%% Generate baseband data
P_data=randi([0 1],1,N_sc*Nd*N_frm);


%% Channel encoding (Convolutional Code & Interweaving)
trellis = poly2trellis(7,[133 171]);       % (2,1,7)Convolutional code
code_data=convenc(P_data,trellis);


%% Modulation - QPSK
data_temp1= reshape(code_data,log2(M),[])';             % Group the data by 2，M=4
data_temp2= bi2de(data_temp1);                             % Transform the binary data into decimal
modu_data=pskmod(data_temp2,M,pi/M);              % QPSK
% figure(1);
scatterplot(modu_data),grid;                  % Constellation Diagram

%% Spread the spectrum
code = mseq(stage,ptap1,regi1,N_sc);     % Generate the spreading code
code = code * 2 - 1;         % Transform 1,0 to 1,-1
modu_data=reshape(modu_data,N_sc,length(modu_data)/N_sc);
spread_data = spread(modu_data,code);        % Spread the spectrum
spread_data=reshape(spread_data,[],1);

%% Insert the pilot frequency (Pattern C)

P_f=3+3*1i;                       % Pilot frequency
P_f_station=[1:P_f_inter:N_fft];
Length_P_f=length(P_f_station);
Length_P_t=floor(length(spread_data)/(N_fft*P_t_inter-Length_P_f)); % Number of pilot frequencies in time domain
remain_P=mod(length(spread_data),(N_fft*P_t_inter-Length_P_f)); % Mod remainder
Length_P_t1=ceil(length(spread_data)/(N_fft*P_t_inter-Length_P_f)); % Maximum length of the matrix in time domain
P_t_station=[1:P_t_inter:(Length_P_t-1)*P_t_inter+1]; % 
% Construct data matrix
if (remain_P==0)
    data=zeros(N_fft,Length_P_t*P_t_inter);
else if remain_P<N_fft-Length_P_f+1
    data=zeros(N_fft,Length_P_t*P_t_inter+1);
    else
        r=ceil((remain_P-N_fft+Length_P_f)/N_fft);
        data=zeros(N_fft,Length_P_t*P_t_inter+1+r);
    end
end

% Insert the pilot frequency into the matrix
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

%% Serial-to-parallel conversion
data=reshape(data,data_row,data_col);

%% IFFT
ifft_data=ifft(data)*sqrt(N_fft); 

%% Insert the guard interval and cyclic prefix
Tx_cd=[ifft_data(N_fft-N_cp+1:end,:);ifft_data];

%% Parallel-to-serial conversion
Tx_data=reshape(Tx_cd,[],1);

%% Channel (Rayleigh Multipath Channel & AWGN Channel)
    
% 'PathDelays',[0 30 150 310 370 710 1090 1730 2510]*1e-9, ...
%     'AveragePathGains',[0.0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9], ...

%  'PathDelays',[0 1 2 3 4 5]/sample_rate, ...
%     'AveragePathGains',[0.0 -1.5  -3.6  -7.0 -12.0 -16.9], ...

% Rayleigh Multipath Channel
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
    [n1,h_pathGains]=RayleighMultiPath(Tx_data);

 
 % Calculate the channel parameters
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
 
 % AWGN channel
 S = RandStream('mt19937ar','Seed',5489);
 Ber=zeros(1,length(SNR));
 Ber2=zeros(1,length(SNR));
 for jj=1:length(SNR)
    rx_channel=awgn(n1,SNR(jj),'measured',S); % Add Gaussian noise

%% Serial-to-parallel conversion
    % Trim the channel delay, and add 0 
    rx_channel_dp = rx_channel(8:end);
    rx_channel_b0 = [rx_channel_dp;zeros(7,1)];
    Rx_data1=reshape(rx_channel_b0,N_fft+N_cp,[]);
    
%% Drop the guard interval and cyclic prefix  
    Rx_data2=Rx_data1(N_cp+1:end,:);
%% FFT
    fft_data=fft(Rx_data2)/sqrt(N_fft);
    
%% Channel estimation & Interpolation

     data3=fft_data(1:N_fft,:); 
     Rx_pilot=data3(P_f_station(1:end),P_t_station(1:end)); % Pilot frequency received
     pilot_seq = data(P_f_station(1:end),P_t_station(1:end)); % Pilot frequency sent 
     h=Rx_pilot./pilot_seq; 
     
     loc_f = 1:N_fft;
     loc_t = 1:data_col;
     % Frequency domain
     H_f=interp1( P_f_station(1:end),h,loc_f,'spline');
     % Time domain
     H_t=interp1( P_t_station(1:end)',H_f',loc_t','spline');
     H = H_t';
%     ERRA=abs((h-H_real(P_f_station(1:end),:))./H_real(P_f_station(1:end),:));
%     MSEA(jj)=sum(sum(ERRA.^2))/(length(P_f_station)*data_col);
    
%% Correct the estimated channel
      data_s = reshape(data3,[],1);
      H_s = reshape(H,[],1);
      data_aftereq=data_s(ind(1:end))./H_s(ind(1:end));
%% Parallel-to-serial conversion
     data_aftereq=data_aftereq(1:length(spread_data));
     data_aftereq=reshape(data_aftereq,N_sc,length(data_aftereq)/N_sc);
    
%% De-spread the frequency
    demspread_data = despread(data_aftereq,code);      
    
%% Demodulation - QPSK
    demodulation_data=pskdemod(demspread_data,M,pi/M);    
    De_data1 = reshape(demodulation_data,[],1);
    De_data2 = de2bi(De_data1);
    De_Bit = reshape(De_data2',1,[]);

%% Decode the channel (Viterbi)
    trellis = poly2trellis(7,[133 171]);
    rx_c_de = vitdec(De_Bit,trellis,tblen,'trunc','hard');   

%% Calculate the bit error rate
    [err,Ber2(jj)] = biterr(De_Bit(1:length(code_data)),code_data); % Before decoding
    [err, Ber(jj)] = biterr(rx_c_de(1:length(P_data)),P_data); %After decoding

end
 figure(2);
 semilogy(SNR,Ber2,'b-s');
 hold on;
 semilogy(SNR,Ber,'r-o');
 hold on;
 legend('4PSK Modulation. Before decoding the Convolutional Code (with spectrum spreaded)','4PSK Modulation. After decoding the Convolutional Code (with spectrum spreaded)');
 hold on;
 xlabel('SNR');
 ylabel('BER');
 title('Bit error rate under Rayleigh Multipath Channel');

 figure(3)
 subplot(2,1,1);
 x=0:1:30;
 stem(x,P_data(1:31));
 ylabel('amplitude');
 title('Data sent');
 legend('4PSK Modulation','With Convolutional Code','With spectrum spreaded');

 subplot(2,1,2);
 x=0:1:30;
 stem(x,rx_c_de(1:31));
 ylabel('amplitude');
 title('Data received');
 legend('4PSK Modulation','With Convolutional Code','With spectrum spreaded');
 



