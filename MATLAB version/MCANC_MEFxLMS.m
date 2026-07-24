clc;clear;close all;

set(groot,'defaultAxesTickLabelInterpreter','latex');

%% configuration
Fs = 16000; % sampling frequency
T  = 30;     % time
t  = 0:1/Fs:T;
N  = length(t);
Nums  = 6;   % number of secondary source
Nume  = 6;   % number of error sensor

%% system parameters
wLen = 512;  % local control filter length
sLen = 256;  % secondary path length
muw = 1e-6; % control filter step size

%% Primary path
p_low  = 80  ;
p_high = 5000 ;
plen  = 512;
Pri_path = fir1(plen-1,[2*p_low/Fs 2*p_high/Fs]);
P = repmat(Pri_path,Nume,1);
PrimaryPath = reshape(P,[Nume, plen]); % set dimension to Nume*plen

% load measured path
% load("simulation path/PrimaryPath_1x6.mat");
% PrimaryPath = Primary_path;

%% secondary path
s_low = 80 ;
s_high = 5000;
Sec_path = fir1(sLen-1,[2*s_low/Fs 2*s_high/Fs]);
S = repmat(Sec_path,Nume*Nums,1);
SecondaryPath = reshape(S,[Nume, Nums, sLen]); % set dimension to Nume*Nums*slen

% load measured path
% load("simulation path/SecondaryPath_6x6.mat");
% SecondaryPath = Secondary_path;

%% noise generation
noise = randn(N,1);  % random noise
low = 100;
high = 1000;
fil = fir1(63,[2*low/Fs 2*high/Fs]);
Ref = filter(fil,1,noise);   % reference

for i = 1:Nume
  Dis(i,:) = filter(PrimaryPath(i,:),1,Ref);   % Disturbance         
end

Ref = awgn(Ref,40,'measured');

%% MEFxLMS control

MEFxLMS = McANC_SRMSE(wLen,SecondaryPath,sLen,Nums,Nume,Dis,Ref);
[e_CMANC,MEFxLMS] = McFxLMS_SRMSE_166(MEFxLMS,muw);
% [e_CMANC,MEFxLMS] = McFxLMS_SRMSE_ANC(MEFxLMS,muw);

%% draw figure

figure;
for i = 1:Nume

    subplot(3,2,i);
    plot(t(1:T*Fs),Dis(i,1:T*Fs));
    hold on;
    plot(t(1:T*Fs),e_CMANC(i,1:T*Fs))
    if i == 1
        legend('Disturbance','MEFxLMS','Interpreter','latex');
    end
    title(['(', char('a' + i - 1), '). Error ', num2str(i)], 'Interpreter', 'latex');
    xlabel('Time (seconds)', 'Interpreter', 'latex');
    ylabel('Amplitude', 'Interpreter', 'latex');
    grid on;
end



figure;
for i = 1:Nume
    dis = smooth((Dis(i,1:T*Fs).^2),2000);
    ecmanc_adaptive = smooth((e_CMANC(i,1:T*Fs).^2),2000);

    mse_adaptive = 10*log10(ecmanc_adaptive./dis);

    subplot(3,2,i);
    plot(t(1:T*Fs),smooth(mse_adaptive(1:T*Fs),5000));
    if i == 1
        legend('MEFxLMS','Interpreter','latex');
    end
    title(['(', char('a' + i - 1), '). Error ', num2str(i)], 'Interpreter', 'latex');
    xlabel('Time (seconds)', 'Interpreter', 'latex');
    ylabel('Normalized squared error (dB)', 'Interpreter', 'latex');
    axis([0 inf -inf 5]);
    grid on;
end
