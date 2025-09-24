%% MCFxLMS fully connected with arbitrary channel [J(Ref) x K(Secondary sources) x M(Error)]

clc;clear;close all;

set(groot,'defaultAxesTickLabelInterpreter','latex');

%% configuration
Fs = 16000; % sampling frequency
T  = 30;     % time
t  = 0:1/Fs:T;
N  = length(t);
J  = 6;   % number of reference sensor
K  = 6;   % number of secondary source
M  = 6;   % number of error sensor

%% Primary path
p_low  = 80  ;
p_high = 5000 ;
Pri_path = fir1(511,[2*p_low/Fs 2*p_high/Fs]);
a=repmat(Pri_path,M*J,1);
PrimaryPath=reshape(a,[M, J, 512]);

% load measured path
% load("path/PrimaryPath_6x6.mat");
% PrimaryPath = Primary_path;


%% Secondary path
s_low = 80 ;
s_high = 5000;
Sec_path = fir1(255,[2*s_low/Fs 2*s_high/Fs]);
b=repmat(Sec_path,M*K,1);
SecondaryPath=reshape(b,[M, K, 256]);

% load measured path
% load("path/SecondaryPath_6x6.mat");
% SecondaryPath = Secondary_path;

%% system parameters
wLen = 512;  % local control filter length
sLen = 256;  % secondary path length
muw = 1e-6; % control filter step size

%% noise generation
noise = randn(N,1);  % random noise
low = 100;
high = 1000;
fil = fir1(63,[2*low/Fs 2*high/Fs]);
X = filter(fil,1,noise);   % noise 

x=repmat(X,J,1);
Ref = reshape(x,[J,length(X)]);  % reference

Dis=zeros(M,N);
for i=1:M
    b=reshape(PrimaryPath(i,:,:),[J,size(PrimaryPath,3)]);
    for ii=1:J
        Dis(i,:)=Dis(i,:)+filter(b(ii,:),1,Ref(ii,:));    %primary response sensed by error transducer
    end
end

%% control
MCFXLMS =  MultiChannelFxLMS(wLen,SecondaryPath,sLen,Ref,Dis,J,K,M);
MCFXLMS =  McFxLMS_controller(MCFXLMS,muw);

%% draw figure

e_CMANC = MCFXLMS.Err;

figure;
for i = 1:M

    subplot(3,2,i);
    plot(t(1:T*Fs),Dis(i,1:T*Fs));
    hold on;
    plot(t(1:T*Fs),e_CMANC(i,1:T*Fs))
    if i == 1
        legend('Disturbance','MCFxLMS','Interpreter','latex');
    end
    title(['(', char('a' + i - 1), '). Error ', num2str(i)], 'Interpreter', 'latex');
    xlabel('Time (seconds)', 'Interpreter', 'latex');
    ylabel('Amplitude', 'Interpreter', 'latex');
    grid on;
end
