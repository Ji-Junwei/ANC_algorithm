%% single channel ANC

clc;clear;close all;

set(groot,'defaultAxesTickLabelInterpreter','latex');


Fs = 16000; % sampling frequency
T  = 30;     % time (s)
t  = 0:1/Fs:T;
N  = length(t);

%% Primary path
p_low  = 20  ;
p_high = 6000 ;
Pri_path = fir1(511,[2*p_low/Fs 2*p_high/Fs]);

%% secondary path
s_low = 20 ;
s_high = 6000;
Sec_path = fir1(255,[2*s_low/Fs 2*s_high/Fs]);

%% noise generation
noise = randn(N,1);  % random noise
low = 100;
high = 1000;
fil = fir1(127,[2*low/Fs 2*high/Fs]); % broadband 100 to 1,000Hz
Ref = filter(fil,1,noise);   % reference
Dis = filter(Pri_path,1,Ref);

%% system parameters
wlen = 512;  % control filter length
slen = 256;  % secondary path length
muw = 0.001; % step size

%% call matlab function https://www.mathworks.com/help/releases/R2025a/dsp/ref/dsp.filteredxlmsfilter-system-object.html?searchPort=54096#bva0bw8
fxlms = dsp.FilteredXLMSFilter(wlen, 'StepSize', muw, 'SecondaryPathCoefficients', Sec_path,'SecondaryPathEstimate',Sec_path);
[y,e] = fxlms(Ref,Dis); % y: control signal; e: error signal; e = Dis + filter(Sec_path,1,y)

figure;
plot(t(1:T*Fs),Dis(1:T*Fs));
hold on;
plot(t(1:T*Fs),e(1:T*Fs))
title('Active Noise Control','Interpreter','latex');
legend('Disturbance','Residual error signal','Interpreter','latex');
xlabel('Time Index'); 
ylabel('Signal Value');

%% FxLMS algorithm
single_ANC = ANC_algorithm(wlen,slen,Sec_path',Dis,Ref);
[err,single_ANC] = ANC_FxLMS(single_ANC,muw);

figure;
plot(t(1:T*Fs),Dis(1:T*Fs));
hold on;
plot(t(1:T*Fs),err(1:T*Fs))
title('Active Noise Control','Interpreter','latex');
legend('Disturbance','Residual error signal','Interpreter','latex');
xlabel('Time Index'); 
ylabel('Signal Value');

% compare control filter
figure;
plot(-fxlms.Coefficients);
hold on;
plot(single_ANC.Wc);
title('Coefficients comparsion','Interpreter','latex');
legend('Matlab function','FxLMS','Interpreter','latex');
xlabel('Tap length'); 
ylabel('Coefficients Value');
