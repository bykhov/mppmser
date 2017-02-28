clearvars
close all

%% Definitions
N = 5e6; % number of simulation trials

M = 12; % total number of slots
k = 2; % blocked slots (0's)
n = M - k; % open slotes (1's)
L = nchoosek(M,k); % number of codewords

x = [ones(1,k) zeros(1,n)]'; % Eq. (5) in the paper
% Assumed: mu_0 = 0; mu_1 = 1;
signal = repmat(x,1,N); % sent signal - generated only one time
noise_vect = randn(M,N); % general noise vector - generated only one time

SNRdB = 5:16; % SNR range
SNR = 10.^(SNRdB/10);
% Assumed: sigma_0 = sigma_1 = sigma
sigma = sqrt(1./SNR); % noise std

%% Simulation
for w = 1:length(SNR)        
    % receiver noise    
    noise = sigma(w)*noise_vect;
    
    % signal at the receiver
    signal_r = signal + noise;
    
    % received signal
    % indices of the k lowest values are compared to the transmitted signal
    % (Sec. II-C in the paper)
    [~, r_idx] = sort(signal_r,'descend');
    r_signal_idx = r_idx(1:k,:);
    difference = sort(r_signal_idx) - repmat(find(x),1,N);
    [~, col] = find(difference ~= 0);
    total_errors = length(unique(col));
    ser(w) = total_errors/N;
end

%% Upper Bound (Eq. (17) in paper)
ser_ub = 0;
for p = 1:k
    ser_ub = ser_ub + nchoosek(k,p)*nchoosek(M-k,p)*qfunc(sqrt(p)./sqrt(2)./sigma);
end

%% Theoretical Value (Eq. (15-16 in paper)
%syms v
Q = @(x) 0.5*erfc(x/sqrt(2));
syms v
for w = 1:length(SNR)
    f = @(v) (1 - Q(v/sigma(w)))^k * exp(-(v-1)^2/(2*sigma(w)^2)) * Q((v-1)/sigma(w))^(n-1);
    p_e(w) = 1- n/(sqrt(2*pi)*sigma(w))*double(vpa(int(f(v),v,0,Inf),5));
end

%% Plot Results
h = semilogy(SNRdB,ser,SNRdB,ser_ub,SNRdB,p_e);
h(2).LineStyle = '-.';
h(3).LineStyle = '--';
h(3).Marker = 'o';

legend('Simulation','UB','Theory')
grid on
xlabel('SNR [dB]')
ylabel('SER')
set(gcf, 'Color', 'w');