clear all
T = 64;
% Set Diracs
amp = [5.755,10.1];
loca = [8.8125,14.9531];
for momen = 5: 10 %Try different degrees, from 5 to 10
% Get the shifted kernels
[shifted_kernels, c_matrix] = obtainC_kernels_dB(momen);

% Generate diracs

signal_diracs = zeros(1, 2048);
for i = 1 : length(amp)
    signal_diracs(ceil(loca(i) * T)) = amp(i);
end

% Generate s and tau
yn =  signal_diracs * shifted_kernels'; % Sampled signal
tau = yn * c_matrix';

% Generate noise
noise_var = 1; % noise variance
noise_mean = 0;
noise = noise_mean + rand(1, momen) * sqrt(noise_var);

% Add noise
tau = tau + noise; % Need to be the same signal

% Get the filter h using TLS
N = momen;
K = 2 + 1;
%prog_K = K + 1;
% Create matrix tau
matrix_tau=[];
for i = 1:N-K+1
    one_line = [];
    for j = 1:K
        one_line = [one_line, tau(K + i - j)];
    end
    matrix_tau = [matrix_tau; one_line];
end
% Conduct the SVD algorithm and solve the matrix to get h
[tau_U, tau_S, tau_V] = svd(matrix_tau);
h = tau_V(:, length(tau_V));
%clear tau_U
%clear tau_S
%clear tau_V
% Get the t_n
% Create equation
syms x;
%decom = 1; % h0 is 1
decom=h(1);
for i = 1:K-1
    decom = decom + (x ^ (-i)) * h(i+1);
end

result=double(solve(decom,x)); % Solve equation to tget tn
%clear x;
t_TLS = sort(result);
%decom_new = factor(decom)

% Calculate a_n
matrix_t=[];
for i = 1:K-1
    matrix_t = [matrix_t, 1]; %The first line
end
    
for i = 1:K-1-1
    one_line = [];
    for j = 1:K-1
        one_line = [one_line, t_TLS(j) ^ i];
    end
    matrix_t = [matrix_t; one_line];
end

y_tau_2 = tau(1:K-1);
a_TLS = matrix_t \ y_tau_2'; % Solve the matrix to get ak

% Get the filter h using Cadzow
%prog_K = K + 1;
% Create matrix tau
matrix_tau=[];
for i = 1:N-K+1
    one_line = [];
    for j = 1:K
        one_line = [one_line, tau(K + i - j)];
    end
    matrix_tau = [matrix_tau; one_line];
end


%h = tau_V(:, length(tau_V));

while 1
    [tau_U, tau_S, tau_V] = svd(matrix_tau); % Do SVD algorithm
    if rank(matrix_tau) == K-1
        break;
    end
    for i = 3: size(tau_S, 2)
        tau_S(i, i) = 0; %Keep the K largest diagonal coefficients and set the others to zero.
    end
    matrix_tau = tau_U * tau_S * tau_V';
    toep = zeros(size(matrix_tau));
    %Make it as toeplitz
    for r = 1: size(matrix_tau, 1)
        for c = 1: size(matrix_tau, 2)
            toep(r, c) = mean(diag(matrix_tau, c - r));
        end
    end
    matrix_tau = toep;
    %[tau_U, tau_S, tau_V] = svd(matrix_tau);
    %disp('hello');
    
end

%[tau_U, tau_S, tau_V] = svd(matrix_tau);
h = tau_V(:, length(tau_V));

% Get the t_n
% Create equation
syms x;
%decom = 1; % h0 is 1
decom=h(1);
for i = 1:K-1
    decom = decom + (x ^ (-i)) * h(i+1);
end

result=double(solve(decom,x)); % Solve the equation
%clear x;
t_Cadzow = sort(result);
%decom_new = factor(decom)

% Calculate a_n
matrix_t=[];
for i = 1:K-1
    matrix_t = [matrix_t, 1]; %The first line
end
    
for i = 1:K-1-1
    one_line = [];
    for j = 1:K-1
        one_line = [one_line, t_Cadzow(j) ^ i];
    end
    matrix_t = [matrix_t; one_line];
end

y_tau_2 = tau(1:K-1);
a_Cadzow = matrix_t \ y_tau_2'; % Solve the matrix to get a

% Plotting all
subplot(3, 2, momen - 4);
h1 = stem(t_TLS, a_TLS, 'filled'); %The results by using TLS
hold on;
h2 = stem(t_Cadzow, a_Cadzow); %The results by using Cadzow
hold on;
h3 = stem(loca, amp);
hold on;
legend([h1, h2, h3], {'TLS','Cadzow', 'Original Diracs'});
grid('on');
xlabel('t');
ylabel('amplitude');
fprintf('TLS reconstruction, dB%d --------\n', momen);
for i = 1: K - 1
    fprintf('location %d is %f, and the amplitude is %f\n', i,t_TLS(i),a_TLS(i));
end
fprintf('Cadzow reconstruction, dB%d --------\n', momen);
for i = 1: K - 1
    fprintf('location %d is %f, and the amplitude is %f\n', i,t_Cadzow(i),a_Cadzow(i));
end
title(sprintf('The reconstruction results when N=%d', momen));
end
