clear all

momen = 4; % Set moments
[shifted_kernels, c_matrix] = obtainC_kernels_dB(momen); % Create shifted kernels
T = 64;
max_len = 2048;
t = 0: max_len - 1;
T = 64;
t = t/T;

%% Load the given data
load('samples.mat');
tau = y_sampled * c_matrix'; % Create tau matrix

%% Get the filter h
N = momen;
K = 2;
%prog_K = K + 1;
% Create matrix tau
matrix_tau=[];
for i = 1:N-K
    one_line = [];
    for j = 1:K
        one_line = [one_line, tau(K + i - j)];
    end
    matrix_tau = [matrix_tau; one_line];
end

y_tau_1 = - tau(K + 1: N);
h = matrix_tau \ y_tau_1'; % Solve the matrix to get filter h

%% Get the t_n
% Create the equation
syms x;
decom = 1; % h0 is 1
for i = 1:K
    decom = decom + (x ^ (-i)) * h(i);
end

t_k=double(solve(decom,x)); %By solving the equation, we can get t_k.

%% Calculate a_n
matrix_t=[];
for i = 1:K
    matrix_t = [matrix_t, 1]; %The first line
end
    
for i = 1:K-1
    one_line = [];
    for j = 1:K
        one_line = [one_line, t_k(j) ^ i];
    end
    matrix_t = [matrix_t; one_line];
end

y_tau_2 = tau(1:K);
a_k = matrix_t \ y_tau_2'; % Solve the matrix to get the ak

%% Plotting recovered Diracs
figure;
stem(t_k,a_k);
xlabel('t');
ylabel('amplitude');
grid('on');
for i = 1:K
    fprintf('location %d is %f, and the amplitude is %f\n', i,t_k(i),a_k(i));
end