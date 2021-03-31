clear all
load('tau.mat') % Load the tau data

max_len = 2048;
t = 0: max_len - 1; % Set time
T = 64; % Set period
t = t/T; % Get the real time

% Get the filter h
K = 2;
degree_L = 3; %degree >= 2K-1
N = degree_L+1; 

%prog_K = K + 1;
% Create matrix tau
matrix_tau=[];
for i = 1:K
    one_line = [];
    for j = 1:K
        one_line = [one_line, tau(K + i - j)];
    end
    matrix_tau = [matrix_tau; one_line];
end

y_tau_1 = - tau(K + 1: N);
h = matrix_tau \ y_tau_1'; % Solve the matrix to get filter h

% Get the t_n

syms x;
decom = 1; % h0 is 1
for i = 1:K
    decom = decom + (x ^ (-i)) * h(i);
end

t_k=double(solve(decom,x)); %By solving the equation, we can get t_k.
%t_k = sort(1./result);
%decom_new = factor(decom)

% Calculate a_n
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
a_k = matrix_t \ y_tau_2'; % Solve the matrix to get ak
% Plot all the diagrams
figure;
stem(t_k,a_k);
xlabel('t');
ylabel('amplitude');
grid('on');
for i = 1:K
    fprintf('location %d is %f, and the amplitude is %f\n', i,t_k(i),a_k(i));
end
