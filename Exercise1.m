clear all
%% Use inner product to calcluate cmn

t = 0: 2048 - 1; % Set time
T = 64; % Set period
t = t/T; % Set the real time

% t constant
t0 = ones(1, 2048); % Get constant

% t linear
t1 = t;

% t square:
t2 = t.^2;

% t^3
t3 = t.^3;

tn = [t0; t1; t2; t3]; % Get a matrix

phi_recon = zeros(1,2048); 
[phi_T, psi_T, xval] = wavefun('dB4', 6);
phi_recon(1:length(phi_T)) = phi_T; % Get the reconstruction scaling function
% Remove zeros
%phi_recon = wshift('1D', phi_recon, T);

phi_sample = phi_recon; % Get the sampling signal

% Calculate the inner product
len = 32;
% Create shifted kernels
shifted_kernels = [];
for i = 1: len
    shifted_kernels = [shifted_kernels; right_shift(phi_sample, (i - 1) * T )];
end
% Sampling (inner product)
for i = 1: len
    c0n = t0* shifted_kernels' / T;
    c1n = t1* shifted_kernels' / T;
    c2n = t2* shifted_kernels' / T;
    c3n = t3* shifted_kernels' / T;
end
cn = [c0n; c1n; c2n; c3n]; % Combine the sampled signals

%% Reconstruction
fx0 = c0n * shifted_kernels;
fx1 = c1n * shifted_kernels;
fx2 = c2n * shifted_kernels;
fx3 = c3n * shifted_kernels;
fx = [fx0; fx1; fx2; fx3];

figure;
%% Plotting all the diagrams
for i = 1: 4
    subplot(2,2,i);
    h1 = plot(t, tn(i,:), 'k');
    hold on
    h2 = plot(t, fx(i,:), 'r');    
    hold on
    h3 = plot(t, (cn(i,:) .* shifted_kernels')', 'g');    
    legend([h1, h2, h3(1)], {'original signal', 'reproduced', 'kernels'});
    xlabel('t');
    ylabel('x(t)');
    title_str = ['Reproduced polynomials (using dB4) degree ' num2str(i-1)];
    title(title_str);
end