function [shifted_kernels, c_matrix] = obtainC_kernels_dB(dBValue)
%OBTAINCMATRIX Summary of this function goes here
%   Detailed explanation goes here
% Use inner product to calcluate cmn


t = 0: 2048-1;
T = 64;
t = t/T;

phi_recon = zeros(1,2048);
[phi_T, psi_T, xval] = wavefun(sprintf('dB%d', dBValue), 6);
phi_recon(1:length(phi_T)) = phi_T;
% Remove zeros
%phi_recon = wshift('1D', phi_recon, T);

phi_sample = phi_recon; % For dB, it is ortho

% Calculate the inner product
len = 32;
% Create shifted kernels
shifted_kernels = [];
for i = 1: len
    shifted_kernels = [shifted_kernels; right_shift(phi_sample, (i - 1) * T )];
end

c_matrix = [];

for i = 1: dBValue
    c_matrix = [c_matrix; t.^(i-1)* shifted_kernels' / T];
end


end

