% The function for right shifting

function [shifted] = right_shift(signal, N)
    len=length(signal);
    shifted = zeros(1, len);
    shifted(N+1: len) = signal(1: len-N); % Shift the signal
end
