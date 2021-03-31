% Generate Bspline scaling functino and its dual

function [dual_phi_Bn, phi_Bn] = BSpline(BSplineN)
%BSPLINE Summary of this function goes here
%   Detailed explanation goes here

per = 64; % set period
B0 = ones(1, per); % B0 is the box function
phi_Bn = B0;
for i = 1:BSplineN
    phi_Bn = conv(phi_Bn, B0)/per; % do the convolution
end
%plot(phi_Bn);
% Calculate the dual function
g = phi_Bn' * phi_Bn; % Gram matrix
g = g / norm(g);
dual_phi_Bn = phi_Bn * g; % Get the dual
%figure
%plot(dual_phi_Bn);

end

