function R_bits = metaconverse_biawgn_fixed_s(eps_vec,n,snr,s,N)
%
% R_bits = metaconverse_biawgn_fixed_s(eps_vec,n,snr,s,N):
% Computation of the metaconverse bound for the binary-input AWGN channel.
%
% INPUTS:
% eps_vec = Values of the error probability at which the bound will be
%             evaluated
% n       = Blocklength
% snr     = Value of the snr in linear scale (no dB!)
% s       = Parameter s of the generalized information density (s>0)
% N       = Number of samples to generate the generalized info. density
%
% OUTPUT:
% R_bits = Rate (bits) obtained as a result of the computation of the bound

R_bits  = nan(size(eps_vec));
i_s=info_dens_biawgn(snr,n,s,N);
mu_f = @(y) (1/(2*sqrt(2*pi)^s) * (exp(-s/2*(y+sqrt(snr)).^2) + exp(-s/2*(y-sqrt(snr)).^2))).^(1/s);
Zmin = -10; Zmax = 10;
logmu = log(integral(mu_f, Zmin, Zmax));
j_s = n*logmu + i_s/s;
[cdf,gamma] = ecdf(j_s); % Compute the empirical cdf and save the values where it was evaluated (gamma).
for ii = 1:length(eps_vec)
    eps = eps_vec(ii);
    mu = interp1(cdf ,gamma, eps);
    R_bits(ii) = -1/n*log2(mean(exp(-j_s).*(j_s>= mu)));
end
end