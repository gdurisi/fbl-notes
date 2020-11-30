function Pe = vh_metaconverse_biawgn_fixed_s(R_vec,n,snr,s,N)
%
% Pe = vh_metaconverse_biawgn_fixed_s(R_vec,n,snr,s,N):
% Computation of the Verdú-Han metaconverse bound for the binary-input AWGN channel.
%
% INPUTS:
% R_vec = Values of Rate at which the bound will be computed
% n     = Blocklength
% snr   = Value of the snr in linear scale (no dB!)
% s     = Parameter s of the generalized information density (s>0)
% N     = Number of samples to generate the generalized info. density
%
% OUTPUT:
% Pe = Error probability obtained as a result of the computation of the bound
Pe  = nan(size(R_vec));
i_s=info_dens_biawgn(snr,n,s,N);
mu_f = @(y) (1/(2*sqrt(2*pi)^s) * (exp(-s/2*(y+sqrt(snr)).^2) + exp(-s/2*(y-sqrt(snr)).^2))).^(1/s);
Zmin = -10; Zmax = 10;
logmu = log(integral(mu_f, Zmin, Zmax));
j_s = n*logmu + i_s/s;
[tail,gamma] = ecdf(j_s); % Compute the empirical cdf and save the values where it was evaluated (gamma).
for ii = 1:length(R_vec)
    R       = R_vec(ii);
    Pe(ii)  = max(0,max(tail-exp(gamma-n*R)));
end
end