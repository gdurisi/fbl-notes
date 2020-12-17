function Pe = rcus_biawgn_fixed_s(R_vec,n,snr,s,N)
%
% Pe = rcus_biawgn_fixed_s(R_vec,n,snr,s,N):
% Computation of the RCUs bound for the binary-input AWGN channel.
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
i_s = info_dens_biawgn(snr,n,s,N);
for ii = 1:length(R_vec)
    R       = R_vec(ii);
    Pe(ii)  = mean( exp( - max(0, i_s - log(exp(n*R)-1))));
end
end