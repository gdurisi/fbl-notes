function id=info_dens_biawgn(snr,n,s,N)
%
% id = info_dens_biawgn(snr,n,s,N):
% Generation of samples of the generalized information density.
%
% INPUTS:
% snr   = Value of the snr in linear scale (no dB!)
% n     = blocklength 
% s     = parameter s of the generalized information density (s>0)
% N     = number of samples to generate
%
% OUTPUT:
% id    = samples of the generalized information density 

normrv=randn(n,N)*sqrt(snr)+snr;
idvec=log(2)-log(1+exp(-2*s*normrv));
id=sum(idvec,1);
end