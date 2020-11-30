function Pe = vh_metaconverse_saddle_biawgn_fixed_s(R_vec,n,snr,s)
%
% Pe = vh_metaconverse_saddle_biawgn_fixed_s(R_vec,n,snr,s):
% Saddlepoint implmentation of the Verdú-Han metaconverse bound for the
% binary-input AWGN channel.
%
% INPUTS:
% R_vec = Values of Rate at which the bound will be computed
% n     = Blocklength
% snr   = Value of the snr in linear scale (no dB!)
% s     = Parameter s of the generalized information density (s>0)
%
% OUTPUT:
% Pe = Error probability obtained as a result of the computation of the bound
N   = 500;
tau = linspace(0,4,N);
Pe  = nan(size(R_vec));
psi = zeros(size(tau)); psi1 = psi; psi2 = psi;

Is_f = @(z) exp(-z.^2/2)/sqrt(2*pi) .* (log(2) - log(1+exp(-2*s*snr-2*s*sqrt(snr)*z)));
Zmin = -15; Zmax = 15;
I_s = integral(Is_f, Zmin, Zmax);
mu_f = @(y) (1/(2*sqrt(2*pi)^s) * (exp(-s/2*(y+sqrt(snr)).^2) + exp(-s/2*(y-sqrt(snr)).^2))).^(1/s);
Zmin = -20; Zmax = 20;
logmu = log(integral(mu_f, Zmin, Zmax));
J_s = logmu + I_s/s;

psi_f = @(t,f) t*(I_s-log(2)) + log(f);
psi1_f = @(t,f,f1) (I_s-log(2)) + f1/f;
psi2_f = @(t,f,f1,f2)  (f2*f-f1^2)/f^2;
psi_int = @(z,t) exp(-z.^2/2)/sqrt(2*pi) .* ((1+exp(-2*s*snr-2*s*sqrt(snr)*z)).^t);
psi1_int = @(z,t) exp(-z.^2/2)/sqrt(2*pi) .* (log(1+exp(-2*s*snr-2*s*sqrt(snr)*z)).*(1+exp(-2*s*snr-2*s*sqrt(snr)*z)).^t);
psi2_int = @(z,t) exp(-z.^2/2)/sqrt(2*pi) .* ((log(1+exp(-2*s*snr-2*s*sqrt(snr)*z))).^2.*(1+exp(-2*s*snr-2*s*sqrt(snr)*z)).^t);

for tt=1:length(tau)
    tausel = tau(tt);
    psi(tt)    = psi_f(tausel,integral(@(z)psi_int(z,tausel),Zmin,Zmax));
    psi1(tt)    = psi1_f(tausel,integral(@(z)psi_int(z,tausel),Zmin,Zmax),integral(@(z)psi1_int(z,tausel),Zmin,Zmax));
    psi2(tt)    = psi2_f(tausel,integral(@(z)psi_int(z,tausel),Zmin,Zmax),integral(@(z)psi1_int(z,tausel),Zmin,Zmax),integral(@(z)psi2_int(z,tausel),Zmin,Zmax));
end

for ii=1:length(R_vec)
    R = R_vec(ii);
    tail_est = exp(log(qfunc(tau.*sqrt(n*psi2))) + n*(psi-tau.*psi1+0.5*tau.^2.*psi2));
    Pe(ii) = max(0,max(tail_est - exp(n*(J_s-psi1/s-R))));
end

end