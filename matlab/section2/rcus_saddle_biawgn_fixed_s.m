function Pe = rcus_saddle_biawgn_fixed_s(R_vec,n,snr,s)
%
% Pe = rcus_saddle_biawgn_fixed_s(R_vec,n,snr,s):
% Saddlepoint implementation of RCUs for the binary-input AWGN
%
% R_vec = Values of Rate at which the bound will be computed
% n     = Blocklength
% snr   = Value of the snr in linear scale (no dB!)
% s     = Parameter s of the generalized information density (s>0)
%
% OUTPUT:
% Pe = Error probability obtained as a result of the computation of the bound
N   = 500;
tau = linspace(-2,4,N);
psi = zeros(size(tau)); psi1 = psi; psi2 = psi;
Pe = nan(size(R_vec));
maxtau=0;

Is_f = @(z) exp(-z.^2/2)/sqrt(2*pi) .* (log(2) - log(1+exp(-2*s*snr-2*s*sqrt(snr)*z)));
Zmin = -15; Zmax = 15;
I_s = integral(Is_f, Zmin, Zmax);

psi_f = @(t,f) t*(I_s-log(2)) + log(f);
psi1_f = @(f,f1) (I_s-log(2)) + f1/f;
psi2_f = @(f,f1,f2)  (f2*f-f1^2)/f^2;
psi_int = @(z,t) exp(-z.^2/2)/sqrt(2*pi) .* ((1+exp(-2*s*snr-2*s*sqrt(snr)*z)).^t);
psi1_int = @(z,t) exp(-z.^2/2)/sqrt(2*pi) .* (log(1+exp(-2*s*snr-2*s*sqrt(snr)*z)).*(1+exp(-2*s*snr-2*s*sqrt(snr)*z)).^t);
psi2_int = @(z,t) exp(-z.^2/2)/sqrt(2*pi) .* ((log(1+exp(-2*s*snr-2*s*sqrt(snr)*z))).^2.*(1+exp(-2*s*snr-2*s*sqrt(snr)*z)).^t);
    
for tt=1:length(tau)
    tausel = tau(tt);
    psi(tt)    = psi_f(tausel,integral(@(z)psi_int(z,tausel),Zmin,Zmax));
    psi1(tt)    = psi1_f(integral(@(z)psi_int(z,tausel),Zmin,Zmax),integral(@(z)psi1_int(z,tausel),Zmin,Zmax));
    psi2(tt)    = psi2_f(integral(@(z)psi_int(z,tausel),Zmin,Zmax),integral(@(z)psi1_int(z,tausel),Zmin,Zmax),integral(@(z)psi2_int(z,tausel),Zmin,Zmax));
end

Rcr = I_s - psi1_f(integral(@(z)psi_int(z,1),Zmin,Zmax),integral(@(z)psi1_int(z,1),Zmin,Zmax));
% Computation of the RCUs
for ii=1:length(R_vec)
    R = R_vec(ii);
    tauopt = interp1(psi1, tau, I_s-R); % Finding tau s.t. psi'(tau) = threshold.
    if tauopt > maxtau
        maxtau = tauopt;
    end
    if min(psi1)> (I_s-R)
        Pe(ii) = 1;
        continue;
    end
    if tauopt > 1 % R < Rcr
        psiopt = psi_f(1,integral(@(z)psi_int(z,1),Zmin,Zmax));
        psi2opt= psi2_f(integral(@(z)psi_int(z,1),Zmin,Zmax),integral(@(z)psi1_int(z,1),Zmin,Zmax),integral(@(z)psi2_int(z,1),Zmin,Zmax));
        PsiTilde = @(x,y)  exp(n*x*(Rcr-R+psi2opt/2) + log(qfunc(x*sqrt(n*psi2opt) + y*n*(Rcr - R)/sqrt(n*psi2opt))));
        Pe(ii) = exp(n*(psiopt-I_s+R)) *(PsiTilde(1,1) + PsiTilde(0,-1));
    elseif tauopt < 0 %R > Is
        psiopt = psi_f(tauopt,integral(@(z)psi_int(z,tauopt),Zmin,Zmax));
        psi2opt = psi2_f(integral(@(z)psi_int(z,tauopt),Zmin,Zmax),integral(@(z)psi1_int(z,tauopt),Zmin,Zmax),integral(@(z)psi2_int(z,tauopt),Zmin,Zmax));
        Pe(ii) = 1 - exp(log(qfunc(-tauopt.*sqrt(n*psi2opt))) + n*(psiopt-tauopt.*(I_s-R)+0.5*tauopt.^2.*psi2opt))...
                - exp(log(qfunc((1-tauopt).*sqrt(n*psi2opt))) + n*(psiopt-tauopt.*(I_s-R)+0.5*(1-tauopt).^2.*psi2opt));
    else %Rcr < R < Is
        psiopt = psi_f(tauopt,integral(@(z)psi_int(z,tauopt),Zmin,Zmax));
        psi2opt = psi2_f(integral(@(z)psi_int(z,tauopt),Zmin,Zmax),integral(@(z)psi1_int(z,tauopt),Zmin,Zmax),integral(@(z)psi2_int(z,tauopt),Zmin,Zmax));
        Pe(ii) = exp(log(qfunc(tauopt.*sqrt(n*psi2opt))) + n*(psiopt-tauopt.*(I_s-R)+0.5*tauopt.^2.*psi2opt))...
             + exp(log(qfunc((1-tauopt).*sqrt(n*psi2opt))) + n*(psiopt-tauopt.*(I_s-R)+0.5*(1-tauopt).^2.*psi2opt));
    end
    
    if isnan(Pe(ii))
        keyboard;
    end
    if ~isnan(Pe(ii))
        Pe(ii) = min(1, Pe(ii));
    end
end
end