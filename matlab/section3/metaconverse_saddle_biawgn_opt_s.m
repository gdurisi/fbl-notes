function epsilon = metaconverse_saddle_biawgn_opt_s(snr,R_vec,n)
%
% epsilon = metaconverse_saddle_biawgn_opt_s(snr,R_vec,n):
% Saddlepoint implementation of metaconverse for the binary-input AWGN, based on:
% G. Vazquez-Vilar et al. "Saddlepoint Approximation of the Error Probability
% of Binary Hypothesis Testing",  ISIT 2018
%
% INPUTS:
% snr       = Value of the snr in linear scale (no dB!)
% R_vec     = Values of Rate at which the bound will be computed
% n         = Blocklength
%
% OUTPUT:
% epsilon    = Error probability obtained as a result of the computation of the bound
%
% IMPORTANT NOTE: The accuracy of the saddlepoint approximation depends on
% the right choice of the parameter rho. If non-sense values are obtained
% after evaluation, consider to increase the upper limit of the variable
% rhorange. 

Zmax = 10; % determines support for numerical integration routines
rhorange = -0.9:0.01:6; % range of values to optimize the bound over; note that this might need to be adjusted depending on the target snr, R, and n.

epsilon=-Inf*ones(size(R_vec));
rho_optimal=-1*ones(size(R_vec));

for ii=1:length(R_vec)
    R = R_vec(ii);
    for jj=1:1:length(rhorange)
        rho=rhorange(jj);
        [E0,~, K1, K2]=Kplus(snr,rho,1/(1+rho),Zmax);
        
        V=K2/(1+rho)^2;
        
        etan=0.5*(erfcx(sqrt(n*V/2))+sign(rho)*erfcx(sqrt(n*rho^2*V/2)));
        
        E1=(E0+K1)/(1+rho);
        epsilon_current= double(rho<0)+(etan-exp(-n*(R-E1)))*exp(-n*(E0-rho*E1));
        
        if (epsilon_current>epsilon(ii))
            epsilon(ii)=epsilon_current;
            rho_optimal(ii)=rho;
        end
    end
end

end

function [E0,K0, K1, K2]=Kplus(snr,rho,s,Zmax)
%
% The function evaluated K_P(s) and its first and second derivatives for a tilted Q distribution tilted by the parameter rho.
% Here, s is the parameter of the moment generating function, whereas rho is the tilting parameter of the output distribution.

% first compute log \mu0
g0 = @(z) exp(-z.^2/2)/sqrt(2*pi) .* exp(-rho*(log(2) - log(1+exp(-(2/(1+rho))*(z*sqrt(snr)+snr)))));

mu0=integral(g0, -Zmax, Zmax);

logmu0=log(mu0);

E0=-logmu0;

% now evaluate K(s) and its derivatives.
f0 = @(z) exp(-z.^2/2)/sqrt(2*pi) .* exp(-(1-s)*(logmu0+ (1+rho)*log(2) - (1+rho)*log(1+exp(-(2/(1+rho))*(z*sqrt(snr)+snr)))));

f1= @(z) exp(-z.^2/2)/sqrt(2*pi) .*(logmu0+ (1+rho)*log(2) - (1+rho)*log(1+exp(-(2/(1+rho))*(z*sqrt(snr)+snr)))).* exp(-(1-s)*(logmu0+ (1+rho)*log(2) - (1+rho)*log(1+exp(-(2/(1+rho))*(z*sqrt(snr)+snr)))));

f2= @(z) exp(-z.^2/2)/sqrt(2*pi) .*((logmu0+ (1+rho)*log(2) - (1+rho)*log(1+exp(-(2/(1+rho))*(z*sqrt(snr)+snr)))).^2).* exp(-(1-s)*(logmu0+ (1+rho)*log(2) - (1+rho)*log(1+exp(-(2/(1+rho))*(z*sqrt(snr)+snr)))));

M0=integral(f0, -Zmax, Zmax);
M1=integral(f1, -Zmax, Zmax);
M2=integral(f2, -Zmax, Zmax);
K0=log(M0);
K1=M1/M0;
K2= (M2*M0-M1^2)/M0^2;
%sqK2/log(2)^2;
end