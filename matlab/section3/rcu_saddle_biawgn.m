function [epsilon,rho_current] = rcu_saddle_biawgn(snr,R_vec,n)
%
% [epsilon,rho_current] = rcu_saddle_biawgn(snr,R_vec,n):
% Saddlepoint implementation of RCU for the binary-input AWGN, based on:
% J. Font-Segura et al. "Saddlepoint Approximations of Lower and Upper
% Bounds to the Error Probability in Channel Coding", CISS 2018
%
% INPUTS:
% snr   = Value of the snr in linear scale (no dB!)
% R_vec = Values of Rate at which the bound will be computed
% n     = Blocklength
%
% OUTPUT:
% epsilon     = Error probability obtained as a result of the computation of the bound
% rho_current = Value of the parameter rho for which stationary equation is
%               satisfied
 
Zmax=10; % defines the range for numerical integration [-Zmax, Zmax]; typically, Zmax=10 sufficies
rho_current=zeros(size(R_vec));
epsilon = ones(size(R_vec));
Rcritical=E0prime(snr,1,1/2,Zmax);
Capacity=E0prime(snr,0,1,Zmax);

rho_vec = linspace(-0.9,4,1000);
s_vec=1./(1+rho_vec);
E0prime_vec = nan(size(rho_vec));
for ii=1:length(rho_vec)
    E0prime_vec(ii) = E0prime(snr,rho_vec(ii),s_vec(ii),Zmax);
end
rem_idx = find(isnan(E0prime_vec));
E0prime_vec(rem_idx) = [];
rho_vec(rem_idx) = [];
for ii = 1:length(R_vec)
    R = R_vec(ii);
    rho_current(ii) = interp1(E0prime_vec,rho_vec,R);
    if (R>Capacity)
        [~,epsilon(ii)] = rcu_fixed(snr,n,Zmax,rho_current(ii));
        epsilon(ii) = 1 + epsilon(ii);
    elseif (R<Rcritical)
        %disp(['The rate R=' num2str(R) ' is below the critical rate. I will set epsilon=1!'])
        %[rho_current(ii),~]=search_rho(snr,Zmax,R); % Obtain the rho satisfying the stationary equation
        
        if isnan(rho_current(ii))
            keyboard
        end
        [~,epsilon(ii)]= rcu_fixed(snr,n,Zmax,rho_current(ii));
        E01=E0plus(snr,1,1/(1+1),Zmax);
        W21=computeW2(snr,1,Zmax,E01);
        thetan1=(1/sqrt(1+1))*((1+1)/(sqrt(2*pi*n*W21)) )^1;
        xi1 = exp(n*(R-E01))*thetan1;
        epsilon(ii) = epsilon(ii) + xi1; 
    else
        %[rho_current_old(ii),~]=search_rho(snr,Zmax,R); % Obtain the rho satisfying the stationary equation
       
        [~,epsilon(ii)]=rcu_fixed(snr,n,Zmax,rho_current(ii));
    end
end
end

function [rate,epsilon]=rcu_fixed(snr,n,Zmax,rho)
  
[Eres,E1res,E2res]=E0plus(snr,rho,1/(1+rho),Zmax);
V=-E2res;
rate=E1res;

W2=computeW2(snr,rho,Zmax,Eres);

thetan=(1/sqrt(1+rho))*((1+rho)/(sqrt(2*pi*n*W2)) )^rho;

epsilon= exp(n*(rho*E1res-Eres))*0.5*( erfcx(abs(rho*sqrt(n*V/2)))*sign(rho) + erfcx(abs((1-rho)*sqrt(n*V/2)))*sign(1-rho))*thetan;

%eps_gallager=exp(n*(rho*E1res-Eres));

end

function W2=computeW2(snr,rho,Zmax,E0)
  
  logmu0=-E0;
  
  f0 = @(z) exp(-z.^2/2)/sqrt(2*pi) .* exp(-(logmu0+ (1+rho)*log(2) - (1+rho)*log(1+exp(-(2/(1+rho))*(z*sqrt(snr)+snr))))).*(2*(z*sqrt(snr)+snr)).^2.*exp(-(2/(1+rho))*(z*sqrt(snr)+snr))./(1+exp(-(2/(1+rho))*(z*sqrt(snr)+snr))).^2;
  
  W2=integral(f0, -Zmax, Zmax);
  
end

function res=E0prime(snr,rho,s,Zmax) 
  
 f0 = @(z) exp(-z.^2/2)/sqrt(2*pi) .* exp(-rho*(log(2) - log(1+exp(-2*s*(z*sqrt(snr)+snr)))));
 
 f1= @(z) exp(-z.^2/2)/sqrt(2*pi) .*(log(2) - log(1+exp(-2*s*(z*sqrt(snr)+snr)))).* exp(-rho*(log(2) - log(1+exp(-2*s*(z*sqrt(snr)+snr)))));
 
 res= integral(f1, -Zmax, Zmax)/integral(f0,-Zmax,Zmax);
 
end

function [E0, E1, E2]=E0plus(snr,rho,s,Zmax)
  
  f0 = @(z) exp(-z.^2/2)/sqrt(2*pi) .* exp(-rho*(log(2) - log(1+exp(-2*s*(z*sqrt(snr)+snr)))));
 
  f1= @(z) exp(-z.^2/2)/sqrt(2*pi) .*(log(2) - log(1+exp(-2*s*(z*sqrt(snr)+snr)))).* exp(-rho*(log(2) - log(1+exp(-2*s*(z*sqrt(snr)+snr)))));

  f2= @(z) exp(-z.^2/2)/sqrt(2*pi) .*((log(2) - log(1+exp(-2*s*(z*sqrt(snr)+snr)))).^2).* exp(-rho*(log(2) - log(1+exp(-2*s*(z*sqrt(snr)+snr)))));
  
  
  K0=integral(f0, -Zmax, Zmax);
  K1=integral(f1, -Zmax, Zmax);
  K2=integral(f2, -Zmax, Zmax);
  E0=-log(K0);
  E1=K1/K0;
  E2= (-K2*K0+K1^2)/K0^2;
end
