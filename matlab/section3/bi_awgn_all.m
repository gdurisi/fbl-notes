% This script presents different implementations of finte-blocklength upper
% and lower bounds on the minimum error probability achievable using
% channel codes of length n and a given snr vs the coding rate R. 

DEBUG = 1;

snr_db  = 0.189;
snr     = 10^(snr_db/10);
R_bits  = linspace(0.1,0.5,50);
R       = R_bits*log(2);
epsilon = logspace(-8,0,50);
n       = 128;
N       = 1e6;

%% Capacity bi-AWGN 
f = @(z) exp(-z.^2/2)/sqrt(2*pi) .* (log2(2) - log2(1+exp(-2*snr-2*sqrt(snr)*z)));
Zmin = -10; Zmax = 10;
C = integral(f, Zmin, Zmax);

%% Normal approximations
fV = @(z) exp(-z.^2/2)/sqrt(2*pi) .* (log2(2) - log2(1+exp(-2*snr-2*sqrt(snr)*z))).^2;
V = integral(fV, Zmin, Zmax) - (C)^2;
R_NA = C - sqrt(V/n)*qfuncinv(epsilon');
R_NA_refined = C - sqrt(V/n)*qfuncinv(epsilon') + 0.5*log2(n)/n;

%% Achievability
s_vec = 0.1:0.1:1;
eps_rcus= nan(length(s_vec),length(R)); eps_saddle_rcus = eps_rcus;
for ss=1:length(s_vec)
    s = s_vec(ss);
    eps_saddle_rcus(ss,:) = rcus_saddle_biawgn_fixed_s(R,n,snr,s);
    eps_rcus(ss,:)        = rcus_biawgn_fixed_s(R,n,snr,s,N);
end
if size(eps_rcus,1)>1
    [eps_saddle_rcus,best_s_rcus_saddle] = min(eps_saddle_rcus);
    [eps_rcus,best_rcus_s] = min(eps_rcus);
end
eps_saddle_rcu        = rcu_saddle_biawgn(snr,R,n);

%% Converse
s_vec = 0.1:0.1:1; % If s=1, capacity achieving output dist.
eps_vh= nan(length(s_vec),length(R)); eps_saddle_vh = eps_vh; % Only needed if optimization over s. 
Rbits_mc = nan(length(s_vec),length(epsilon)); %eps_mc = Rbits_mc;
for ss=1:length(s_vec)
    s = s_vec(ss);
    eps_saddle_vh(ss,:) = vh_metaconverse_saddle_biawgn_fixed_s(R,n,snr,s);  
    eps_vh(ss,:)        = vh_metaconverse_biawgn_fixed_s(R,n,snr,s,N); 
    Rbits_mc(ss,:) = metaconverse_biawgn_fixed_s(epsilon,n,snr,s,N);
    %eps_mc(ss,:) = metaconverse_rate(R,n,snr,s,N);
end
if size(Rbits_mc,1)>1
    [eps_saddle_vh,best_s_vh_saddle] = max(eps_saddle_vh);
    [eps_vh,best_s_vh] = max(eps_vh);
    [Rbits_mc,best_s_mc] = min(Rbits_mc);
    %[eps_mc,best_s_mc_rate] = max(eps_mc);
end
eps_saddle_mc = metaconverse_saddle_biawgn_opt_s(snr,R,n);

if DEBUG == 1
%% Figures:
semilogy(R_bits,eps_rcus,'.b'); hold on
semilogy(R_bits,eps_saddle_rcus,'b');hold on
plot(R_bits,eps_saddle_rcu,'c')
semilogy(R_bits,eps_vh,'.r');hold on
plot(R_bits,eps_saddle_vh,'r')
plot(Rbits_mc,epsilon,'.m')
%plot(R_bits,eps_mc,'y')
plot(R_bits,eps_saddle_mc,'m')
plot(R_NA, epsilon, '--k')
plot(R_NA_refined, epsilon, 'k')
plot(ones(size(epsilon))*C,epsilon,'g')
ylim([1e-5 1])
xlim([0.1 0.6])
legend('RCUs','RCUs saddle','RCU saddle','Verdu-Han','Verdu-Han saddle','metaconverse','metaconverse saddle','normal approx','normal approx refined','capacity')
else
    %% Save Files
end

