% Generate Figure 2.7: Normal approximation vs finite-blocklength upper and lower bounds
% on the error probability epsilon as a function of the rate R for different
% values of the blocklength n.

DEBUG = 1;

snr_db  = 0.189;
snr     = 10^(snr_db/10);
R_bits  = linspace(0.1,0.5,50);
R       = R_bits*log(2);
epsilon = logspace(-9,-0.3010,50);
n_vec   = [128 256 512 1024 1e4];

for ii = 1:length(n_vec)
    n = n_vec(ii);
    
    %% Normal approximation
    f = @(z) exp(-z.^2/2)/sqrt(2*pi) .* (log2(2) - log2(1+exp(-2*snr-2*sqrt(snr)*z)));
    Zmin = -9; Zmax = 9;
    C = integral(f, Zmin, Zmax);
    fV = @(z) exp(-z.^2/2)/sqrt(2*pi) .* (log2(2) - log2(1+exp(-2*snr-2*sqrt(snr)*z))).^2;
    V = integral(fV, Zmin, Zmax) - (C)^2;
    eps_NA_refined = qfunc((C-R_bits+0.5*log2(n)/n)/sqrt(V/n));
    
    %% Achievability
    eps_saddle_rcu        = rcu_saddle_biawgn(snr,R,n);
    
    %% Converse
    eps_saddle_mc = metaconverse_saddle_biawgn_opt_s(snr,R,n);
    
    if DEBUG == 1
        semilogy(R_bits, eps_NA_refined, 'k'); hold on
        plot(R_bits,eps_saddle_rcu,'b');
        plot(R_bits,eps_saddle_mc,'r');
        ylim([1e-8 1])
        xlim([0.1 0.6])
    else
        fileName_NA = ['NA_n' num2str(n) '_snrdB' num2str(snr_db) '.txt'];
        fileName_RCU = ['RCU_n' num2str(n) '_snrdB' num2str(snr_db) '.txt'];
        fileName_MC = ['Metaconverse_n' num2str(n) '_snrdB' num2str(snr_db) '.txt'];
        T_NA = table(R', eps_NA_refined');
        T_RCU = table(R', eps_saddle_rcu');
        T_MC = table(R', eps_saddle_mc');
        writetable(T_NA, fileName_NA, 'WriteVariableNames',false, 'Delimiter',  ' ')
        writetable(T_RCU, fileName_RCU, 'WriteVariableNames',false, 'Delimiter',  ' ')
        writetable(T_MC, fileName_MC, 'WriteVariableNames',false, 'Delimiter',  ' ')
    end
end
