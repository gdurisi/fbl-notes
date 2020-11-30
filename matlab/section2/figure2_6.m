% Generate Figure 2.6: SNR (dB) vs  error probability for normal approximations,
% achievability and converse bounds for the bi-AWGN channel. Here, R = 1/2 and n = 128.

DEBUG = 1;

% System parameters:
snr_db_vec  = 1:0.25:5;
snr_vec     = 10.^(snr_db_vec/10);
R_bits      = 1/2;
R           = R_bits*log(2);
epsilon     = logspace(-8,0,50);
n           = 128;

%Initializations:
s_vec = 0.1:0.1:1;
eps_NA = nan(size(snr_vec));
eps_NA_ref = nan(size(snr_vec));
eps_saddle_rcu = nan(size(snr_vec));
eps_saddle_mc = nan(size(snr_vec));

for ii = 1:length(snr_vec)
    snr = snr_vec(ii);
    
    %% Normal approximation (Everything computed in bits!)
    f = @(z) exp(-z.^2/2)/sqrt(2*pi) .* (log2(2) - log2(1+exp(-2*snr-2*sqrt(snr)*z)));
    Zmin = -10; Zmax = 10;
    C = integral(f, Zmin, Zmax);
    
    fV = @(z) exp(-z.^2/2)/sqrt(2*pi) .* (log2(2) - log2(1+exp(-2*snr-2*sqrt(snr)*z))).^2;
    V = integral(fV, Zmin, Zmax) - (C)^2;
    eps_NA(ii) = qfunc((C-R_bits)/sqrt(V/n));
    eps_NA_ref(ii) = qfunc((C-R_bits+0.5*log2(n)/n)/sqrt(V/n));
    
    %% Achievability
    eps_saddle_rcu(ii) = rcu_saddle_biawgn(snr,R,n);
    
    %% Converse
    eps_saddle_mc(ii) = metaconverse_saddle_biawgn_opt_s(snr,R,n);
    
end

if DEBUG == 1
    %% Figures:
    semilogy(snr_db_vec,eps_NA,'--k');hold on
    plot(snr_db_vec,eps_NA_ref,'k')
    plot(snr_db_vec,eps_saddle_rcu,'c')
    plot(snr_db_vec,eps_saddle_mc,'m')
    ylim([1e-8 1])
    legend('RCUs saddle','RCU saddle','Verdu-Han saddle','metaconverse saddle')
else
    %% Save Files
    fileName_NA = ['NA_eps_vs_SNR_n' num2str(n) '_R' num2str(R_bits) '.txt'];
    fileName_NA_ref = ['NA_refined_eps_vs_SNR_n' num2str(n) '_R' num2str(R_bits) '.txt'];
    fileName_RCU = ['RCU_saddle_eps_vs_SNR_n' num2str(n) '_R' num2str(R_bits) '.txt'];
    fileName_MC = ['Metaconverse_saddle_eps_vs_SNR_n' num2str(n) '_R' num2str(R_bits) '.txt'];
    T_NA = table(snr_db_vec', eps_NA');
    T_NA_ref = table(snr_db_vec', eps_NA_ref');
    T_RCU = table(snr_db_vec', eps_saddle_rcu');
    T_MC = table(snr_db_vec', eps_saddle_mc');
    writetable(T_NA, fileName_NA, 'WriteVariableNames',false, 'Delimiter',  ' ')
    writetable(T_NA_ref, fileName_NA_ref, 'WriteVariableNames',false, 'Delimiter',  ' ')
    writetable(T_RCU, fileName_RCU, 'WriteVariableNames',false, 'Delimiter',  ' ')
    writetable(T_MC, fileName_MC, 'WriteVariableNames',false, 'Delimiter',  ' ')
end
