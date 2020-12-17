% Generate Figure 2.5: SNR (dB) vs  error probability for achievability
% and converse bounds for the bi-AWGN channel. Here, R = 1/2 and n = 128.

DEBUG = 1;

% System parameters:
snr_db_vec  = -1:0.25:5;
snr_vec     = 10.^(snr_db_vec/10);
R_bits      = 1/2;
R           = R_bits*log(2);
epsilon     = logspace(-8,0,50);
n           = 128;

%Initializations:
s_vec = 0.1:0.1:1;
eps_saddle_rcus = nan(size(snr_vec));
eps_saddle_rcu = nan(size(snr_vec));
eps_saddle_vh = nan(size(snr_vec));
eps_saddle_mc = nan(size(snr_vec));

for ii = 1:length(snr_vec)
    snr = snr_vec(ii);
    
    %% Achievability
    
    eps_saddle_rcus_s= nan(size(s_vec));
    for ss = 1:length(s_vec)
        s = s_vec(ss);
        eps_saddle_rcus_s(ss) = rcus_saddle_biawgn_fixed_s(R,n,snr,s);
    end
    [eps_saddle_rcus(ii),best_s_rcus_saddle] = min(eps_saddle_rcus_s);
    eps_saddle_rcu(ii) = rcu_saddle_biawgn(snr,R,n);
    
    %% Converse
    eps_saddle_vh(ii) = vh_metaconverse_saddle_biawgn_fixed_s(R,n,snr,1);
    eps_saddle_mc(ii) = metaconverse_saddle_biawgn_opt_s(snr,R,n);
    
end

if DEBUG == 1
    %% Figures:
    semilogy(snr_db_vec,eps_saddle_rcus,'b');hold on
    plot(snr_db_vec,eps_saddle_rcu,'c')
    plot(snr_db_vec,eps_saddle_vh,'r')
    plot(snr_db_vec,eps_saddle_mc,'m')
    ylim([1e-8 1])
    legend('RCUs saddle','RCU saddle','Verdu-Han saddle','metaconverse saddle')
else
    %% Save Files
    fileName_RCUs = ['RCUs_saddle_eps_vs_SNR_n' num2str(n) '_R' num2str(R_bits) '.txt'];
    fileName_RCU = ['RCU_saddle_eps_vs_SNR_n' num2str(n) '_R' num2str(R_bits) '.txt'];
    fileName_VH = ['Verdu-Han_saddle_eps_vs_SNR_n' num2str(n) '_R' num2str(R_bits) '.txt'];
    fileName_MC = ['Metaconverse_saddle_eps_vs_SNR_n' num2str(n) '_R' num2str(R_bits) '.txt'];
    T_RCUs = table(snr_db_vec', eps_saddle_rcus');
    T_RCU = table(snr_db_vec', eps_saddle_rcu');
    T_VH = table(snr_db_vec', eps_saddle_vh');
    T_MC = table(snr_db_vec', eps_saddle_mc');
    writetable(T_RCUs, fileName_RCUs, 'WriteVariableNames',false, 'Delimiter',  ' ')
    writetable(T_RCU, fileName_RCU, 'WriteVariableNames',false, 'Delimiter',  ' ')
    writetable(T_VH, fileName_VH, 'WriteVariableNames',false, 'Delimiter',  ' ')
    writetable(T_MC, fileName_MC, 'WriteVariableNames',false, 'Delimiter',  ' ')
end
