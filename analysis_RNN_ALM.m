clear
close all

%% load trained RNN data

input_file_path = '~/Dropbox (HHMI)/Force learning/current_code/code_for_sharing/trained_RNN_output_data/';

file_name = 'N=1000_g=1.5_a=1_sn=0.1497_ss=0.1_sr=0.05_hin=0_Wc=0.1_Ws=1_Wr=1_22-Jan-2020 13:03:19';

load([input_file_path, file_name,'.mat'])

%%
params_name = ['params_', file_name];
% ramp = 'long';
% ramp = 'short';

cd('/Users/fontolanl10/Dropbox (HHMI)/Force learning/current_code/code_for_sharing/RNN_analysis')


%% params from mat file

if isfile(params_name)
    load([params_name,'.mat'])
    save_figs = ['Figures_',params_name];

else
    t_stim_start = ((1000-3500)/1000)./dt;
    t_stim_end = ((1400-3500)/1000)./dt;
    
    t_sample_start = ((500-3500)/1000)./dt;
    t_sample_end = ((1500-3500)/1000)./dt;
    
    t_dist_early_start = ((1900-3500)/1000)./dt;
    t_dist_early_end = ((2300-3500)/1000)./dt;
    
    t_dist_late_start = ((2700-3500)/1000)./dt;
    t_dist_late_end = ((3100-3500)/1000)./dt;
            
    params.tau = tau;
    params.dt = dt;
    params.T = T;
    params.N = N;
    params.inp_chirp_temp = inp_chirp_temp;
    params.inp_stim_temp = inp_stim_temp;
    params.ramp_train = ramp_train;
    params.f0 = f0;
    params.beta0 = beta0;
    params.theta0 = theta0;
    params.W = W;
    params.simtime = simtime;
    params.simtime_len = simtime_len;
    params.b = b;
    params.eff_dt = eff_dt*dt;
    params.des_r_left_norm = des_r_left_norm;
    params.des_r_right_norm = des_r_right_norm;
    % params.des_r_left_norm_cell = des_r_left_norm_cell;
    % params.des_r_right_norm_cell = des_r_right_norm_cell;
    params.des_r_left_norm_cell = des_r_left;
    params.des_r_right_norm_cell = des_r_right;
%     params.des_r_left_norm = des_rt_left;
%     params.des_r_right_norm = des_rt_right;
    params.des_out_left = des_out_left;
    params.des_out_right = des_out_right;
    params.idx_sorted = idx_sorted;
    params.idx_sorted_l = idx_sorted_l;
    params.fr_smooth = fr_smooth;    
    params.fr_win = fr_win;
    params.stim_sigma = stim_sigma;
    params.t_sample_start = t_sample_start;
    params.t_sample_end = t_sample_end;    
    
    % obtain firing rates
    
    des_r_left_temp = f0./(1.0 + exp(-beta0.*(des_out_left-theta0)));
    des_r_right_temp = f0./(1.0 + exp(-beta0.*(des_out_right-theta0)));
    
    params.des_r_left = des_r_left_temp;
    params.des_r_right = des_r_right_temp;
    
    % create directory where to save figures
    mkdir(['Figures_',params_name]);
    save_figs_folder = ['Figures_',params_name];
    
    save([params_name,'.mat'],'params')   

end

%% Plot distribution of spike rates and average spike rate

figure
histogram([mean(des_r_left(:,340:end),2);mean(des_r_right(:,340:end),2)])
ylabel('No. of cells')
xlabel('Spike rate (Hz)')
set(gca,'fontname','Arial','color','w','fontsize',18)

figure
scatter(mean(des_r_left(:,340:end),2)-mean(des_r_left(:,1:60),2),mean(des_r_right(:,340:end),2)-mean(des_r_right(:,1:60),2))
xlabel({'Spike rate lick left trials'; 'baseline subtracted)'})
ylabel({'Spike rate lick right trials'; 'baseline subtracted)'})
set(gca,'fontname','Arial','color','w','fontsize',18)

%%
figure
hold on
scatter([mean(des_psth_left(:,340:end),2);mean(des_psth_right(:,340:end),2)], [mean(des_r_left(:,340:end),2);mean(des_r_right(:,340:end),2)])
plot(0:90,0:90,'k')
xlabel({'Spike rate lick left trials'; 'baseline subtracted)'})
ylabel({'Spike rate lick right trials'; 'baseline subtracted)'})
set(gca,'fontname','Arial','color','w','fontsize',18)


%% Compute coding direction (Choice mode) 
% output_struct = f_cd(number of trials, ramp type, ramp slope, ramp duration,
%                     ramp standard deviation, fast noise amplitude, shape of stimulus, 
%                     parameters struct, trial duration (ms), cd_span)

ramp = 'delay';

params.ramp_bsln = 0; 

params.stim_amp = 1;
params.chirp_amp = 1;

ramp_amp = 1.0;
ramp_sigma = 0.0;

stim_shape_in = 'square';
ramp_dur = 3000/dt;
noise_sigma_cd = 0./N;

N_trials_cd = 4;

cd_span = 200/dt; % Choice mode window (always ending at the Go cue)

[ struct_out_cd ] = f_cd(N_trials_cd, ramp, ramp_amp, ramp_dur, ramp_sigma, noise_sigma_cd, stim_shape_in, params, 5000, cd_span);

%% Plot all trials (used to calculate coding direction)

figure
hold on
for i = 1:N_trials_cd
    if (i<=N_trials_cd/2)&&(ismember(i,struct_out_cd.correct_tri_left_nd))
        plot(mean(struct_out_cd.rp_nd_mat_all_cd(:,:,i)),'r');
    elseif (i>N_trials_cd/2)&&(ismember(i,struct_out_cd.correct_tri_right_nd + N_trials_cd/2))
        plot(mean(struct_out_cd.rp_nd_mat_all_cd(:,:,i)),'b');
    elseif (i<=N_trials_cd/2)&&(ismember(i,struct_out_cd.error_tri_left_nd))
        plot(mean(struct_out_cd.rp_nd_mat_all_cd(:,:,i)),'m');
    elseif (i>N_trials_cd/2)&&(ismember(i,struct_out_cd.error_tri_right_nd + N_trials_cd/2))
        plot(mean(struct_out_cd.rp_nd_mat_all_cd(:,:,i)),'c');
    else
        plot(mean(struct_out_cd.rp_nd_mat_all_cd(:,:,i)),'k');
        
    end
    
end
ylabel('Spike Rate (Hz)')
xlabel('Time to Go cue (s)')
title('Mean Firing rate (Network)')
% xlim([0,3500])
title('All trials')

%% Generate trials without distractors
%   f_unperturbed(number of trials, ramp type, ramp slope, ramp duration, ramp std,
%                 parameters struct, trial duration (ms), output struct from f_cd,...
%                 stimulus amplitude standard deviation, fast noise amplitude,
%                 plot figures flag,  shape of stimulus)

params.ramp_bsln = 0; 
params.stim_amp = 1.0;
params.chirp_amp = 1.0;
N_trials_unperturbed = 10;
noise_sigma = 100./N;

ramp_amp = 1.0;
ramp_sigma = 0.05;

dist_sigma = 0.1;

[ struct_out_unpert ] = f_unperturbed(N_trials_unperturbed, 'delay', ramp_amp, ramp_dur, ramp_sigma, params, 5000,...
    struct_out_cd, dist_sigma, noise_sigma, 1, stim_shape_in);

%% Generate trials with distractors
% f_distractors(saving folder name, save figures flag, number of trials, parameters struct,...
%     output struct from f_cd, output struct from f_unperturbed, ramp type, ramp duration, ramp std, trial duration (ms),...
%     ramp sope, stimulus amplitude mean, stimulus amplitude standard deviation, stimulus duration, fast noise amplitude,...
%     time of first distractor, time of second distractor, endpoint (ms) to calculate performance,...
%     full vs mini distractor flag, amplitude of full distractor, duration of full distractor,...
%     amplitude of mini distractor, distractor amplitude std, shape of stimulus, input (stim and distr) direction type, vector of input direction (optional))

N_trials_distr = 100;

ramp_slope = 1.0;
ramp_dur = 3000/dt;
% sigma_ramp  = 0.1;
ramp_sigma  = 0.05;

stim_amp = 1.0;
% sigma_stim = 0.4;
stim_sigma = 0.1;

stim_dur = 400/dt;

% sigma_noise = 100./N;
noise_sigma = 100./N;

t_dist_1 = 1900/dt;
t_dist_2 = 2700/dt;

full = 1;
full_dist_amp = 1.0;
full_dist_dur = 400/dt;
mini_amp = 1.0;
% sigma_dist = 0.4;
dist_sigma = 0.1;

endpoint = 3500/dt;

[ struct_out_distr ] = f_distractors(save_figs_folder, 1, N_trials_distr, params, struct_out_cd, struct_out_unpert,...
    'delay', ramp_dur, ramp_sigma, 5000,...
    ramp_slope, stim_amp, stim_sigma, stim_dur, noise_sigma, t_dist_1, t_dist_2, endpoint,...
    full, full_dist_amp, full_dist_dur, mini_amp, dist_sigma, stim_shape_in , 's');

%% Explore proportion of switching trials as a function of ramping slope
% f_distractors(saving folder name, save figures flag, number of trials, parameters struct,...
%     output struct from f_cd, output struct from f_unperturbed, ramp type, ramp duration, ramp std, trial duration (ms),...
%     ramp sope, stimulus amplitude mean, stimulus amplitude standard deviation, stimulus duration, fast noise amplitude,...
%     time of first distractor, time of second distractor, endpoint (ms) to calculate performance,...
%     full vs mini distractor flag, amplitude of full distractor, duration of full distractor,...
%     amplitude of mini distractor, distractor amplitude std, shape of stimulus, input (stim and distr) direction type, vector of input direction (optional))

N_trials_distr = 50;

ramp_dur = 3000;
ramp_sigma  = 0.05;
ramp_slope = 1.00;

stim_amp = 1.0;
stim_sigma = 0.0;
stim_dur = 400;

noise_sigma = 100./N;


full = 1;
% amp_full = 1.0;
amp_full = 1.0;

full_dist_dur = 400;
mini_amp = 1.0;
dist_sigma = 0.0;

endpoint = 3500;
%
t_dist_1 = 1900;
tic
ramp_vec = [0.9:0.05:1.2];
for k = 1:length(ramp_vec)
    ramp_slope = ramp_vec(k);
    [ struct_out_distr ] = f_distractors(save_figs_folder, 1, params, struct_out_cd, struct_out_unpert,...
        'delay', ramp_dur, ramp_sigma, 5000,...
        ramp_slope, stim_amp, stim_sigma, stim_dur, noise_sigma, t_dist_1, t_dist_2, endpoint,...
        N_trials_distr, full, full_dist_amp, full_dist_dur, mini_amp, dist_sigma, stim_shape_in , 's');
    
    perf_ramp(k) = length(struct_out_distr.error_tri_left{1});
    
end

%% Ramping vs proportion of switching trials

figure
hold on
plot(ramp_vec, perf_ramp./N_trials_distr,'Linewidth',2)
xlabel('Ramping input slope')
ylabel('Proportion of switching trials')
