%% load trained RNN data

% input_file_path = 'D:\Dropbox (HHMI)\Force learning\current_code\code_for_sharing\git_paper_ALM_S1\RNN_ALM_gating\input_data\';
% 
% output_file_path = 'D:\Dropbox (HHMI)\Force learning\current_code\code_for_sharing\git_paper_ALM_S1\RNN_ALM_gating\output_data';

input_file_path = 'C:\Users\natne\Dropbox (HHMI)\Force learning\current_code\code_for_sharing\git_paper_ALM_S1\RNN_ALM_gating\input_data\';

output_file_path = 'C:\Users\natne\Dropbox (HHMI)\Force learning\current_code\code_for_sharing\git_paper_ALM_S1\RNN_ALM_gating\output_data\';


file_name = 'input_data_wramp';

load([input_file_path, file_name,'.mat'])

%% params from mat file

params_file_name = 'params_data_wramp';

load([input_file_path, params_file_name,'.mat'])
    
N = params.N;

params.t_ramp_start = 500; %% time at which the ramping input begins
params.t_stim_interval = [1000+1:1400]; %% stimulus presentation epoch

%% Compute coding direction using:
%  f_cd(N_trials, ramp_mean, ramp_sigma, ramp_dur, sigma_noise_cd, params, T_trial, cd_span);

params.T_test = 5000;   % trial duration in ms
params.ramp_dur = 3000;  % duration of ramp input
params.sigma_noise_cd = 100./N;  % std of fast noise
params.ramp_mean = 1.0;  % mean slope of ramping input
params.ramp_sigma = 0.05;  % std of slope of ramping input
params.amp_stim = 1;  % amplitude of stimulus

params.sigma_stim = 0.1;  % amplitude of stimulus
params.endpoint = 3500; % decision time

params.amp_chirp = 1;  % amplitude of chirps

cd_span = 200;  % window used to compuite the cd: last 200ms before the end of the delay period (3.5s)

N_trials_cd = 200;  % total number of trials (left + right)

[ cd_struct ] = f_cd(N_trials_cd, params, cd_span);

%% Plot all trials (used to calculate coding direction)

figure
hold on
for i = 1:N_trials_cd
    if (i<=N_trials_cd/2)&&(ismember(i,cd_struct.correct_tri_left_cd))
        plot(mean(cd_struct.RNN_fr_cd(:,:,i)),'r');
    elseif (i>N_trials_cd/2)&&(ismember(i,cd_struct.correct_tri_right_cd + N_trials_cd/2))
        plot(mean(cd_struct.RNN_fr_cd(:,:,i)),'b');
    elseif (i<=N_trials_cd/2)&&(ismember(i,cd_struct.error_tri_left_cd))
        plot(mean(cd_struct.RNN_fr_cd(:,:,i)),'m');
    elseif (i>N_trials_cd/2)&&(ismember(i,cd_struct.error_tri_right_cd + N_trials_cd/2))
        plot(mean(cd_struct.RNN_fr_cd(:,:,i)),'c');
    else
        plot(mean(cd_struct.RNN_fr_cd(:,:,i)),'k');    
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

% f_unperturbed(ramp, ramp_prefactor, ramp_dur, sigma_ramp, p, T_test, s_cd, N_trials_unperturbed,...
%    sigma_stim, sigma_noise, plot_figs,  stim_shape_in)

params.T_test = 5000;   % trial duration in ms
params.ramp_dur = 3000;  % duration of ramp input
params.sigma_noise = 100./N;  % std of fast noise
params.ramp_mean = 1.0;  % mean slope of ramping input
params.ramp_sigma = 0.05;  % std of slope of ramping input
params.amp_stim = 1; % mean amplitude of stimulus
params.sigma_stim = 0.4; % std amplitude of stimulus
N_trials_up = 200; % Number of trials

plot_fig = 1;

[ up_struct ] = f_up(N_trials_up, params, cd_struct, plot_fig);

%% Generate trials with distractors
% f_distractors(saving folder name, save figures flag, number of trials, parameters struct,...
%     output struct from f_cd, output struct from f_unperturbed, ramp type, ramp duration, ramp std, trial duration (ms),...
%     ramp sope, stimulus amplitude mean, stimulus amplitude standard deviation, stimulus duration, fast noise amplitude,...
%     time of first distractor, time of second distractor, endpoint (ms) to calculate performance,...
%     full vs mini distractor flag, amplitude of full distractor, duration of full distractor,...
%     amplitude of mini distractor, distractor amplitude std, shape of stimulus, input (stim and distr) direction type, vector of input direction (optional))

params.T_test = 5000;   % trial duration in ms
params.ramp_dur = 3000;  % duration of ramp input
params.sigma_noise = 100./N;  % std of fast noise
params.ramp_mean = 1.0;  % mean slope of ramping input
params.ramp_sigma = 0.05;  % std of slope of ramping input

params.amp_dist = 1;
params.sigma_dist = 0.4;
params.dur_dist = 400;

params.t_dist_e = 1900;
params.t_dist_l = 2700;

params.endpoint = 3500;

N_trials_dist = 100;

[ dist_struct ] = f_dist(N_trials_dist, params, cd_struct, up_struct, 's');

%% save structs
save(output_file_path, 'cd_struct', 'up_struct', 'dist_struct', '-v7.3')
