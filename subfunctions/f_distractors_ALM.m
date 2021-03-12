function [ struct_out ] = f_distractors_ALM(saving_folder, save_fig_flag, N_trials_distr, p, s_cd, s_unpert,...
    ramp, ramp_dur, ramp_sigma, T_test,...
    ramp_prefactor, stim_amp, stim_sigma, stim_dur, noise_sigma, t_dist_1, t_dist_2, endpoint,...
    full, amp_full, dur_full, amp_mini, dist_sigma, stim_shape_in, dist_pj, dist_vec )

%% Params
dt = p.dt;
N = p.N;
tau = p.tau;
% inp_chirp_temp = p.inp_chirp_temp;
% inp_stim_temp = p.inp_stim_temp;
ramp_train = p.ramp_train;
f0 = p.f0;
beta0 = p.beta0;
theta0 = p.theta0;
W = p.W;
eff_dt = p.eff_dt;
simtime_len = p.simtime_len;
b = p.b;
% des_rt_left = p.des_rt_left;
% des_rt_right = p.des_rt_right;
des_out_left = p.des_out_left;
des_out_right = p.des_out_right;

fr_smooth = p.fr_smooth;
ramp_bsln = p.ramp_bsln;
chirp_amp = p.chirp_amp;

%% from cd mode
cd_late_delay = s_cd.cd_late_delay;
cd_late_delay_neg = s_cd.cd_late_delay_neg;
cd_sample = s_cd.cd_sample;
cd_sample_orth = s_cd.cd_sample_orth;
rp_nd_mat_all_cd = s_cd.rp_nd_mat_all_cd;
% endpoint = s_cd.endpoint;
N_trials_cd = s_cd.N_trials_cd;
% proj_cd = s_cd.proj_cd;

%% from unperturbed
proj_mean = s_unpert.proj_mean;
mean_left_cd = s_unpert.mean_left_cd;
mean_right_cd = s_unpert.mean_right_cd;

std_left_cd = s_unpert.std_left_cd;
std_right_cd = s_unpert.std_right_cd;

rp_nd_mat_all = s_unpert.rp_nd_mat_all;

rp_nd_proj = s_unpert.rp_nd_proj;

N_trials_unperturbed = s_unpert.N_trials_unperturbed;

correct_tri_left_nd = s_unpert.correct_tri_left_nd;
correct_tri_right_nd = s_unpert.correct_tri_right_nd;
error_tri_left_nd = s_unpert.error_tri_left_nd;
error_tri_right_nd = s_unpert.error_tri_right_nd;
aberrant_left_nd = s_unpert.aberrant_left_nd;
aberrant_right_nd = s_unpert.aberrant_right_nd;

right_trials_nd = s_unpert.right_trials_nd;
left_trials_nd = s_unpert.left_trials_nd;

%rp_nd_proj_neg = s_unpert.rp_nd_proj_neg;
%rp_nd_proj_s = s_unpert.rp_nd_proj_s;

%% arguments check

if nargin < 25
    disp('ERROR: not enough arguments')
    return
elseif nargin > 26
    disp('ERROR: too many arguments')
    return
elseif (nargin == 26) && isempty(dist_vec)
    disp('WARNING: distractor input weights not specified, reverted to Sample mode')
    dist_vec = cd_sample;
end

%% specify projection of distractor input
W_cd = W;

switch dist_pj
    case 's'
        disp('Input along SAMPLE mode')
    case 'cd'
        disp('Input along CHOICE mode')
        W_cd(1:N,N+2) = cd_late_delay;
    case 'cdx'
        disp('Input along currents CHOICE mode')
        x_left = p.invtransFun(mean(rp_nd_mat_all_cd(:,:,1:N_trials_cd/2),3));
        x_right = p.invtransFun(mean(rp_nd_mat_all_cd(:,:,N_trials_cd/2+1:end),3));        
        proj_vec_x = x_right - x_left;
%         cd_sample = mean(proj_cd(:,1001:1400),2);
%         cd_sample_x = mean(proj_vec_x(:,1001:1400),2);
        W_cd(1:N,N+2) = mean(proj_vec_x(:,1500/dt:3500/dt),2);
    case 'cstm'
        disp('Input along CUSTOM mode')
        W_cd(1:N,N+2) = dist_vec;
end

%% stimulus shape
switch stim_shape_in
    case 'sinewave'
        stim_shape = abs(sin(0.005*2*pi*[1:400]));
    case 'square'
        stim_shape = 1;      
end

%% test with distractors along choice mode

simtime_test = [0:dt:T_test-dt];
simtime_test_len = length(simtime_test);

t_vec_temp = (simtime_test-3500)./1000;
t_vec = t_vec_temp(1:1/dt:end);

t_stim_start = ((1000-3500)/1000/dt);
t_stim_end = ((1400-3500)/1000)/dt;

t_sample_start = ((500-3500)/1000)/dt;
t_sample_end = ((1500-3500)/1000)/dt;

%
%% Initial distribution of input weights

inp_chirp_temp = zeros(simtime_test_len,1);
inp_chirp_temp([500/dt+1:650/dt, 1350/dt+1:1500/dt]) = 1;
inp_chirp_temp = smooth(inp_chirp_temp, fr_smooth);

r_in = zeros(3,simtime_test_len);

r_in(1,:) = chirp_amp.*[inp_chirp_temp];

save_fig = [saving_folder];


%%

full_dur = dur_full;
mini_dur = 100/dt;

dist_amp = amp_mini*(1-full) + amp_full*full;

%% params to run RNN
params.simtime_test_len = simtime_test_len;
params.simtime_test = simtime_test;

params.simtime_len = simtime_len;
params.N_trials_distr = N_trials_distr;
params.full_dur = full_dur;
params.full = full;

params.mini_dur = mini_dur;
params.dist_amp = dist_amp;

params.stim_shape = stim_shape;
params.stim_amp = stim_amp;

params.ramp_dur = ramp_dur;

params.ramp_bsln = ramp_bsln;

params.des_out_left = des_out_left;
params.des_out_right = des_out_right;

params.W = W;
params.W_cd = W_cd;

params.ramp = ramp;
params.ramp_sigma = ramp_sigma;
params.ramp_prefactor = ramp_prefactor;
params.stim_dur = stim_dur;
params.noise_sigma = noise_sigma;
params.dist_sigma = dist_sigma;
params.ramp_train = ramp_train;
params.fr_smooth = fr_smooth;
params.f0 = f0;
params.theta0 = theta0;
params.beta0 = beta0;
params.dt = dt;
params.b = b;
params.eff_dt = eff_dt;
params.stim_sigma = stim_sigma;

%% effective sigma noise
noise_sigma_eff = sqrt(dt).*noise_sigma./tau;

params.noise_sigma_eff = noise_sigma_eff;

%% early distractor input
params.t_dist = t_dist_1;

[ rp_vec_early, r_in_early ] = run_RNN_left( params, r_in );


%% average

left_trials = 1:N_trials_distr;
right_trials = [];

% right_trials = (N_trials_distr/2 + 1):N_trials_distr;

rp_mat_all_early = cat(3,rp_vec_early{:});

rp_vec_early = [];

%% average

%% early distractor input
params.t_dist = t_dist_2;

[ rp_vec_late, r_in_late ] = run_RNN_left( params, r_in );

rp_mat_all_late = cat(3,rp_vec_late{:});

rp_vec_late = [];
%%
rp_proj_early = [];
rp_proj_early_neg = [];

for i = 1:N_trials_distr
    rp_proj_early(i,:) = squeeze(rp_mat_all_early(:,:,i))'*cd_late_delay;
    rp_proj_early_neg(i,:) = squeeze(rp_mat_all_early(:,:,i))'*cd_late_delay_neg;
end


%% classify trials early delay

correct_tri_left_pj = find(rp_proj_early(left_trials,endpoint)<proj_mean(endpoint));
correct_tri_right_pj = find(rp_proj_early(right_trials,endpoint)>proj_mean(endpoint));

error_tri_left_pj = find(rp_proj_early(left_trials,endpoint)>proj_mean(endpoint));
error_tri_right_pj = find(rp_proj_early(right_trials,endpoint)<proj_mean(endpoint));

%%% find aberrant trials

aberrant_right_early = find((mean(rp_mat_all_early(:,endpoint,right_trials),1) >...
    (mean_right_cd + 4*std_right_cd)) | (mean(rp_mat_all_early(:,endpoint,right_trials),1) <...
    (mean_left_cd - 4*std_left_cd)));

aberrant_left_early = find((mean(rp_mat_all_early(:,endpoint,left_trials),1) >...
    (mean_right_cd + 4*std_right_cd)) | (mean(rp_mat_all_early(:,endpoint,left_trials),1) <...
    (mean_left_cd - 4*std_left_cd)));

correct_tri_left{1} = setdiff(correct_tri_left_pj, aberrant_left_early);
correct_tri_right{1} = setdiff(correct_tri_right_pj, aberrant_right_early);
error_tri_left{1} = setdiff(error_tri_left_pj, aberrant_left_early);
error_tri_right{1} = setdiff(error_tri_right_pj, aberrant_right_early);


%% store trajectories
rp_proj_late = [];
rp_proj_late_neg = [];

for i = 1:N_trials_distr
    rp_proj_late(i,:) = squeeze(rp_mat_all_late(:,:,i))'*cd_late_delay;
    rp_proj_late_neg(i,:) = squeeze(rp_mat_all_late(:,:,i))'*cd_late_delay_neg;
end


%% classify trials late delay

%%% find aberrant trials


correct_tri_left_pj_l = find(rp_proj_late(left_trials,endpoint)<proj_mean(endpoint));
correct_tri_right_pj_l = find(rp_proj_late(right_trials,endpoint)>proj_mean(endpoint));

error_tri_left_pj_l = find(rp_proj_late(left_trials,endpoint)>proj_mean(endpoint));
error_tri_right_pj_l = find(rp_proj_late(right_trials,endpoint)<proj_mean(endpoint));

aberrant_right_late = find((mean(rp_mat_all_late(:,endpoint,right_trials),1) >...
    (mean_right_cd + 4*std_right_cd)) | (mean(rp_mat_all_late(:,endpoint,right_trials),1) <...
    (mean_left_cd - 4*std_left_cd)));

aberrant_left_late = find((mean(rp_mat_all_late(:,endpoint,left_trials),1) >...
    (mean_right_cd + 4*std_right_cd)) | (mean(rp_mat_all_late(:,endpoint,left_trials),1) <...
    (mean_left_cd - 4*std_left_cd)));

correct_tri_left{2} = setdiff(correct_tri_left_pj_l, aberrant_left_late);
correct_tri_right{2} = setdiff(correct_tri_right_pj_l, aberrant_right_late);
error_tri_left{2} = setdiff(error_tri_left_pj_l, aberrant_left_late);
error_tri_right{2} = setdiff(error_tri_right_pj_l, aberrant_right_late);


%% Color scheme

contra_colors{1} = [0, 0, 255]./255;
contra_colors{2} = [140, 140, 254]./255;
contra_colors{3} = [0, 230, 255]./255;
contra_colors{4} = [0, 51, 102]./255;

ipsi_colors{1} = [255, 0, 0]./255;
% ipsi_colors{2} = [255, 128, 77]./255;
ipsi_colors{2} = [129, 129, 129]./255;
ipsi_colors{3} = [9, 15, 16]./255;


% ipsi_colors{3} = [231, 207, 13]./255;
ipsi_colors{4} = [179, 51, 26]./255;
ipsi_colors{5} = [255, 0, 255]./255;

shade_color = [151, 216, 252]./255;


y_axis = [-0.2, 1.2];

y_axis_neg = [-5, 2];

%%
correct_trials_pj = [correct_tri_left_pj; correct_tri_right_pj + left_trials(end)];

correct_trials_pj_l = [correct_tri_left_pj_l; correct_tri_right_pj_l + left_trials(end)];

%% projection onto Choice Mode - all weights

if save_fig_flag
    
    figure
    hold on
    
    plot(t_vec, (mean(rp_proj_early(error_tri_left{1},:))-mean(rp_proj_early(correct_trials_pj,200)))...
        ./(mean(rp_nd_proj(correct_tri_right_nd + N_trials_unperturbed/2,endpoint))-mean(rp_proj_early(correct_trials_pj,200))),...
        'Color',ipsi_colors{2},'Linewidth',3, 'Linestyle','--');
    plot(t_vec, (mean(rp_proj_late(error_tri_left{2},:))-mean(rp_proj_early(correct_trials_pj,200)))...
        ./(mean(rp_nd_proj(correct_tri_right_nd + N_trials_unperturbed/2,endpoint))-mean(rp_proj_early(correct_trials_pj,200))),...
        'Color',ipsi_colors{3},'Linewidth',3, 'Linestyle','--');
    
    plot(t_vec, (mean(rp_nd_proj(correct_tri_right_nd + N_trials_unperturbed/2,1:length(t_vec)))-mean(rp_proj_early(correct_trials_pj,200)))...
        ./(mean(rp_nd_proj(correct_tri_right_nd + N_trials_unperturbed/2,endpoint))-mean(rp_proj_early(correct_trials_pj,200))),...
        'b','Linewidth',3);
    plot(t_vec, (mean(rp_nd_proj(correct_tri_left_nd,1:length(t_vec)))-mean(rp_proj_early(correct_trials_pj,200)))...
        ./(mean(rp_nd_proj(correct_tri_right_nd + N_trials_unperturbed/2,endpoint))-mean(rp_proj_early(correct_trials_pj,200))),...
        'r','Linewidth',3);
    
    plot(t_vec, (mean(rp_proj_early(correct_tri_left{1},:))-mean(rp_proj_early(correct_trials_pj,200)))...
        ./(mean(rp_nd_proj(correct_tri_right_nd + N_trials_unperturbed/2,endpoint))-mean(rp_proj_early(correct_trials_pj,200))),'Color',ipsi_colors{2},'Linewidth',3);

    plot(t_vec, (mean(rp_proj_late(correct_tri_left{2},:))-mean(rp_proj_early(correct_trials_pj,200)))...
        ./(mean(rp_nd_proj(correct_tri_right_nd + N_trials_unperturbed/2,endpoint))-mean(rp_proj_early(correct_trials_pj,200))),'Color',ipsi_colors{3},'Linewidth',3);
    
    plot(t_sample_start*ones(2,1), y_axis, 'k--','Linewidth',2)
    plot(t_sample_end*ones(2,1), y_axis, 'k--','Linewidth',2)

    plot(zeros(2,1), y_axis, 'k--','Linewidth',2)
    xlim([-3.4,1.6])
    ylabel('Choice mode projection (a.u.)')
    xlabel('Time to Go cue (s)')
    title('Robustness of RNN trajectories')
%     print(fullfile(save_fig,['avg_proj_choice_mode_error']), '-dpdf'); %<-Save as jpeg
%     savefig(fullfile(save_fig,'avg_proj_choice_mode_error')); %<-Save as Fig
end

    
%% All trials projected with distractor (sample mode) early

rp_proj_s_early = [];

for i = 1:N_trials_distr
    rp_proj_s_early(i,:) = squeeze(rp_mat_all_early(:,:,i))'*cd_sample_orth;
    
end


%% All trials projected with distractor (sample mode) late

rp_proj_s_late = [];

for i = 1:N_trials_distr
    rp_proj_s_late(i,:) = squeeze(rp_mat_all_late(:,:,i))'*cd_sample_orth;

end

%% Average projection correct trials
% y_axis = [min(mean(rp_nd_proj_s(right_trials_nd,100:length(t_vec))))-0.5,...
%     max(mean(rp_nd_proj_s(right_trials_nd,1:length(t_vec))))+0.2];
% 
% y_axis = [-1,2];
% 
% if save_fig_flag
%     
%     figure
%     hold on
%     
%     plot(t_vec, mean(rp_proj_s_early(correct_tri_left_pj,:)),'Color',ipsi_colors{2},'Linewidth',2);
%     
%     plot(t_vec, mean(rp_proj_s_late(correct_tri_left_pj_l,:)),'Color',ipsi_colors{3},'Linewidth',2);
%     
%     plot(t_vec, mean(rp_proj_s_late(correct_tri_left_pj_l,:)),'Color',ipsi_colors{3},'Linewidth',2);
%     
%     plot(t_vec, mean(rp_nd_proj_s(N_trials_unperturbed/2 + correct_tri_right_nd,1:length(t_vec))),'b','Linewidth',2);
%     plot(t_vec, mean(rp_nd_proj_s(correct_tri_left_nd,1:length(t_vec))),'r','Linewidth',2);
%     
%     plot(t_sample_start*ones(2,1), y_axis, 'k--','Linewidth',2)
%     plot(t_sample_end*ones(2,1), y_axis, 'k--','Linewidth',2)
%     % plot(t_stim_start*ones(2,1), y_axis, 'k--','Linewidth',2)
%     % plot(t_dist_early_start*ones(2,1), y_axis, 'k--','Linewidth',2)
%     % plot(t_dist_early_end*ones(2,1), y_axis, 'k--','Linewidth',2)
%     plot(zeros(2,1), y_axis, 'k--','Linewidth',2)
%     xlim([-3.4,0.1])
%     ylim([-0.6,2.0])
%     xlabel('Time to Go cue (s)')
%     ylabel('Choice mode projection (a.u.)')
%     title('Projection onto CD Sample - Correct trials')
%     % print(fullfile(save_fig,['avg_proj_sample_mode_correct']), '-dpdf'); %<-Save as jpeg
%     % savefig(fullfile(save_fig,'avg_proj_sample_mode_correct')); %<-Save as Fig
% end


%% Average projection error trials
% y_axis = [min(mean(rp_nd_proj_s(right_trials_nd,100:length(t_vec))))-0.5,...
%     max(mean(rp_nd_proj_s(right_trials_nd,1:length(t_vec))))+0.2];
% 
% y_axis = [-1,2];
% 
% if save_fig_flag
%     
%     figure
%     hold on
%     
%     plot(t_vec, mean(rp_proj_s_early(left_trials,:)),'Color',ipsi_colors{2},'Linewidth',2);
% 
%     
%     plot(t_vec, mean(rp_proj_s_late(left_trials,:)),'Color',ipsi_colors{3},'Linewidth',2);
%     
%     plot(t_vec, mean(rp_nd_proj_s(N_trials_unperturbed/2 + correct_tri_right_nd,1:length(t_vec))),'b','Linewidth',2);
%     plot(t_vec, mean(rp_nd_proj_s(correct_tri_left_nd,1:length(t_vec))),'r','Linewidth',2);
%     
%     plot(t_sample_start*ones(2,1), y_axis, 'k--','Linewidth',2)
%     plot(t_sample_end*ones(2,1), y_axis, 'k--','Linewidth',2)
% 
%     plot(zeros(2,1), y_axis, 'k--','Linewidth',2)
%     xlim([-3.4,0.1])
%     ylim([-0.6,2.0])
%     xlabel('Time to Go cue (s)')
%     ylabel('Choice mode projection (a.u.)')
%     title('Projection onto CD Sample - All trials')
%     % print(fullfile(save_fig,['avg_proj_sample_mode_all']), '-dpdf'); %<-Save as jpeg
%     % savefig(fullfile(save_fig,'avg_proj_sample_mode_all')); %<-Save as Fig
% end

%% Performance

correct_vec_left = [length(correct_tri_left_nd)/(N_trials_unperturbed/2), length(correct_tri_left_pj)/(length(left_trials)),...
    length(correct_tri_left_pj_l)/(length(left_trials))];
correct_vec_right = [length(correct_tri_right_nd)/(N_trials_unperturbed/2), length(correct_tri_right_pj)/(length(right_trials)),...
    length(correct_tri_right_pj_l)/(length(right_trials))];

error_vec_left = [length(error_tri_left_nd)/(N_trials_unperturbed/2), length(error_tri_left_pj)/(length(left_trials)),...
    length(error_tri_left_pj_l)/(length(left_trials))];
error_vec_right = [length(error_tri_right_nd)/(N_trials_unperturbed/2), length(error_tri_right_pj)/(length(right_trials)),...
    length(error_tri_right_pj_l)/(length(right_trials))];

aberrant_vec_left = [length(aberrant_left_nd)/(N_trials_unperturbed/2), length(aberrant_left_early)/(length(left_trials)),...
    length(aberrant_left_late)/(length(left_trials))];
aberrant_vec_right = [length(aberrant_right_nd)/(N_trials_unperturbed/2), length(aberrant_right_early)/(length(right_trials)),...
    length(aberrant_right_late)/(length(right_trials))];

%%
if save_fig_flag
    
    figure
    subplot(3,1,1)
    bar([correct_vec_left, correct_vec_right],'k')
    ylim([-0.01,1.01])
    ylabel('Fraction of trials')
    title('Correct')
    xticklabels({'Left unperturbed', 'Left early', 'Left late', 'Right unperturbed', 'Right early', 'Right late'})
    
    subplot(3,1,2)
    bar([error_vec_left, error_vec_right],'g')
    ylim([-0.01,1.01])
    ylabel('Fraction of trials')
    title('Error')
    xticklabels({'Left unperturbed', 'Left early', 'Left late', 'Right unperturbed', 'Right early', 'Right late'})
    
    
    subplot(3,1,3)
    bar([aberrant_vec_left, aberrant_vec_right],'y')
    ylim([-0.01,1.01])
    ylabel('Fraction of trials')
    title('Aberrant')
    xticklabels({'Left unperturbed', 'Left early', 'Left late', 'Right unperturbed', 'Right early', 'Right late'})
    
%     print(fullfile(save_fig,['performance_all']), '-dpdf');
%     savefig(fullfile(save_fig,'performance_all'));
end
%% Performance Arseny's plot
if save_fig_flag
    
    t_stim_vec_contra = [-5, -3.5, -1.6, -0.8];
    performance_vec_left = [length(correct_tri_left_nd)/(N_trials_unperturbed/2), length(error_tri_left{1})/(length(left_trials)),...
        length(error_tri_left{2})/(length(left_trials))];
    performance_vec_right = [length(correct_tri_right_nd)/(N_trials_unperturbed/2), length(correct_tri_right{1})/(length(right_trials)),...
        length(correct_tri_right{2})/(length(right_trials))];
    
    figure
    hold on
    
    plot([-2.5, t_stim_vec_contra(3:end)],[performance_vec_right(1), performance_vec_left(2:end)],'k','Linewidth',3)

    scatter([-2.5, t_stim_vec_contra(3:end)],[performance_vec_right(1), performance_vec_left(2:end)], 10,'k','Linewidth',3)
    
    plot(ones(2,1)*t_stim_end,[-1,105],'k-','Linewidth',2)

    xticks([-2.5, t_stim_vec_contra(3:4), 0])
    xticklabels({num2str(-2.5), num2str(t_stim_vec_contra(3)), num2str(t_stim_vec_contra(4)), num2str(0)})
    ylabel('Proportion of lick right responses')
    xlabel('Time to Go cue (s)')
    ylim([0,1.05])
    xlim([-3,0])
    box off
    set(gca,'fontname','Arial','color','w','fontsize',18)
    
%     print(fullfile(save_fig,['performance_arseny']), '-dpdf');
%     savefig(fullfile(save_fig,'performance_arseny'));
end

%% save output space

struct_out.rp_nd_mat_all = rp_nd_mat_all;

struct_out.correct_vec_left = correct_vec_left;
struct_out.correct_vec_right = correct_vec_right;
struct_out.error_vec_left = error_vec_left;
struct_out.error_vec_right = error_vec_right;
struct_out.aberrant_vec_left = aberrant_vec_left;
struct_out.aberrant_vec_right = aberrant_vec_right;

struct_out.error_tri_left_pj = error_tri_left_pj;
struct_out.error_tri_right_pj = error_tri_right_pj;
struct_out.error_tri_left_pj_l = error_tri_left_pj_l;
struct_out.error_tri_right_pj_l = error_tri_right_pj_l;

struct_out.correct_tri_left_pj = correct_tri_left_pj;
struct_out.correct_tri_right_pj = correct_tri_right_pj;
struct_out.correct_tri_left_pj_l = correct_tri_left_pj_l;
struct_out.correct_tri_right_pj_l = correct_tri_right_pj_l;

struct_out.aberrant_right_nd = aberrant_right_nd;
struct_out.aberrant_left_nd = aberrant_left_nd;
struct_out.aberrant_right_early = aberrant_right_early;
struct_out.aberrant_left_early = aberrant_left_early;
struct_out.aberrant_right_late = aberrant_right_late;
struct_out.aberrant_left_late = aberrant_left_late;

struct_out.correct_tri_left = correct_tri_left;
struct_out.correct_tri_right = correct_tri_right;
struct_out.error_tri_left = error_tri_left;
struct_out.error_tri_right = error_tri_right;

struct_out.correct_trials_pj = correct_trials_pj;
struct_out.correct_trials_pj = correct_trials_pj;

struct_out.rp_proj_early = rp_proj_early;
struct_out.rp_proj_early_neg = rp_proj_early_neg;
struct_out.rp_proj_late = rp_proj_late;
struct_out.rp_proj_late_neg = rp_proj_late_neg;
struct_out.rp_proj_s_early = rp_proj_s_early;
struct_out.rp_proj_s_late = rp_proj_s_late;

struct_out.rp_mat_all_early = rp_mat_all_early;
struct_out.rp_mat_all_late = rp_mat_all_late;

struct_out.r_in_early = r_in_early;
struct_out.r_in_late = r_in_late;

