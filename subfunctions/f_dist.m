function [ s_out ] = f_dist(N_trials, p, cd_s, up_s, dist_pj, dist_vec )
%f_dist Summary of this function goes here
%   Detailed explanation goes here
%% arguments check

if nargin < 5
    disp('ERROR: not enough arguments')
    return
elseif nargin > 6
    disp('ERROR: too many arguments')
    return
elseif (nargin == 6) && isempty(dist_vec)
    disp('WARNING: distractor input weights not specified, reverted to Sample mode')
    dist_vec = cd_sample;
end

%% Unpack Params
dt = p.dt;
N = p.N;
ramp_train = p.ramp_train;
tau = p.tau;
f0 = p.f0;
beta0 = p.beta0;
theta0 = p.theta0;
M = p.M;
eff_dt = p.eff_dt;
h = p.h;
trg_left = p.des_out_left;
trg_right = p.des_out_right;
fr_left_norm = p.des_r_left_norm;
fr_right_norm = p.des_r_right_norm;
amp_chirp = p.amp_chirp;
fr_smooth = p.fr_smooth;
ramp_bsln = p.ramp_bsln;
t_ramp_start = p.t_ramp_start;
t_stim_interval = p.t_stim_interval;
endpoint = p.endpoint;
T_test = p.T_test;
ramp_dur = p.ramp_dur;
sigma_noise = p.sigma_noise;
ramp_mean = p.ramp_mean;
ramp_sigma = p.ramp_sigma;
amp_dist = p.amp_dist;
sigma_dist = p.sigma_dist;
dur_dist = p.dur_dist;
t_dist_e = p.t_dist_e;
t_dist_l = p.t_dist_l;
t_sample_start = p.t_sample_start;
t_sample_end = p.t_sample_end;

%% from coding direction
cd_late_delay = cd_s.cd_late_delay;

RNN_fr_cd = cd_s.RNN_fr_cd;
N_trials_cd = cd_s.N_trials_cd;

mean_left_cd = cd_s.mean_left_cd;
mean_right_cd = cd_s.mean_right_cd;
std_left_cd = cd_s.std_left_cd;
std_right_cd = cd_s.std_right_cd;

%% from unperturbed
proj_mean = up_s.proj_mean;
r_up_proj = up_s.r_up_proj;
correct_tri_right_up = up_s.correct_tri_right_up;
correct_tri_left_up = up_s.correct_tri_left_up;
N_trials_up = up_s.N_trials_up;

%%
disp(['Simulate the RNN with distractors.']);

simtime_test = [0:dt:T_test-dt];
simtime_test_len = length(simtime_test);

%% specify projection of distractor input

switch dist_pj
    case 's'
        disp('Input along SAMPLE mode')
        w_pert = M(:,N+2);
    case 'cd'
        disp('Input along CHOICE mode')
        w_pert = cd_late_delay;
    case 'cstm'
        disp('Input along CUSTOM mode')
        w_pert = dist_vec;
end

%% from cd mode
cd_late_delay = cd_s.cd_late_delay;

%% Duration and noise

noise_sigma_eff = sqrt(dt).*sigma_noise./tau; %% define the effective variance of the fast noise, divided by tau

simtime_test = [0:dt:T_test-dt];
simtime_test_len = length(simtime_test);

%% Initial distribution of input weights

inp_chirp_temp = zeros(simtime_test_len,1);
inp_chirp_temp([500/dt+1:650/dt, 1350/dt+1:1500/dt]) = 1;
inp_chirp_temp = smooth(inp_chirp_temp, fr_smooth);

r_in = zeros(3,simtime_test_len);
r_in(1,:) = amp_chirp.*[inp_chirp_temp];

x_init = mean([mean(trg_left(:,1:10),2), mean(trg_right(:,1:10),2)],2); %% mean initial conditions
init_sigma = 0.05; %% variance of initial conditions noise

%% early distractor
inp_dist_vec = zeros(T_test,1);
inp_dist_vec(t_dist_e + 1:t_dist_e + dur_dist) = amp_dist;

%% test with distractors

r_vec_early = [];

for i=1:N_trials
    disp(['trial # ', num2str(i)])
    
    % Ramping input %
    inp_ramp_test = zeros(T_test,1);
    inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./ramp_dur;
    r_in(3,1:t_ramp_start+ramp_dur-1) = ramp_mean.*ramp_train.*inp_ramp_test(1:t_ramp_start+ramp_dur-1).*...
        (1+ramp_sigma*randn) + ramp_bsln;
    r_in(3,t_ramp_start+ramp_dur:end) = r_in(3,t_ramp_start+ramp_dur-1);
    r_in(3,:) = smooth(r_in(3,:), fr_smooth);
    
    % Stimulus input %
    inp_stim_temp = zeros(simtime_test_len,1);
    r_in(2,:) = inp_stim_temp;
    
    inp_dist_temp = smooth(inp_dist_vec.*(1+sigma_dist*randn),fr_smooth);
    inp_dist_test = circshift(inp_dist_temp,randi([-10 10],1)/dt);
    
    % Firing rate for storage
    r_vec_early{i} = zeros(N,simtime_test_len);
    
    
    % initial conditions currents & firing rates
    x = x_init.*(1 + init_sigma.*randn);
    x(N+1:N+3) = 0;
    r = f0./(1.0 + exp(-beta0.*(x(1:N)-theta0)));
    
    % Integrate the dynamics
    ti = 0;
    for up_s = simtime_test
        ti = ti+1;
        
        x = x + (-x + M*[r;r_in(:,ti)] + w_pert.*inp_dist_test(ti) + h)*eff_dt + [noise_sigma_eff*randn(N,1);zeros(3,1)];
        r = f0./(1.0+exp(-beta0.*(x(1:N)-theta0)));
        
        r_vec_early{i}(:,ti) = r; % save the firing rates
    end
end  


%% divide trials and store matrix of firing rates

RNN_fr_dist_e = cat(3,r_vec_early{:});

%% Project trials on CD mode
r_proj_early = [];

for i = 1:N_trials
    r_proj_early(i,:) = squeeze(RNN_fr_dist_e(:,:,i))'*cd_late_delay;
end

%% find aberrant trials
for k = 1:endpoint-1500
    kk = k+1500;

    aberrant_e_dist_temp{k} = find((vecnorm(squeeze(RNN_fr_dist_e(:,kk,:))) >...
    (mean_right_cd(k) + 6*std_right_cd(k))) & (vecnorm(squeeze(RNN_fr_dist_e(:,kk,:))) >...
    (mean_left_cd(k) + 6*std_left_cd(k))));

end
%%
[aberrant_e_dist,~] = unique([aberrant_e_dist_temp{:}]);

%% classify trials early delay

correct_tri_e_temp = find(r_proj_early(:,endpoint)<proj_mean(endpoint));

error_tri_e_temp = find(r_proj_early(:,endpoint)>proj_mean(endpoint));

correct_tri_e = setdiff(correct_tri_e_temp, aberrant_e_dist);
error_tri_e = setdiff(error_tri_e_temp, aberrant_e_dist);

%% late distractor
inp_dist_vec = zeros(T_test,1);
inp_dist_vec(t_dist_l + 1:t_dist_l + dur_dist) = amp_dist;

%% test with late distractor

r_vec_late = [];

for i=1:N_trials
    disp(['trial # ', num2str(i)])
    
    % Ramping input %
    inp_ramp_test = zeros(T_test,1);
    inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./ramp_dur;
    r_in(3,1:t_ramp_start+ramp_dur-1) = ramp_mean.*ramp_train.*inp_ramp_test(1:t_ramp_start+ramp_dur-1).*...
        (1+ramp_sigma*randn) + ramp_bsln;
    r_in(3,t_ramp_start+ramp_dur:end) = r_in(3,t_ramp_start+ramp_dur-1);
    r_in(3,:) = smooth(r_in(3,:), fr_smooth);
    
    % Stimulus input %
    inp_stim_temp = zeros(simtime_test_len,1);
    r_in(2,:) = inp_stim_temp;
    
    inp_dist_temp = smooth(inp_dist_vec.*(1+sigma_dist*randn),fr_smooth);
    inp_dist_test = circshift(inp_dist_temp,randi([-10 10],1)/dt);
    
    % Firing rate for storage
    r_vec_late{i} = zeros(N,simtime_test_len);
    
    
    % initial conditions currents & firing rates
    x = x_init.*(1 + init_sigma.*randn);
    x(N+1:N+3) = 0;
    r = f0./(1.0 + exp(-beta0.*(x(1:N)-theta0)));
    
    % Integrate the dynamics
    ti = 0;
    for up_s = simtime_test
        ti = ti+1;
        
        x = x + (-x + M*[r;r_in(:,ti)] + w_pert.*inp_dist_test(ti) + h)*eff_dt + [noise_sigma_eff*randn(N,1);zeros(3,1)];
        r = f0./(1.0+exp(-beta0.*(x(1:N)-theta0)));
        
        r_vec_late{i}(:,ti) = r; % save the firing rates
    end
end  


%% divide trials and store matrix of firing rates

RNN_fr_dist_l = cat(3,r_vec_late{:});

%%
r_proj_late = [];

for i = 1:N_trials
    r_proj_late(i,:) = squeeze(RNN_fr_dist_l(:,:,i))'*cd_late_delay;
end

%% find aberrant trials
for k = 1:endpoint-1500
    kk = k+1500;

    aberrant_l_dist_temp{k} = find((vecnorm(squeeze(RNN_fr_dist_l(:,kk,:))) >...
    (mean_right_cd(k) + 6*std_right_cd(k))) & (vecnorm(squeeze(RNN_fr_dist_l(:,kk,:))) >...
    (mean_left_cd(k) + 6*std_left_cd(k))));

end
%%
[aberrant_l_dist,~] = unique([aberrant_l_dist_temp{:}]);

%% classify trials early delay

correct_tri_l_temp = find(r_proj_late(:,endpoint)<proj_mean(endpoint));
error_tri_l_temp = find(r_proj_late(:,endpoint)>proj_mean(endpoint));

correct_tri_l = setdiff(correct_tri_l_temp, aberrant_l_dist);
error_tri_l = setdiff(error_tri_l_temp, aberrant_l_dist);

%% Average projection error trials

contra_colors{1} = [0, 0, 255]./255;
contra_colors{2} = [140, 140, 254]./255;
contra_colors{3} = [0, 230, 255]./255;
contra_colors{4} = [0, 51, 102]./255;

ipsi_colors{1} = [255, 0, 0]./255;
ipsi_colors{2} = [255, 128, 77]./255;
ipsi_colors{3} = [231, 207, 13]./255;
ipsi_colors{4} = [179, 51, 26]./255;
ipsi_colors{5} = [255, 0, 255]./255;

y_axis = [-0.2, 1.2];

%% projection error trials left - positive weights
t_vec = ([1:5000]-3500)/1000;

figure
hold on

plot(t_vec, (mean(r_proj_early(error_tri_e,:))-mean(r_proj_early(correct_tri_e,200)))...
    ./(mean(r_up_proj(correct_tri_right_up + N_trials_up/2, endpoint))-mean(r_proj_early(correct_tri_e,200))),...
    'Color',ipsi_colors{2},'Linewidth',3, 'Linestyle',':');
plot(t_vec, (mean(r_proj_late(correct_tri_l,:))-mean(r_proj_early(correct_tri_e,200)))...
    ./(mean(r_up_proj(correct_tri_right_up + N_trials_up/2,endpoint))-mean(r_proj_early(correct_tri_e,200))),...
    'Color',ipsi_colors{3},'Linewidth',3);

plot(t_vec, (mean(r_up_proj(correct_tri_right_up + N_trials_up/2, 1:length(t_vec)))-mean(r_proj_early(correct_tri_e,200)))...
    ./(mean(r_up_proj(correct_tri_right_up + N_trials_up/2, endpoint))-mean(r_proj_early(correct_tri_e,200))),...
    'b','Linewidth',3);
plot(t_vec, (mean(r_up_proj(correct_tri_left_up,1:length(t_vec)))-mean(r_proj_early(correct_tri_e,200)))...
    ./(mean(r_up_proj(correct_tri_right_up + N_trials_up/2,endpoint))-mean(r_proj_early(correct_tri_e,200))),...
    'r','Linewidth',3);

plot(t_sample_start*ones(2,1), y_axis, 'k--','Linewidth',2)
plot(t_sample_end*ones(2,1), y_axis, 'k--','Linewidth',2)

plot(zeros(2,1), y_axis, 'k--','Linewidth',2)
xlim([-3.4,1.6])
title('CD proj. - Early dist. - error tri left')


%% Performance Arseny's plot
    
t_stim_vec_contra = [-5, -3.5, -1.6, -0.8];
performance_vec_left = [length(correct_tri_left_up)/(N_trials_up/2), length(error_tri_e)/(N_trials),...
    length(error_tri_l)/(N_trials)];
performance_vec_right = [length(correct_tri_right_up)/(N_trials_up/2), length(correct_tri_e)/(N_trials),...
    length(correct_tri_l)/(N_trials)];

figure
hold on
plot([-2.5, t_stim_vec_contra(3:end)],[100*performance_vec_right(1), 100*performance_vec_left(2:end)],'k','Linewidth',3)
scatter([-2.5, t_stim_vec_contra(3:end)],[100*performance_vec_right(1), 100*performance_vec_left(2:end)], 20,'k','Linewidth',3)
plot(ones(2,1)*t_sample_end,[-1,105],'k-','Linewidth',2)
xticks([t_stim_vec_contra(1:2), -2.5, t_stim_vec_contra(3:4), 0])
xticklabels({'Control', num2str(t_stim_vec_contra(2)), num2str(-2.5), num2str(t_stim_vec_contra(3)), num2str(t_stim_vec_contra(4)), num2str(0)})
ylabel('Performance (%)')
ylim([0,105])
xlim([-3,0])
box off
set(gca,'fontname','Arial','color','w','fontsize',18)


%% save output space

s_out.RNN_fr_dist_e = single(RNN_fr_dist_e);
s_out.RNN_fr_dist_l = single(RNN_fr_dist_l);

s_out.error_tri_e = error_tri_e;
s_out.error_tri_l = error_tri_l;
s_out.correct_tri_e = correct_tri_e;
s_out.correct_tri_l = correct_tri_l;
s_out.correct_tri_e = correct_tri_e;
s_out.correct_tri_l = correct_tri_l;

s_out.r_proj_early = r_proj_early;
s_out.r_proj_late = r_proj_late;
