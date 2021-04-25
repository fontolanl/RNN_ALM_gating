function [ s_out ] = f_up(N_trials, p, cd_s, plot_figs)
%f_unperturbed Summary of this function goes here
%   Detailed explanation goes here

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
amp_stim = p.amp_stim;
sigma_stim = p.sigma_stim;
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

disp(['Simulate the RNN - stimulus only (no distractors).']);

%% from cd mode
cd_late_delay = cd_s.cd_late_delay;

RNN_fr_cd = cd_s.RNN_fr_cd;
N_trials_cd = cd_s.N_trials_cd;
correct_tri_left_cd = cd_s.correct_tri_left_cd;
correct_tri_right_cd = cd_s.correct_tri_right_cd;
mean_left_cd = cd_s.mean_left_cd;
mean_right_cd = cd_s.mean_right_cd;
std_left_cd = cd_s.std_left_cd;
std_right_cd = cd_s.std_right_cd;

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

stm_trials = [zeros(N_trials./2,1); amp_stim.*ones(N_trials./2,1)]; %% Stimulus mask for trials, first half left lick, second half right lick

%% test without distractor

r_vec_up = [];

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
    inp_stim_temp(t_stim_interval) = stm_trials(i).*(1+sigma_stim*randn);
    inp_stim_temp = smooth(inp_stim_temp, fr_smooth);
    r_in(2,:) = inp_stim_temp;
    
    % Firing rate for storage
    r_vec_up{i} = zeros(N,simtime_test_len);
    
    
    % initial conditions currents & firing rates
    x = x_init.*(1 + init_sigma.*randn);
    x(N+1:N+3) = 0;
    r = f0./(1.0 + exp(-beta0.*(x(1:N)-theta0)));
    
    % Integrate the dynamics
    ti = 0;
    for t = simtime_test
        ti = ti+1;
        
        x = x + (-x + M*[r;r_in(:,ti)] + h)*eff_dt + [noise_sigma_eff*randn(N,1);zeros(3,1)];
        r = f0./(1.0+exp(-beta0.*(x(1:N)-theta0)));
        
        r_vec_up{i}(:,ti) = r; % save the firing rates
    end
end  

%% divide trials
left_trials_up = 1:N_trials/2;

right_trials_up = (N_trials/2 + 1):N_trials;

RNN_fr_up = cat(3,r_vec_up{:});

%% All trials projected no distractor

r_up_proj = [];

for i = 1:N_trials
    r_up_proj(i,:) = squeeze(RNN_fr_up(:,:,i))'*cd_late_delay;
end


%% find aberrant trials
for k = 1:endpoint-1500
    kk = k+1500;

    aberrant_right_up_temp{k} = find((vecnorm(squeeze(RNN_fr_up(:,kk,right_trials_up))) >...
    (mean_right_cd(k) + 6*std_right_cd(k))) & (vecnorm(squeeze(RNN_fr_up(:,kk,right_trials_up))) >...
    (mean_left_cd(k) + 6*std_left_cd(k))));

    aberrant_left_up_temp{k} = find((vecnorm(squeeze(RNN_fr_up(:,kk,left_trials_up))) >...
    (mean_right_cd(k) + 6*std_right_cd(k))) & (vecnorm(squeeze(RNN_fr_up(:,kk,left_trials_up))) >...
    (mean_left_cd(k) + 6*std_left_cd(k))));

end
%%
[aberrant_right_up,~] = unique([aberrant_right_up_temp{:}]);

[aberrant_left_up,~] = unique([aberrant_left_up_temp{:}]);


%% find correct trials

proj_mean = mean([RNN_fr_cd(:,:,1)'*cd_late_delay, RNN_fr_cd(:,:,N_trials_cd/2 + 1)'*cd_late_delay],2);

correct_tri_left_up = setdiff(find(r_up_proj(left_trials_up,endpoint)<proj_mean(endpoint)),aberrant_left_up);
correct_tri_right_up = setdiff(find(r_up_proj(right_trials_up,endpoint)>proj_mean(endpoint)),aberrant_right_up);

error_tri_left_up = setdiff(find(r_up_proj(left_trials_up,endpoint)>proj_mean(endpoint)),aberrant_left_up);
error_tri_right_up = setdiff(find(r_up_proj(right_trials_up,endpoint)<proj_mean(endpoint)),aberrant_right_up);

rp_up_left_mean = mean(RNN_fr_up(:,:,correct_tri_left_up),3);

rp_up_right_mean = mean(RNN_fr_up(:,:,N_trials/2 + correct_tri_right_up),3);


%% Correct and error trials projected
if plot_figs
    figure
    hold on
    for i = 1:N_trials
        
        if ismember(i,N_trials/2 + correct_tri_right_up)
            plot(r_up_proj(i,:),'b');
        elseif ismember(i,correct_tri_left_up)
            plot(r_up_proj(i,:),'r');
        end
        
        if ismember(i,N_trials/2 + error_tri_right_up)
            plot(r_up_proj(i,:),'m');
        elseif ismember(i,error_tri_left_up)
            plot(r_up_proj(i,:),'c');
        end
        
        if ismember(i,N_trials/2 + aberrant_right_up)
            plot(r_up_proj(i,:),'k');
        end
        
    end
    plot(proj_mean,'k--','Linewidth',2)
    plot(mean(r_up_proj(N_trials/2 + correct_tri_right_up,:)),'c','Linewidth',2);
    plot(mean(r_up_proj(correct_tri_left_up,:)),'m','Linewidth',2);
    plot(mean(RNN_fr_cd(:,:,correct_tri_left_cd),3)'*cd_late_delay,'m--','Linewidth',3)
    plot(mean(RNN_fr_cd(:,:,N_trials_cd/2 + correct_tri_right_cd),3)'*cd_late_delay,'c--','Linewidth',3)
    plot((500-0)*ones(2,1),[-1,max(r_up_proj(:))+0.1],'k--','Linewidth',2')
    plot((1000-0)*ones(2,1),[-1,max(r_up_proj(:))+0.1],'k--','Linewidth',2')
    plot((1500-0)*ones(2,1),[-1,max(r_up_proj(:))+0.1],'k--','Linewidth',2')
    title('Projection onto CD Late Delay All Trials')
end

%% save output space

s_out.RNN_fr_up = single(RNN_fr_up);

s_out.correct_tri_left_up = correct_tri_left_up;
s_out.correct_tri_right_up = correct_tri_right_up;
s_out.error_tri_left_up = error_tri_left_up;
s_out.error_tri_right_up = error_tri_right_up;
s_out.aberrant_left_up = aberrant_left_up;
s_out.aberrant_right_up = aberrant_right_up;
s_out.proj_mean = proj_mean;

s_out.mean_left_cd = mean_left_cd;
s_out.mean_right_cd = mean_right_cd;

s_out.rp_up_left_mean = rp_up_left_mean;
s_out.rp_up_right_mean = rp_up_right_mean;

s_out.std_left_cd = std_left_cd;
s_out.std_right_cd = std_right_cd;

s_out.r_up_proj = r_up_proj;

s_out.N_trials_up = N_trials;
s_out.right_trials_up = right_trials_up;
s_out.left_trials_up = left_trials_up;

end

