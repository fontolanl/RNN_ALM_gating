function [ s_out ] = f_cd(N_trials_cd, p, cd_span)
%f_cd Summary of this function goes here
%   Detailed explanation goes here

%% Unpack Parameters
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
ramp_mean = p.ramp_mean;
ramp_sigma = p.ramp_sigma;
ramp_dur = p.ramp_dur;
sigma_noise = p.sigma_noise_cd;
T_test = p.T_test;
endpoint = p.endpoint;

noise_sigma_eff = sqrt(dt).*sigma_noise./tau; %% define the effective variance of the fast noise, divided by tau

disp(['Running a few trials to compute the CD mode']);

simtime_test = [0:dt:T_test-dt]; %% trial time vector
simtime_test_len = length(simtime_test);


%% Initial conditions

x_init = mean([mean(trg_left(:,1:10),2), mean(trg_right(:,1:10),2)],2); %% mean initial conditions
init_sigma = 0.05; %% variance of initial conditions noise

stm_trials = [zeros(N_trials_cd./2,1); amp_stim.*ones(N_trials_cd./2,1)]; %% Stimulus mask for trials, first half left lick, second half right lick

%% Define the Chirps input

inp_chirp_temp = zeros(simtime_test_len,1);
inp_chirp_temp([500/dt+1:650/dt, 1350/dt+1:1500/dt]) = 1;
inp_chirp_temp = smooth(inp_chirp_temp, fr_smooth);

r_in_cd(1,:) = amp_chirp.*inp_chirp_temp;

%% dynamics

rp_vec_nd = [];

for i=1:N_trials_cd
    disp(['trial # ', num2str(i)])        
    
    % Ramping input %
    inp_ramp_test = zeros(T_test,1);
    inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./ramp_dur;
    r_in_cd(3,1:t_ramp_start+ramp_dur-1) = ramp_mean.*ramp_train.*inp_ramp_test(1:t_ramp_start+ramp_dur-1).*...
        (1+ramp_sigma*randn) + ramp_bsln;
    r_in_cd(3,t_ramp_start+ramp_dur:end) = r_in_cd(3,t_ramp_start+ramp_dur-1);
    r_in_cd(3,:) = smooth(r_in_cd(3,:), fr_smooth);    
    
    % Stimulus input %
    inp_stim_temp = zeros(simtime_test_len,1);
    inp_stim_temp(t_stim_interval) = stm_trials(i).*(1+sigma_stim*randn);
    inp_stim_temp = smooth(inp_stim_temp, fr_smooth);
    r_in_cd(2,:) = inp_stim_temp;

    % Firing rate for storage
    rp_vec_nd{i} = zeros(N,simtime_test_len);
      

    % initial conditions currents & firing rates
    x = x_init.*(1 + init_sigma.*randn);       
    x(N+1:N+3) = 0;  
    r = f0./(1.0 + exp(-beta0.*(x(1:N)-theta0)));
    
    % Integrate the dynamics
    ti = 0;
    for t = simtime_test       
        ti = ti+1;
                
        x = x + (-x + M*[r;r_in_cd(:,ti)] + h)*eff_dt + [noise_sigma_eff*randn(N,1);zeros(3,1)];
        r = f0./(1.0+exp(-beta0.*(x(1:N)-theta0))); 
        
        rp_vec_nd{i}(:,ti) = r; % save the firing rates     
    end
    
end

%% divide trials
left_trials_cd = 1:N_trials_cd/2;

right_trials_cd = N_trials_cd/2 + 1:N_trials_cd;

RNN_fr_cd = cat(3,rp_vec_nd{:});

%% find aberrant trials
mean_left_cd = zeros(endpoint-500,1);
mean_right_cd = zeros(endpoint-500,1);
std_left_cd = zeros(endpoint-500,1);
std_right_cd = zeros(endpoint-500,1);

for k = 1:endpoint-1500
    kk = k+1500;
    mean_left_cd(k) = mean(vecnorm(squeeze(RNN_fr_cd(:,kk, 1:N_trials_cd/2))));
    mean_right_cd(k) = mean(vecnorm(squeeze(RNN_fr_cd(:,kk, N_trials_cd/2+1:end))));
   
    std_left_cd(k) = std(vecnorm(squeeze(RNN_fr_cd(:,kk, 1:N_trials_cd/2))));
    std_right_cd(k) = std(vecnorm(squeeze(RNN_fr_cd(:,kk, N_trials_cd/2+1:end))));    

    aberrant_right_cd_temp{k} = find((vecnorm(squeeze(RNN_fr_cd(:,kk,right_trials_cd))) >...
    (mean_right_cd(k) + 6*std_right_cd(k))) & (vecnorm(squeeze(RNN_fr_cd(:,kk,right_trials_cd))) >...
    (mean_left_cd(k) + 6*std_left_cd(k))));

    aberrant_left_cd_temp{k} = find((vecnorm(squeeze(RNN_fr_cd(:,kk,left_trials_cd))) >...
    (mean_right_cd(k) + 6*std_right_cd(k))) & (vecnorm(squeeze(RNN_fr_cd(:,kk,left_trials_cd))) >...
    (mean_left_cd(k) + 6*std_left_cd(k))));
end

%%
[aberrant_right_cd,~] = unique([aberrant_right_cd_temp{:}]);

[aberrant_left_cd,~] = unique([aberrant_left_cd_temp{:}]);

%% find correct and error trials

correct_tri_left_cd = setdiff(find(mean(RNN_fr_cd(:,endpoint,left_trials_cd),1)<mean(mean(RNN_fr_cd(:,endpoint,:),3),1)), aberrant_left_cd);
correct_tri_right_cd = setdiff(find(mean(RNN_fr_cd(:,endpoint,right_trials_cd),1)>mean(mean(RNN_fr_cd(:,endpoint,:),3),1)), aberrant_right_cd);

error_tri_left_cd = setdiff(find(mean(RNN_fr_cd(:,endpoint,left_trials_cd),1)>mean(mean(RNN_fr_cd(:,endpoint,:),3),1)), aberrant_left_cd);
error_tri_right_cd = setdiff(find(mean(RNN_fr_cd(:,endpoint,right_trials_cd),1)<mean(mean(RNN_fr_cd(:,endpoint,:),3),1)), aberrant_right_cd);

RNN_fr_left_mean_cd = mean(RNN_fr_cd(:,:,correct_tri_left_cd),3);
RNN_fr_right_mean_cd = mean(RNN_fr_cd(:,:,N_trials_cd/2 + correct_tri_right_cd),3);

%% compute CD mode from Data
proj_cd_data = fr_right_norm - fr_left_norm;
cd_late_delay_data = mean(proj_cd_data(:,end-cd_span/5:end),2);

%% compute CD mode from RNN
proj_cd = RNN_fr_right_mean_cd - RNN_fr_left_mean_cd;
cd_late_delay = mean(proj_cd(:,3500-cd_span:3500),2);

%% compute Sample mode from RNN
cd_sample = mean(proj_cd(:,1001:1400),2);

%% Angle between CD Data and CD RNN

angle = (cd_late_delay_data./norm(cd_late_delay_data))'*cd_late_delay./norm(cd_late_delay);
disp(['Coding directions - Data vs RNN = ', num2str(angle)])

%% save output struct

s_out.cd_sample = cd_sample;
s_out.cd_late_delay = cd_late_delay;
s_out.proj_cd = proj_cd;
s_out.proj_cd_data = proj_cd_data;
s_out.RNN_fr_cd = single(RNN_fr_cd);
s_out.correct_tri_left_cd = correct_tri_left_cd;
s_out.correct_tri_right_cd = correct_tri_right_cd;
s_out.error_tri_left_cd = error_tri_left_cd;
s_out.error_tri_right_cd = error_tri_right_cd;
s_out.aberrant_left_cd = aberrant_left_cd;
s_out.aberrant_right_cd = aberrant_right_cd;
s_out.endpoint = endpoint;
s_out.N_trials_cd = N_trials_cd;
s_out.mean_left_cd = mean_left_cd;
s_out.mean_right_cd = mean_right_cd;
s_out.std_left_cd = std_left_cd;
s_out.std_right_cd = std_right_cd;

