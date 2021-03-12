function [ struct_out ] = f_cd_ALM(N_trials_cd, ramp, ramp_prefactor, ramp_dur, ramp_sigma, noise_sigma, stim_shape_in, p, T_test, cd_span)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Unpack Params struct
dt = p.dt;
N = p.N;
tau = p.tau;
ramp_train = p.ramp_train;
f0 = p.f0;
beta0 = p.beta0;
theta0 = p.theta0;
W = p.W;
eff_dt = p.eff_dt;
% simtime_len = p.simtime_len;
b = p.b;
des_out_left = p.des_out_left;
des_out_right = p.des_out_right;
des_r_right = p.des_r_right;
des_r_left = p.des_r_left;
des_r_left_norm = p.des_r_left_norm;
des_r_right_norm = p.des_r_right_norm;
des_r_left_norm_cell = p.des_r_left_norm_cell;
des_r_right_norm_cell = p.des_r_right_norm_cell;
simtime = p.simtime;
idx_sorted = p.idx_sorted;
idx_sorted_l = p.idx_sorted_l;
stim_amp = p.stim_amp;
chirp_amp = p.chirp_amp;
fr_smooth = p.fr_smooth;
ramp_bsln = p.ramp_bsln;

noise_sigma_eff = sqrt(dt).*noise_sigma./tau; %% effective sigma noise

% time vectors
simtime_test = [0:dt:T_test-dt];
simtime_test_len = length(simtime_test);

% t_vec_coarse = [1:5/dt:T_test/dt];

%% Chirps time series

inp_chirp_temp = zeros(simtime_test_len,1);
inp_chirp_temp([500/dt+1:650/dt, 1350/dt+1:1500/dt]) = 1;
inp_chirp_temp = smooth(inp_chirp_temp, fr_smooth);

%% type of stimulus (square pulse or sine wave)
switch stim_shape_in
    case 'sinewave'
        stim_shape = abs(sin(0.005*2*pi*[1:400]./dt));
    case 'square'
        stim_shape = 1;
end

%% test without distractor

disp(['Running a few trials to compute the CD mode']);

% inp_ext = zeros(N,simtime_test_len);
% inp_rec = zeros(N,simtime_test_len);

rp_vec_nd = [];

r_in_cd(1,:) = chirp_amp*[inp_chirp_temp];
t_ramp_start = 500/dt+1;

for i=1:N_trials_cd
    i
    
    % determine ramping input
    switch ramp
        case 'endpoint'
            
            inp_ramp_test = zeros(T_test,1);
            inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./ramp_dur;
            %             inp_ramp_test = smooth(inp_ramp_test, fr_smooth);
            
            %             inp_ramp_test = zeros(T,1);
            %             inp_ramp_test(501:500 + 1100) = [1:1100]./1100;
            
            r_in_cd(3,1:t_ramp_start+ramp_dur-1) = ramp_prefactor.*ramp_train.*inp_ramp_test(1:t_ramp_start+ramp_dur-1).*...
                (1+ramp_sigma*randn) + ramp_bsln;
            r_in_cd(3,t_ramp_start+ramp_dur:end) = r_in_cd(3,t_ramp_start+ramp_dur-1);
            %         r_in_cd(3,2301:end) = r_in_cd(3,2300);
            
            r_in_cd(3,:) = smooth(r_in_cd(3,:), fr_smooth);
            
        case 'delay'
            
            inp_ramp_test = zeros(T_test,1);
            inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./(3000/dt);
            %             inp_ramp_test = smooth(inp_ramp_test, fr_smooth);
            
            %             inp_ramp_test = zeros(T,1);
            %             inp_ramp_test(501:500 + 1100) = [1:1100]./1100;
            
            r_in_cd(3,1:t_ramp_start+ramp_dur-1) = ramp_prefactor.*ramp_train.*inp_ramp_test(1:t_ramp_start+ramp_dur-1).*...
                (1+ramp_sigma*randn) + ramp_bsln;
            r_in_cd(3,t_ramp_start+ramp_dur:end) = r_in_cd(3,t_ramp_start+ramp_dur-1);
            %         r_in_cd(3,2301:end) = r_in_cd(3,2300);
            
            r_in_cd(3,:) = smooth(r_in_cd(3,:), fr_smooth);
            
        case 'no_ramp'
            r_in_cd(3,:) = ramp_bsln;
            
        case 'descending_endpoint'
            inp_ramp_test = zeros(T_test,1);
            inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./ramp_dur;
            %             inp_ramp_test = smooth(inp_ramp_test, fr_smooth);
            
            %             inp_ramp_test = zeros(T,1);
            %             inp_ramp_test(501:500 + 1100) = [1:1100]./1100;
            
            r_in_cd(3,:) = ramp_prefactor.*ramp_train.*inp_ramp_test + ramp_bsln;
            r_in_cd(3,t_ramp_start+ramp_dur:end) = r_in_cd(3,t_ramp_start+ramp_dur-1);
            
            r_in_cd(3,3000:3500) =  [500:-1:0]./500  + ramp_bsln;
            r_in_cd(3,3501:end) = r_in_cd(3,3500);
            
            r_in_cd(3,:) = smooth(r_in_cd(3,:), fr_smooth);
    end
    
    rp_vec_nd{i} = zeros(N,simtime_test_len);
    xp_vec_nd{i} = zeros(N,simtime_test_len);
    
    
    % determine stimulus instruction
    if i<=N_trials_cd/2      % first half of all trials are lick left (no stimulus)
        r_in_cd(2,:) = 0;
        des_out = des_out_left;
    else                     % second half of all trials are lick right (stimulus)
        inp_stim_temp = zeros(simtime_test_len,1);
        t_stim_interval = [1000/dt+1:1400/dt];
        
        inp_stim_temp(t_stim_interval) = stim_amp*stim_shape;
        inp_stim_temp = smooth(inp_stim_temp, fr_smooth);
        
        r_in_cd(2,:) = inp_stim_temp;
        des_out = des_out_right;

    end
    
    % initial conditions
%    xp = mean(des_out(:,1:10),2);       
    xp = mean([mean(des_out_left(:,1:20),2), mean(des_out_right(:,1:20),2)],2);       

    xp(N+1:N+3) = 0;  
    rp = f0./(1.0 + exp(-beta0.*(xp(1:N)-theta0)));
    
    % dynamics
    ti = 0;
    for t = simtime_test
        
        ti = ti+1;

        xp = xp + (-xp + W*[rp;r_in_cd(:,ti)] + b)*eff_dt + [noise_sigma_eff*randn(N,1);zeros(3,1)];

        rp = f0./(1.0+exp(-beta0.*(xp(1:N)-theta0)));
        
        rp_vec_nd{i}(:,ti) = rp;

        xp_vec_nd{i}(:,ti) = xp(1:N);

        % save inputs for inspection
%         inp_rec(1:N,ti) = W(1:N,1:N)*rp;
%         inp_ext(1:N,ti) = W(1:N,N+1:N+3)*r_in_cd(:,ti) + b(1:N);
        
    end
    
end


%% divide trials according to instruction
left_trials_nd = 1:N_trials_cd/2;

right_trials_nd = N_trials_cd/2 + 1:N_trials_cd;

rp_nd_mat_all_cd = cat(3,rp_vec_nd{:});

xp_nd_mat_all_cd = cat(3,xp_vec_nd{:});


%% Plot averages across trials

figure
hold on
for i = 1:N_trials_cd
    
    if i<=N_trials_cd/2
        plot(simtime_test,mean(rp_nd_mat_all_cd(:,:,i)),'r');
        
    else
        plot(simtime_test,mean(rp_nd_mat_all_cd(:,:,i)),'b');
    end
end

plot(simtime_test(1:5:3300), mean(des_r_left),'m')
plot(simtime_test(1:5:3300), mean(des_r_right),'c')

ylabel('Spike Rate (Hz)')
xlabel('Time to Go cue (s)')
title('Mean Firing rate (Network)')

%% find aberrant trials

endpoint = 3500/dt;

% aberrant_right_nd = find(mean(rp_nd_mat_all_cd(:,endpoint,right_trials_nd),1) - mean(mean(rp_nd_mat_all_cd(:,endpoint,right_trials_nd),3),1) >...
%     mean(mean(rp_nd_mat_all_cd(:,endpoint,right_trials_nd),3),1) + 6*std(mean(rp_nd_mat_all_cd(:,endpoint,right_trials_nd),1),[],3));
%
% aberrant_left_nd = find(mean(rp_nd_mat_all_cd(:,endpoint,left_trials_nd),1) - mean(mean(rp_nd_mat_all_cd(:,endpoint,left_trials_nd),3),1) <...
%     mean(mean(rp_nd_mat_all_cd(:,endpoint,left_trials_nd),3),1) - 6*std(mean(rp_nd_mat_all_cd(:,endpoint,left_trials_nd),1),[],3));

aberrant_right_nd = find(mean(rp_nd_mat_all_cd(:,endpoint,right_trials_nd),1) >...
    (mean(mean(rp_nd_mat_all_cd(:,endpoint,right_trials_nd),3),1) + 6*std(mean(rp_nd_mat_all_cd(:,endpoint,right_trials_nd),1),[],3)));

aberrant_left_nd = find(mean(rp_nd_mat_all_cd(:,endpoint,left_trials_nd),1) <...
    (mean(mean(rp_nd_mat_all_cd(:,endpoint,left_trials_nd),3),1) - 6*std(mean(rp_nd_mat_all_cd(:,endpoint,left_trials_nd),1),[],3)));

%% find correct and error trials

var_vec = mean(squeeze(var(rp_nd_mat_all_cd(:,end-1300/dt:end,:),0,2)));
aberrant_var = find(var_vec>0.01);

correct_tri_left_nd = setdiff(find(mean(rp_nd_mat_all_cd(:,endpoint,left_trials_nd),1)<mean(mean(rp_nd_mat_all_cd(:,endpoint,:),3),1)), aberrant_var);
correct_tri_right_nd = setdiff(find(mean(rp_nd_mat_all_cd(:,endpoint,right_trials_nd),1)>mean(mean(rp_nd_mat_all_cd(:,endpoint,:),3),1)), aberrant_var);

error_tri_left_nd = setdiff(find(mean(rp_nd_mat_all_cd(:,endpoint,left_trials_nd),1)>mean(mean(rp_nd_mat_all_cd(:,endpoint,:),3),1)), aberrant_var);
error_tri_right_nd = setdiff(find(mean(rp_nd_mat_all_cd(:,endpoint,right_trials_nd),1)<mean(mean(rp_nd_mat_all_cd(:,endpoint,:),3),1)), aberrant_var);

rp_nd_left_mean = mean(rp_nd_mat_all_cd(:,:,correct_tri_left_nd),3);

rp_nd_right_mean = mean(rp_nd_mat_all_cd(:,:,N_trials_cd/2 + correct_tri_right_nd),3);

%% compute CD mode Data
% proj_vector = computeModeWeights(psth_t_u_tr, trials1, trials2, tint1, tint2, psth_t_vector, mintrials_modeweights)

proj_cd_data = des_r_right_norm - des_r_left_norm;

rate_f1 = proj_cd_data;
rate_f2 = proj_cd_data;

[rho_cd_data,~] = corr(rate_f1,rate_f2);

%% Plot cross-correlation (data)
figure
hold on
uimagesc(simtime,simtime,single(rho_cd_data))
plot(1000*ones(2,1)/dt,[0,3500]/dt,'k--','Linewidth',2)
plot(1400*ones(2,1)/dt,[0,3500]/dt,'k--','Linewidth',2)
plot([0,3500]/dt,1000*ones(2,1)/dt,'k--','Linewidth',2)
plot([0,3500]/dt,1400*ones(2,1)/dt,'k--','Linewidth',2)
plot(1500*ones(2,1)/dt,[0,3500]/dt,'r--','Linewidth',2)
plot([0,3500]/dt,1500*ones(2,1)/dt,'r--','Linewidth',2)
plot(500*ones(2,1)/dt,[0,3500]/dt,'r--','Linewidth',2)
plot([0,3500]/dt,500*ones(2,1)/dt,'r--','Linewidth',2)
xlim([-0.5,3500])
ylim([-0.5,3500])
ax = gca;
ax.CLim = [0, 1];
colormap(jet)
colorbar
xlabel('Time (ms)')
ylabel('Time (ms)')
title('cd vector (data) autocorrs')
% print(['Cross_corr', num2str(k), ' vs ', num2str(l)], '-djpeg'); %<-Save as EPS with 300 DPI

%% CD mode RNN
proj_cd = rp_nd_right_mean(:,1:1/dt:end) - rp_nd_left_mean(:,1:1/dt:end);

% consider only negative weights
proj_cd_neg = proj_cd;
proj_cd_neg(proj_cd>0) = 0;

%%
rate_f1 = proj_cd;
rate_f2 = proj_cd;

[rho_cd,~] = corr(rate_f1,rate_f2);

%% Plot cross-correlation (model)
figure
hold on
uimagesc(simtime_test(1:1/dt:end),simtime_test(1:1/dt:end),single(rho_cd))
plot(1000*ones(2,1),[0,3500],'k--','Linewidth',2)
plot(1400*ones(2,1),[0,3500],'k--','Linewidth',2)
plot([0,3500],1000*ones(2,1),'k--','Linewidth',2)
plot([0,3500],1400*ones(2,1),'k--','Linewidth',2)
plot(1500*ones(2,1),[0,3500],'r--','Linewidth',2)
plot([0,3500],1500*ones(2,1),'r--','Linewidth',2)
plot(500*ones(2,1),[0,3500],'r--','Linewidth',2)
plot([0,3500],500*ones(2,1),'r--','Linewidth',2)
xlim([-0.5,3500])
ylim([-0.5,3500])
ax = gca;
ax.CLim = [0, 1];
% xlim([0, simtime_test(3400)])
% ylim([0, simtime_test(3400)])
colormap(jet)
colorbar
xlabel('Time (ms)')
ylabel('Time (ms)')
title('cd vector (model) autocorrs')
%
%% compare CD mode data vs RNN

rate_f1 = proj_cd_data;
rate_f2 = proj_cd(:,1:5:3300);

[rho_cd_vs,~] = corr(rate_f1,rate_f2);

%% Plot cross-correlation Data vs Model
figure
hold on
uimagesc(simtime,simtime,single(rho_cd_vs))
plot(1000*ones(2,1),[0,3500],'k--','Linewidth',2)
plot(1400*ones(2,1),[0,3500],'k--','Linewidth',2)
plot([0,3500],1000*ones(2,1),'k--','Linewidth',2)
plot([0,3500],1400*ones(2,1),'k--','Linewidth',2)
plot(1500*ones(2,1),[0,3500],'r--','Linewidth',2)
plot([0,3500],1500*ones(2,1),'r--','Linewidth',2)
plot(500*ones(2,1),[0,3500],'r--','Linewidth',2)
plot([0,3500],500*ones(2,1),'r--','Linewidth',2)
xlim([-0.5,3500])
ylim([-0.5,3500])
ax = gca;
ax.CLim = [0, 1];
colormap(jet)
colorbar
xlabel('Time (ms)')
ylabel('Time (ms)')
title('cd vector RNN vs Data cross-corrs')

%% projections onto CD

cd_late_delay = mean(proj_cd(:,3500-cd_span:3500),2);

cd_late_delay_neg = mean(proj_cd_neg(:,3500-cd_span:3500),2);


cd_sample = mean(proj_cd(:,1001:1400),2);

%% Compute orthogonal basis
orth_basis = null(cd_late_delay'./norm(cd_late_delay));

[~,max_vec_idx] = max(abs(orth_basis'*(cd_sample./norm(cd_sample))));
cd_sample_orth = orth_basis(:,max_vec_idx);

%% save output struct

struct_out.cd_sample = cd_sample;
struct_out.cd_sample_orth = cd_sample_orth;

struct_out.cd_late_delay = cd_late_delay;
struct_out.cd_late_delay_neg = cd_late_delay_neg;
struct_out.proj_cd = proj_cd;
struct_out.proj_cd_neg = proj_cd_neg;
struct_out.proj_cd_data = proj_cd_data;

struct_out.rp_nd_mat_all_cd = rp_nd_mat_all_cd;
struct_out.xp_nd_mat_all_cd = xp_nd_mat_all_cd;


struct_out.correct_tri_left_nd = correct_tri_left_nd;
struct_out.correct_tri_right_nd = correct_tri_right_nd;
struct_out.error_tri_left_nd = error_tri_left_nd;
struct_out.error_tri_right_nd = error_tri_right_nd;
struct_out.aberrant_left_nd = aberrant_left_nd;
struct_out.aberrant_right_nd = aberrant_right_nd;
struct_out.rho_cd_vs = rho_cd_vs;
struct_out.endpoint = endpoint;
struct_out.N_trials_cd = N_trials_cd;

struct_out.r_in_cd = r_in_cd;

%% Plot spike rates vs data

figure
hold on
uimagesc(simtime_test(200/dt:3400/dt)/1000-1.5,1:N,rp_nd_mat_all_cd(idx_sorted_l,200/dt:3400/dt,1)...
    ./repmat(max(rp_nd_mat_all_cd(idx_sorted_l,200/dt:3400/dt,1),[],2),[1,size(rp_nd_mat_all_cd(idx_sorted_l,200/dt:3400/dt,1),2)]))
plot(zeros(2,1),[1,N-1],'k--','Linewidth',2)
plot(-1.*ones(2,1),[1,N-1],'k--','Linewidth',2)
ylim([0,N+1])
xlim([-1.4,1.8])
ax = gca;
ax.CLim = [0 1];
xlabel('Time to Go cue (s)')
ylabel('Neuron No.')
title('Lick left RNN')


figure
hold on
uimagesc(simtime_test(200/dt+1:5/dt:3400/dt)/1000-1.5,1:N,des_r_left_norm_cell(idx_sorted_l,200/5:end)...
    ./repmat(max(des_r_left_norm_cell(idx_sorted_l,200/5:end),[],2),[1,size(des_r_right_norm_cell(:,200/5:end),2)]))
plot(zeros(2,1),[1,N-1],'k--','Linewidth',2)
plot(-1.*ones(2,1),[1,N-1],'k--','Linewidth',2)
ylim([0,N+1])
xlim([-1.5,1.8])
ax = gca;
ax.CLim = [0 1];
xlabel('Time to Go cue (s)')
ylabel('Neuron No.')
title('Lick left Data')


figure
hold on
uimagesc(simtime_test(200/dt:3400/dt)/1000-1.5,1:N,rp_nd_mat_all_cd(idx_sorted,200/dt:3400/dt, N_trials_cd/2 + 1)...
    ./repmat(max(rp_nd_mat_all_cd(idx_sorted,200/dt:3400/dt, N_trials_cd/2 + 1),[],2),[1,size(rp_nd_mat_all_cd(idx_sorted,200/dt:3400/dt,3),2)]))
plot(zeros(2,1),[1,N-1],'k--','Linewidth',2)
plot(-1.*ones(2,1),[1,N-1],'k--','Linewidth',2)
ylim([0,N+1])
xlim([-1.5,1.8])
xlabel('Time to Go cue (s)')
ylabel('Neuron No.')
title('Lick right RNN')


figure
hold on
uimagesc(simtime_test(200/dt+1:5/dt:3400/dt)/1000-1.5,1:N,des_r_right_norm_cell(idx_sorted,200/5:end)...
    ./repmat(max(des_r_right_norm_cell(idx_sorted,200/5:end),[],2),[1,size(des_r_right_norm_cell(idx_sorted,200/5:end),2)]))
plot(zeros(2,1),[1,N-1],'k--','Linewidth',2)
plot(-1.*ones(2,1),[1,N-1],'k--','Linewidth',2)
ylim([0,N+1])
xlim([-1.5,1.8])
xlabel('Time to Go cue (s)')
ylabel('Neuron No.')
title('Lick right Data')


