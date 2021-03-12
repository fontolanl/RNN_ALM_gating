function [ struct_out ] = f_unperturbed_ALM(N_trials_unperturbed, ramp, ramp_prefactor, ramp_dur, sigma_ramp, p, T_test, s_cd,...
    sigma_stim, sigma_noise, plot_figs,  stim_shape_in)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%% Params
dt = p.dt;
N = p.N;
% inp_chirp_temp = p.inp_chirp_temp;
% inp_stim_temp = p.inp_stim_temp;
tau = p.tau;
ramp_train = p.ramp_train;
% tau_noise = p.tau_noise;
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
simtime = p.simtime;
idx_sorted = p.idx_sorted;
idx_sorted_l = p.idx_sorted_l;
fr_smooth = p.fr_smooth;
ramp_bsln = p.ramp_bsln;
chirp_amp = p.chirp_amp;
stim_amp = p.stim_amp;

% sigma_stim = p.sigma_stim;
% sigma_stim/2; % changed Feb 22 2019
T = p.T;

%% effective sigma noise
sigma_noise_eff = sqrt(dt).*sigma_noise./tau;

%% from cd mode
cd_late_delay = s_cd.cd_late_delay;
cd_late_delay_neg = s_cd.cd_late_delay_neg;
cd_sample = s_cd.cd_sample;
cd_sample_orth = s_cd.cd_sample_orth;
rp_nd_mat_all_cd = s_cd.rp_nd_mat_all_cd;
% endpoint = s_cd.endpoint;
endpoint = 3500/dt;
N_trials_cd = s_cd.N_trials_cd;

%%
switch stim_shape_in
    case 'sinewave'
        stim_shape = abs(sin(0.005*2*pi*[1:400]));
    case 'square'
        stim_shape = 1;      
end

%% Run trials without distractors
disp(['Now testing...without distractors.']);

simtime_test = [0:dt:T_test-dt];
simtime_test_len = length(simtime_test);

%% Chirp time series

inp_chirp_temp = zeros(simtime_test_len,1);
inp_chirp_temp([500/dt+1:650/dt, 1350/dt+1:1500/dt]) = 1;
inp_chirp_temp = smooth(inp_chirp_temp, fr_smooth);

r_in = zeros(3,simtime_test_len);
r_in(1,:) = chirp_amp.*[inp_chirp_temp];

%% test unperturbed trials
inp_ext = zeros(3,N,simtime_test_len);
inp_rec = zeros(N,simtime_test_len);
inp_noise = zeros(N,simtime_test_len);

rp_vec_nd = [];

simtime_test_len_coarse = length(simtime_test(1:1/dt:end));

for i=1:N_trials_unperturbed
    i
    rp_vec_nd{i} = zeros(N,simtime_test_len_coarse);
    
    xp = mean([mean(des_out_left(:,1:20),2), mean(des_out_right(:,1:20),2)],2) + sigma_noise*randn(N,1);       
    xp(N+1:N+3) = 0;
    rp = f0./(1.0+exp(-beta0.*(xp(1:N)-theta0)));
        
    if i<=N_trials_unperturbed/2
        r_in(2,:) = 0;
    else        
        t_stim_interval = [1000/dt+1:1400/dt] + randi([-10 10],1)/dt;
        
        inp_stim_temp = zeros(simtime_test_len,1);
        inp_stim_temp(t_stim_interval) = stim_amp.*stim_shape;
        inp_stim_temp = smooth(inp_stim_temp, fr_smooth);
        
        inp_stim_test = [inp_stim_temp];
                
        r_in(2,:) = inp_stim_test.*(1+sigma_stim*randn);
    end
    
    t_ramp_start = 500/dt + 1 + randi([-10 10],1)/dt;

    
    switch ramp

        case 'clamped'
            
            ramp_dur = ramp_prefactor.*3000/dt;
            inp_ramp_test = zeros(T,1);
            inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./ramp_dur;
     
            r_in(3,1:t_ramp_start+ramp_dur-1) = ramp_prefactor.*ramp_train.*inp_ramp_test(1:t_ramp_start+ramp_dur-1).*...
                (1+sigma_ramp*randn) + ramp_bsln;
            
            r_in(3,t_ramp_start+ramp_dur:end) = r_in(3,t_ramp_start+ramp_dur-1);

            r_in(3,:) = smooth(r_in(3,:), fr_smooth);
        
        case 'endpoint'
            
            inp_ramp_test = zeros(T,1);
            inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./ramp_dur;


            
            r_in(3,1:t_ramp_start+ramp_dur-1) = ramp_prefactor.*ramp_train.*inp_ramp_test(1:t_ramp_start+ramp_dur-1).*...
                (1+sigma_ramp*randn) + ramp_bsln;
            
            r_in(3,t_ramp_start+ramp_dur:end) = r_in(3,t_ramp_start+ramp_dur-1);

            r_in(3,:) = smooth(r_in(3,:), fr_smooth);
            
        case 'delay'
            inp_ramp_test = zeros(T,1);
            inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./(3000/dt);
            %             inp_ramp_test = smooth(inp_ramp_test, fr_smooth);
            
            %             inp_ramp_test = zeros(T,1);
            %             inp_ramp_test(501:500 + 1100) = [1:1100]./1100;
            
            r_in(3,1:t_ramp_start+ramp_dur-1) = ramp_prefactor.*ramp_train.*inp_ramp_test(1:t_ramp_start+ramp_dur-1).*...
                (1+sigma_ramp*randn) + ramp_bsln;
            r_in(3,t_ramp_start+ramp_dur:end) = r_in(3,t_ramp_start+ramp_dur-1);
            %         r_in_cd(3,2301:end) = r_in_cd(3,2300);
            
            r_in(3,:) = smooth(r_in(3,:), fr_smooth);
   
        case 'flat'
            
            inp_ramp_test = zeros(T,1);
%             inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./ramp_dur;
            inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./3000;

            %             inp_ramp_test = smooth(inp_ramp_test, fr_smooth);
            
            %             inp_ramp_test = zeros(T,1);
            %             inp_ramp_test(501:500 + 1100) = [1:1100]./1100;
            
            r_in(3,1:t_ramp_start+ramp_dur-1) = ramp_prefactor.*ramp_train.*inp_ramp_test(1:t_ramp_start+ramp_dur-1).*...
                (1+sigma_ramp*randn) + ramp_bsln;
%             r_in(3,1:t_ramp_start + ramp_dur -1) = ramp_prefactor.*ramp_train.*inp_ramp_test.*(1+sigma_ramp*randn) + ramp_bsln;
            
            r_in(3,t_ramp_start+ramp_dur:end) = r_in(3,t_ramp_start+ramp_dur-1);
%             r_in(3,2301:end) = r_in(3,2300);

            r_in(3,:) = smooth(r_in(3,:), fr_smooth);
            
        case 'no_ramp'
            r_in(3,:) = ramp_bsln;
            
        case 'descending_endpoint'
            inp_ramp_test = zeros(T,1);
            inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./ramp_dur;
            %             inp_ramp_test = smooth(inp_ramp_test, fr_smooth);
            
            %             inp_ramp_test = zeros(T,1);
            %             inp_ramp_test(501:500 + 1100) = [1:1100]./1100;
            
            r_in(3,:) = ramp_prefactor.*ramp_train.*inp_ramp_test.*(1+sigma_ramp*randn) + ramp_bsln;
            r_in(3,t_ramp_start+ramp_dur:end) = r_in(3,t_ramp_start+ramp_dur-1);
            
            r_in(3,3000:3500) =  [500:-1:0]./500  + ramp_bsln;
            r_in(3,3501:end) = r_in(3,3500);
            
            r_in(3,:) = smooth(r_in(3,:), fr_smooth);
    end
    %     figure
    ti = 0;
    tii = 0;
    for t = simtime_test
        
        ti = ti+1;
        
        xp = xp + (-xp + W*[rp; r_in(:,ti)] + b)*eff_dt + [sigma_noise_eff*randn(N,1);zeros(3,1)];      
        
        rp = f0./(1.0+exp(-beta0.*(xp(1:N)-theta0)));
        
        if mod(ti,1/dt)==0
            tii = tii + 1;
            rp_vec_nd{i}(:,tii) = single(rp);
        end
        
        if i ==1
            inp_rec(1:N,ti) = W(1:N,1:N)*rp;
            inp_ext(1,1:N,ti) = W(1:N,N+1)*r_in(1,ti);
            inp_ext(2,1:N,ti) = W(1:N,N+2)*r_in(2,ti);
            inp_ext(3,1:N,ti) = W(1:N,N+3)*r_in(3,ti);
        end  
        
    end
  
end

%% divide trials
left_trials_nd = 1:N_trials_unperturbed/2;

right_trials_nd = (N_trials_unperturbed/2 + 1):N_trials_unperturbed;

rp_nd_mat_all = cat(3,rp_vec_nd{:});

clear rp_vec_nd

%% All trials projected no distractor

rp_nd_proj = [];

rp_nd_proj_neg = [];

for i = 1:N_trials_unperturbed
    rp_nd_proj(i,:) = squeeze(rp_nd_mat_all(:,:,i))'*cd_late_delay;
    rp_nd_proj_neg(i,:) = squeeze(rp_nd_mat_all(:,:,i))'*cd_late_delay_neg;
    
end


%% find aberrant trials

mean_left_cd = mean(rp_nd_mat_all_cd(:,endpoint,1));
mean_right_cd = mean(rp_nd_mat_all_cd(:,endpoint,N_trials_cd/2 + 1));
std_left_cd = 0.01;
std_right_cd = 0.01;


aberrant_right_nd = find((mean(rp_nd_mat_all(:,endpoint,right_trials_nd),1) >...
    (mean_right_cd + 4*std_right_cd)) | (mean(rp_nd_mat_all(:,endpoint,right_trials_nd),1) <...
    (mean_left_cd - 4*std_left_cd)));

aberrant_left_nd = find((mean(rp_nd_mat_all(:,endpoint,left_trials_nd),1) >...
    (mean_right_cd + 4*std_right_cd)) | (mean(rp_nd_mat_all(:,endpoint,left_trials_nd),1) <...
    (mean_left_cd - 4*std_left_cd)));

%% find correct trials

proj_mean = mean([rp_nd_mat_all_cd(:,1:1/dt:end,1)'*cd_late_delay, rp_nd_mat_all_cd(:,1:1/dt:end,N_trials_cd/2 + 1)'*cd_late_delay],2);

correct_tri_left_nd = setdiff(find(rp_nd_proj(left_trials_nd,endpoint)<proj_mean(endpoint)),aberrant_left_nd);
correct_tri_right_nd = setdiff(find(rp_nd_proj(right_trials_nd,endpoint)>proj_mean(endpoint)),aberrant_right_nd);

error_tri_left_nd = setdiff(find(rp_nd_proj(left_trials_nd,endpoint)>proj_mean(endpoint)),aberrant_left_nd);
error_tri_right_nd = setdiff(find(rp_nd_proj(right_trials_nd,endpoint)<proj_mean(endpoint)),aberrant_right_nd);

rp_nd_left_mean = mean(rp_nd_mat_all(:,:,correct_tri_left_nd),3);

rp_nd_right_mean = mean(rp_nd_mat_all(:,:,N_trials_unperturbed/2 + correct_tri_right_nd),3);


%% Correct and error trials projected
if plot_figs
    t_vec_plot = (simtime_test - 3500/dt)./1000;
    figure
    hold on
    for i = 1:N_trials_unperturbed
        
        if ismember(i,N_trials_unperturbed/2 + correct_tri_right_nd)
            plot(t_vec_plot, rp_nd_proj(i,:),'b');
        elseif ismember(i,correct_tri_left_nd)
            plot(t_vec_plot, rp_nd_proj(i,:),'r');
        end
        
        if ismember(i,N_trials_unperturbed/2 + error_tri_right_nd)
            plot(t_vec_plot, rp_nd_proj(i,:),'m');
        elseif ismember(i,error_tri_left_nd)
            plot(t_vec_plot, rp_nd_proj(i,:),'c');
        end
        
    end
        
    plot(t_vec_plot, proj_mean, 'Color', [128, 128, 128]./255,'Linewidth',2)
    plot(t_vec_plot, mean(rp_nd_proj(N_trials_unperturbed/2 + correct_tri_right_nd,:)),'c','Linewidth',2);
    plot(t_vec_plot, mean(rp_nd_proj(correct_tri_left_nd,:)),'m','Linewidth',2);
    plot(t_vec_plot, mean(rp_nd_mat_all_cd(:,1:1/dt:end,1:N_trials_cd/2),3)'*cd_late_delay,'m--','Linewidth',3)
    plot(t_vec_plot, mean(rp_nd_mat_all_cd(:,1:1/dt:end,N_trials_cd/2+1:end),3)'*cd_late_delay,'c--','Linewidth',3)
    plot(zeros(2,1),[0,max(rp_nd_proj(:))+0.1],'k--','Linewidth',2')
    plot((500-3500)*ones(2,1)/1000./dt,[0,max(rp_nd_proj(:))+0.1],'k--','Linewidth',2')
    plot((1000-3500)*ones(2,1)/1000./dt,[0,max(rp_nd_proj(:))+0.1],'k--','Linewidth',2')
    plot((1500-3500)*ones(2,1)/1000./dt,[0,max(rp_nd_proj(:))+0.1],'k--','Linewidth',2')
    ylabel('Spike rate (Hz)')
    xlabel('Time to Go cue (s)')
    title('Projection onto Choice mode (No aberrant trials)')
    
    
    figure
    hold on
    for i = 1:N_trials_unperturbed/2       
            plot(t_vec_plot, rp_nd_proj(i,:),'r');        
    end
    for i = N_trials_unperturbed/2 + 1:N_trials_unperturbed       
        plot(t_vec_plot, rp_nd_proj(i,:),'b');
    end
    
    plot(t_vec_plot, proj_mean, 'Color', [128, 128, 128]./255,'Linewidth',2)
    plot(t_vec_plot, mean(rp_nd_proj(N_trials_unperturbed/2 + correct_tri_right_nd,:)),'c','Linewidth',2);
    plot(t_vec_plot, mean(rp_nd_proj(correct_tri_left_nd,:)),'m','Linewidth',2);
    plot(t_vec_plot, mean(rp_nd_mat_all_cd(:,1:1/dt:end,1:N_trials_cd/2),3)'*cd_late_delay,'m--','Linewidth',3)
    plot(t_vec_plot, mean(rp_nd_mat_all_cd(:,1:1/dt:end,N_trials_cd/2+1:end),3)'*cd_late_delay,'c--','Linewidth',3)
    plot(zeros(2,1),[0,max(rp_nd_proj(:))+0.1],'k--','Linewidth',2')
    plot((500-3500)*ones(2,1)/1000./dt,[0,max(rp_nd_proj(:))+0.1],'k--','Linewidth',2')
    plot((1000-3500)*ones(2,1)/1000./dt,[0,max(rp_nd_proj(:))+0.1],'k--','Linewidth',2')
    plot((1500-3500)*ones(2,1)/1000./dt,[0,max(rp_nd_proj(:))+0.1],'k--','Linewidth',2')
    ylabel('Spike rate (Hz)')
    xlabel('Time to Go cue (s)')
    title('Projection onto Choice mode (All Trials)')
    
end

%% All trials projected no distractor - Sample mode orthogonalized
rp_nd_proj_s = [];


for i = 1:N_trials_unperturbed
    rp_nd_proj_s(i,:) = squeeze(rp_nd_mat_all(:,:,i))'*cd_sample_orth;

end


%% save output space

struct_out.rp_nd_mat_all = single(rp_nd_mat_all);

struct_out.correct_tri_left_nd = correct_tri_left_nd;
struct_out.correct_tri_right_nd = correct_tri_right_nd;
struct_out.error_tri_left_nd = error_tri_left_nd;
struct_out.error_tri_right_nd = error_tri_right_nd;
struct_out.aberrant_left_nd = aberrant_left_nd;
struct_out.aberrant_right_nd = aberrant_right_nd;
struct_out.proj_mean = proj_mean;

struct_out.mean_left_cd = mean_left_cd;
struct_out.mean_right_cd = mean_right_cd;

struct_out.std_left_cd = std_left_cd;
struct_out.std_right_cd = std_right_cd;

struct_out.rp_nd_proj = single(rp_nd_proj);
% struct_out.rp_nd_proj_neg = rp_nd_proj_neg;
% struct_out.rp_nd_proj_s = rp_nd_proj_s;

struct_out.N_trials_unperturbed = N_trials_unperturbed;
struct_out.right_trials_nd = right_trials_nd;
struct_out.left_trials_nd = left_trials_nd;
struct_out.inp_rec = inp_rec;
struct_out.inp_ext = inp_ext;
struct_out.inp_noise = inp_noise;

struct_out.r_in = r_in;

end

