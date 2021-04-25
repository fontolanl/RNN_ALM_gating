function [ rp_vec, r_in_s ] = run_RNN( p, r_in )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%%

T = p.simtime_test_len;
simtime_test = p.simtime_test;
T_old = p.simtime_len;
t_dist = p.t_dist;
N_trials_distr = p.N_trials_distr;
full_dur = p.full_dur;
full = p.full;

mini_dur = p.mini_dur;
amp_dist = p.amp_dist;

stim_shape = p.stim_shape;
amp_stim = p.amp_stim;

ramp_dur = p.ramp_dur;

ramp_bsln = p.ramp_bsln;

des_out_left = p.des_out_left;

M = p.M;
M_cd = p.M_cd;

ramp = p.ramp;
sigma_ramp = p.sigma_ramp;
ramp_prefactor = p.ramp_prefactor;
stim_dur = p.stim_dur;
sigma_noise = p.sigma_noise;
sigma_dist = p.sigma_dist;
ramp_train = p.ramp_train;
fr_smooth = p.fr_smooth;
tau_noise = p.tau_noise;
dt = p.dt;
eff_dt = p.eff_dt;
h = p.h;
sigma_stim = p.sigma_stim;

N = size(des_out_left,1);

f0 = p.f0;
theta0 = p.theta0;
beta0 = p.beta0;

if isfield(p, 'init_conds')
    init_conds = p.init_conds;
else
    init_conds = des_out_left(:,1) + sigma_noise.*randn(N,1);
end

if isfield(p, 't_sim_start')
    t_sim_start = p.t_sim_start;
else
    t_sim_start = 1;
end

% if isfield(p, 'epsilon')
%     epsilon = p.epsilon;
% else
%     epsilon = 1;
% end

%% late distractor input

inp_dist_vec = zeros(T,1);
inp_dist_vec(t_dist + 1:t_dist + full_dur*full + mini_dur*(1-full)) = amp_dist.*stim_shape;

% inp_dist_temp = smooth(inp_dist_temp, fr_smooth);
% inp_dist = amp_dist.*cd_late_delay.*inp_dist_temp';
% inp_dist = amp_dist.*amp_stim.*inp_dist_temp';

% r_in_s = NaN([size(r_in), N_trials_distr]);


%% test with distractors

simtime_test_len_coarse = length(simtime_test(1:1/dt:end));

for i=1:N_trials_distr
    i
    rp_vec{i} = zeros(N,simtime_test_len_coarse);
    
    if i<=N_trials_distr/2
        r_in(2,:) = 0;
    else
        
        t_stim_interval = [1000/dt+1:1000/dt+stim_dur] + randi([-10 10],1)/dt;
        
        inp_stim_temp = zeros(T,1);
        inp_stim_temp(t_stim_interval) = amp_stim.*(1+sigma_stim*randn).*stim_shape;
        inp_stim_test = smooth(inp_stim_temp, fr_smooth);
        
        r_in(2,:) = inp_stim_test;
    end
    
    t_ramp_start = 500/dt + 1 + randi([-10 10],1)/dt;
    
    switch ramp
        case 'short'
            
            inp_ramp_test = zeros(T,1);
            inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./ramp_dur;
            %             inp_ramp_test = smooth(inp_ramp_test, fr_smooth);
            
            %             inp_ramp_test = zeros(T,1);
            %             inp_ramp_test(501:500 + 1100) = [1:1100]./1100;
            
            r_in(3,:) = ramp_prefactor.*ramp_train.*inp_ramp_test.*(1+sigma_ramp*randn) + ramp_bsln;
            r_in(3,t_ramp_start+ramp_dur:end) = r_in(3,t_ramp_start+ramp_dur-1);
            %             r_in(3,2301:end) = r_in(3,2300);
            
            r_in(3,:) = smooth(r_in(3,:), fr_smooth);
            
        case 'long'
            %             inp_ramp_temp = zeros(simtime_len,1);
            %             inp_ramp_temp(t_ramp_start:end) = [1:T-(t_ramp_start-1)]./(simtime_len-500);
            
            inp_ramp_test = zeros(T,1);
            inp_ramp_test(t_ramp_start:t_ramp_start + 3000 -1) = [1:3000]./3000;
            
            %             inp_ramp_temp = smooth(inp_ramp_temp, fr_smooth);
            
            %             inp_ramp_test = [inp_ramp_temp; inp_ramp_temp(end)...
            %                 + ([0:simtime_test(end)-simtime_len]./(simtime_len-500))'];
            
            r_in(3,:) = ramp_prefactor.*ramp_train.*inp_ramp_test.*(1+sigma_ramp*randn) + ramp_bsln;
            r_in(3,t_ramp_start + 3000:end) = r_in(3,t_ramp_start + 3000 -1);
            %             r_in(3,3501:end) = 0;
            
            r_in(3,:) = smooth(r_in(3,:), fr_smooth);
            
        case 'no_ramp'
            r_in(3,:) = ramp_bsln;
            
        case 'descending'
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
    
    inp_dist_temp = smooth(inp_dist_vec.*(1+sigma_dist*randn),fr_smooth);
    inp_dist_test = circshift(inp_dist_temp,randi([-10 10],1)/dt);
    
    %     inp_dist_test = smooth(inp_dist_test, fr_smooth);
    
    eta = randn(N,1)/tau_noise;
    xp = init_conds;
    xp(N+1:N+3) = r_in(:,t_sim_start);
    
    rp = f0./(1.0+exp(-beta0.*(xp(1:N)-theta0)));
    
    %     figure
    ti = 0;
    tii = 0;
    for t = 1:T
        ti = ti+1;
        
        eta = eta + (-eta.*dt + sigma_noise.*randn(N,1))/tau_noise;
        xp = xp + (-xp + M*[rp;r_in(:,ti)] + h + M_cd(:,N+2).*inp_dist_test(ti) + [eta;zeros(3,1)])*eff_dt;
        
        rp = f0./(1.0+exp(-beta0.*(xp(1:N)-theta0)));
        
        if mod(ti,1/dt)==0
            tii = tii + 1;
            rp_vec{i}(:,tii) = rp;
        end
        
    end
    %     r_in_s(:,:,i) = r_in;
end
r_in_s = [];
end

