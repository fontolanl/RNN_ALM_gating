function [ rp_vec, r_in_s ] = run_RNN( p, r_in )
%run_RNN Runs trained RNN forward
%   Detailed explanation goes here
%% unpack parameters

T = p.simtime_test_len;
simtime_test = p.simtime_test;
%T_old = p.simtime_len;
t_dist = p.t_dist;
N_trials_distr = p.N_trials_distr;
full_dur = p.full_dur;
full = p.full;

mini_dur = p.mini_dur;
dist_amp = p.dist_amp;

stim_shape = p.stim_shape;
stim_amp = p.stim_amp;

ramp_dur = p.ramp_dur;

ramp_bsln = p.ramp_bsln;

des_out_left = p.des_out_left;

W = p.W;
W_cd = p.W_cd;

ramp = p.ramp;
ramp_sigma = p.ramp_sigma;
ramp_prefactor = p.ramp_prefactor;
stim_dur = p.stim_dur;
noise_sigma = p.noise_sigma;
dist_sigma = p.dist_sigma;
ramp_train = p.ramp_train;
fr_smooth = p.fr_smooth;
% tau_noise = p.tau_noise;
dt = p.dt;
eff_dt = p.eff_dt;
b = p.b;
stim_sigma = p.stim_sigma;
noise_sigma_eff = p.noise_sigma_eff;

N = size(des_out_left,1);

f0 = p.f0;
theta0 = p.theta0;
beta0 = p.beta0;

if isfield(p, 'init_conds')
    init_conds = p.init_conds;
else
    init_conds = des_out_left(:,1) + noise_sigma.*randn(N,1);
end

if isfield(p, 't_sim_start')
    t_sim_start = p.t_sim_start;
else
    t_sim_start = 1;
end

%% late distractor input

inp_dist_vec = zeros(T,1);
inp_dist_vec(t_dist + 1:t_dist + full_dur*full + mini_dur*(1-full)) = dist_amp.*stim_shape;

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
        inp_stim_temp(t_stim_interval) = stim_amp.*(1+stim_sigma*randn).*stim_shape;
        inp_stim_test = smooth(inp_stim_temp, fr_smooth);

        r_in(2,:) = inp_stim_test;
    end
    
    t_ramp_start = 500/dt + 1 + randi([-10 10],1)/dt;
    
    switch ramp
        case 'endpoint'
            
            inp_ramp_test = zeros(T,1);
            inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./ramp_dur;
            
            r_in(3,:) = ramp_prefactor.*ramp_train.*inp_ramp_test.*(1+ramp_sigma*randn) + ramp_bsln;
            r_in(3,t_ramp_start+ramp_dur:end) = r_in(3,t_ramp_start+ramp_dur-1);
            
            r_in(3,:) = smooth(r_in(3,:), fr_smooth);
            
        case 'delay'
            inp_ramp_test = zeros(T,1);
            inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./(3000/dt);
            
            r_in(3,1:t_ramp_start+ramp_dur-1) = ramp_prefactor.*ramp_train.*inp_ramp_test(1:t_ramp_start+ramp_dur-1).*...
                (1+ramp_sigma*randn) + ramp_bsln;
            r_in(3,t_ramp_start+ramp_dur:end) = r_in(3,t_ramp_start+ramp_dur-1);
            
            r_in(3,:) = smooth(r_in(3,:), fr_smooth);
            
        case 'no_ramp'
            r_in(3,:) = ramp_bsln;
            
        case 'descending_delay'
            inp_ramp_test = zeros(T,1);
            inp_ramp_test(t_ramp_start:t_ramp_start + ramp_dur -1) = [1:ramp_dur]./ramp_dur;
            
            r_in(3,:) = ramp_prefactor.*ramp_train.*inp_ramp_test.*(1+ramp_sigma*randn) + ramp_bsln;
            r_in(3,t_ramp_start+ramp_dur:end) = r_in(3,t_ramp_start+ramp_dur-1);
            
            r_in(3,3000:3500) =  [500:-1:0]./500  + ramp_bsln;
            r_in(3,3501:end) = r_in(3,3500);
            
            r_in(3,:) = smooth(r_in(3,:), fr_smooth);
    end
    
    inp_dist_temp = smooth(inp_dist_vec.*(1+dist_sigma*randn),fr_smooth);
    inp_dist_test = circshift(inp_dist_temp,randi([-10 10],1)/dt);
    

    xp = init_conds;
    xp(N+1:N+3) = r_in(:,t_sim_start);
    
    rp = f0./(1.0+exp(-beta0.*(xp(1:N)-theta0)));
    
    ti = 0;
    tii = 0;
    for t = 1:T
        ti = ti+1;
        
        
        xp = xp + (-xp + W*[rp;r_in(:,ti)] + b + W_cd(:,N+2).*inp_dist_test(ti))*eff_dt + [noise_sigma_eff*randn(N,1);zeros(3,1)];

        
        rp = f0./(1.0+exp(-beta0.*(xp(1:N)-theta0)));
        
        if mod(ti,1/dt)==0
            tii = tii + 1;
            rp_vec{i}(:,tii) = rp;
        end
        
    end
end

r_in_s = [];

end

