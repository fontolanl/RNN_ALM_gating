function [output_file_path, data_name] = train_RNN_ALM(N_trials_train, training, ramp, test_every, N_trials_test, n_distr, varargin)
% train_RNN_ALM.m Trains RNN with force
%               N_trials_train = max number of training trials
%               
%               training = type of behavioral trainings. 
%                          'exp_nd' = expert no distractors
%                          'exp' = expert, including PSTHs of trials with
%                                  distractors
%                          'reg' = regular mice (naive) with distractors
%                          'reg_slow_ts' = only regular mice which show
%                                          persistent graded activity
%                                          during correct trials
%
%               ramp = train with or without external ramping signal ('y'
%                      or 'n')
%
%               test_every = every N trials run a test
%               
%               N_trials_test = number of test trials
%
%               n_distr = number of PSTHs from trials containing
%                         distractors (2 if both late and early distractor 
%                         traces are considered, 3 if we include small distractor during sample)
% written by Lorenzo Fontolan

%% Read out the number of distractors and set the dt

% cd('/Users/fontolanl10/Dropbox (HHMI)/Force learning/input_data_to_RNN/')
% load('time_bins_pos.mat')

close all

% N_trials_train = str2double(N_trials_train);
%
% n_distr = str2double(n_distr);

distr = zeros(1,n_distr);

if n_distr>0
    for i=1:n_distr
        %         distr(i) = 2 + str2double(varargin{i});
        distr(i) = 2 + varargin{i};
    end
end

n_conds = numel(varargin)+2;

dt = 1; % simulation time step

dt_ratio = 5/dt; % ratio between Delta t from data (5ms) and simulation dt

%% load data

% input_file_path = '~/Dropbox (HHMI)/Force learning/current_code/code_for_sharing/input_data_to_RNN/';

input_file_path = 'D:\Dropbox (HHMI)\Force learning\current_code\code_for_sharing\input_data_to_RNN\';

switch training
    
    case 'exp_nd'
        load([input_file_path, 'v2_data_ALM_for_training_exp_wdist_0.8_3_668_100.mat'])
        des_out_list_all = {des_out_right, des_out_left};
        time_list_all = {[1001:1400]./dt,[]};
        amp_list_all = {1,0};
        
    case 'exp'
        load([input_file_path,'data_ALM_for_training_exp_wdist.mat'])
        des_out_list_all = {des_out_right, des_out_left, des_out_left_16,...
            des_out_left_08, des_out_right_16, des_out_right_08};
        time_list_all = {[1001:1400]./dt,[],[1901:2300]./dt,[2701:3100]./dt,{[1001:1400]./dt,[1901:2300]./dt},{[1001:1400]./dt,[2701:3100]./dt}};
        amp_list_all = {1,0,1,1,{1,1},{1,1}};
        
    case 'reg'
        load([input_file_path,'data_ALM_for_training_reg_wdist.mat'])
        des_out_list_all = {des_out_right, des_out_left, des_out_left_25, des_out_left_16,...
            des_out_left_08, des_out_right_16, des_out_right_08};
        time_list_all = {[1001:1400]./dt,[],[1001:1100]./dt,[1901:2000]./dt,[2701:2800]./dt,{[1001:1400]./dt,[1901:2000]./dt},{[1001:1400]./dt,[2701:2900]./dt}};
        amp_list_all = {1,0,1/3,1/3,1/3,{1,1/3},{1,1/3}};
        
    case 'reg_slow_ts'
        load([input_file_path,'data_ALM_for_training_reg_wdist_slow_ts.mat'])
        des_out_list_all = {des_out_right, des_out_left, des_out_left_25, des_out_left_16,...
            des_out_left_08, des_out_right_16, des_out_right_08};
        time_list_all = {[1001:1400]./dt,[],[1001:1100]./dt,[1901:2000]./dt,[2701:2800]./dt,{[1001:1400]./dt,[1901:2000]./dt},{[1001:1400]./dt,[2701:2900]./dt}};
        amp_list_all = {1,0,1/3,1/3,1/3,{1,1/3},{1,1/3}};
end

switch ramp
    case 'y'
        ramp_train = 1;
    case 'n'
        ramp_train = 0;
end

output_file_path = 'D:\Dropbox (HHMI)\Force learning\current_code\code_for_sharing\trained_RNN_output_data\';

% dir_name = ['Data_ALM_for_training_', training, '_', ramp, '_', num2str(distr(:))];
% 
% mkdir(dir_name)
% cd(dir_name)

des_out_list = des_out_list_all([1,2,distr]);
time_list = time_list_all([1,2,distr]);
amp_list = amp_list_all([1,2,distr]);

if n_distr == 0
    prob = [0,1/2,1];
elseif n_distr == 1
    prob = [0,1/2,3/4,1];
elseif n_distr == 2
    prob = [0,3/8,6/8,7/8,1];
elseif n_distr == 3
    prob = [0,3/9,6/9,7/9,8/9,1];
end

%cd('/Users/fontolanl10/Dropbox (HHMI)/Force learning')


%% set network parameters
rng('shuffle')

N = size(des_out_left,1); % number of recorded neurons;
disp([' Number of units: ', num2str(N)])
p = 1.0;   % fraction of connections undergoing learning
g = 1.5;	% spectral radius, g greater than 1 leads to chaotic networks (before training).
alpha = 1.0;  % initial constant to enter the covariance estiamtion (learning rate)

fr_smooth = 400/dt; %smoothing window length in ms
T = 3500/dt - fr_smooth/2; % simulation length in ms

avg_learning_tbin = 70/dt_ratio; % average number of timebins at which earning occurs

tau = 10; % integration time constant

% p_dropout = 0.0;
% p_dropconnect = 0.0;

%% Generate connectivity matrix

scaling = 1.0/sqrt(N);
W = randn(N).*g.*scaling;  % Random connectivity matrix

%% print vars on screen
disp([' N: ', num2str(N)]);
disp([' g: ', num2str(g)]);
disp([' p: ', num2str(p)]);
disp([' alpha: ', num2str(alpha,3)]);
disp([' T: ', num2str(T)]);
disp([' learn_every: ', num2str(avg_learning_tbin)]);


%% Chirp input time trace

inp_chirp_temp = zeros(simtime_len,1);
inp_chirp_temp([501:650, 1351:1500]./dt) = 1;
inp_chirp_temp = smooth(inp_chirp_temp, fr_smooth);

%% Initialize vectors and matrices

% target functions inizialization
zt = zeros(N+3,simtime_len);

%%% DC currents (bias)
b_amp = 0.0;
b_mean = b_amp*randn(N,1);
% h = zeros(N,1);
b_mean(N+1:N+3) = 0.0;
b_sigma = 0.0;

% initialization of inverse of correlation matrix
P = (1.0/alpha).*eye(N+3);
P(N+1:N+3,N+1:N+3)=0;

% chirp input
r_in = zeros(3,simtime_len);
r_in(1,:) = inp_chirp_temp;

% chirp inputs tests
r_in_t = zeros(3,simtime_len);
r_in_t(1,:) = inp_chirp_temp;

% input weights
W_avg_chirp = 0.1;
W_avg_stim = 1.0;
W_avg_ramp = 1.0;

% set sparsity of input vectors
p_chirp = 1;
p_stim = 1;
p_ramp = 1;

W_chirp = randperm(N, round(p_chirp.*N));
W_stim = randperm(N, round(p_stim.*N));
W_ramp = randperm(N, round(p_ramp.*N));

W(W_chirp, N+1) = W_avg_chirp*randn(length(W_chirp),1);
% M(1:N,N+3) = M(1:N,N+1)/max(M(1:N,N+1));
W(W_stim, N+2) = W_avg_stim*randn(length(W_stim),1);
% M(1:N,N+2) = cd_sample_data./norm(cd_sample_data); % 21 Feb 2019 changed to test Sample mode from data
% M(1:N,N+2) = cd_sample_data./max(cd_sample_data); % 21 Feb 2019 changed to test Sample mode from data
W(W_ramp, N+3) = W_avg_ramp*randn(length(W_ramp),1);
% M(1:N,N+3) = M(1:N,N+3)/max(M(1:N,N+3));
% M((M(1:N,N+3)<0)) = 0.0; %%% ramping input

% set to zero feedback weights to input units
W(N+1,:) = zeros(1,N+3);
W(N+2,:) = zeros(1,N+3);
W(N+3,:) = zeros(1,N+3);

% save initial vectors and weight matrix
W_init = W;

chirp_init_vec = W(1:N,N+1);
stim_init_vec = W(1:N,N+2);
ramp_init_vec = W(1:N,N+3);

% standard deviations of fast noise
noise_sigma = 100/N; %0.25; %100/N; % 20/N; %10/N

% standard deviations of stimulus amplitude
stim_sigma = 0.01;

% standard deviations of ramping input slope
ramp_sigma = 0.01;

% effective tau
eff_dt = dt/tau;

noise_sigma_eff = sqrt(dt).*noise_sigma./tau; %include tau in the noise amplitude to speed up

% N_syn = numel(M);
 
i_test = 0;
W_temp = W;

%% Dynamics

error_avg_train = zeros(N_trials_train,1);
trial_train_idx = zeros(N_trials_train,1);

for i = 1:N_trials_train
    tic
    i
    
    %     N_dropout = randperm(N, round(p_dropout.*N));
    %     N_dropconnect = randperm(N_syn, round(p_dropconnect.*N_syn));
    %     N_syn_nonzero = setdiff([1:N_syn],N_dropconnect);
    
    l_idx = 1; % reset learning timebin counter
    
    % get timebins for training in each trial (randomized)
    learn_every = cumsum(randi([avg_learning_tbin/2,3*avg_learning_tbin/2],...
        2*round(T./avg_learning_tbin),1))*dt_ratio;
    learn_every(learn_every<50) = [];
    
    % determine which trial type
    rand_n = rand(1);
    
    % determine trial type and generate the correct selective input (to be
    % simplified)
    gg = 1;
    while gg<n_conds+1
        
        if  (prob(gg)< rand_n) && (rand_n <= prob(gg+1))
            trial_train_idx(i) = gg;
            
            des_out = des_out_list{gg};
            inp_stim_temp = zeros(simtime_len,1);
            
            if size(time_list{gg},2)==2
                t_stim_interval_1 = time_list{gg}{1} + randi([-10 10],1);
                t_stim_interval_2 = time_list{gg}{2} + randi([-10 10],1);
                inp_stim_temp(t_stim_interval_1) = amp_list{gg}{1};
                inp_stim_temp(t_stim_interval_2) = amp_list{gg}{2};
                
            else
                t_stim_interval_1 = time_list{gg} + randi([-10 10],1);
                inp_stim_temp(t_stim_interval_1) = amp_list{gg};
            end
            
            inp_stim_temp = smooth(inp_stim_temp, fr_smooth);
            
            r_in(2,:) = inp_stim_temp.*(1+stim_sigma*randn); % selective input
            
            gg = gg+1;
            continue
        else
            gg = gg+1;
        end
    end
    
    % time of ramp start with jitter
    t_ramp_start = 501 + randi([-10 10],1);
    
    inp_ramp_temp = zeros(simtime_len,1);
    inp_ramp_temp(t_ramp_start:end) = [1:T-(t_ramp_start-1)]./(simtime_len-500);
    inp_ramp_temp = smooth(inp_ramp_temp, fr_smooth);
    
    r_in(3,:) = ramp_train.*inp_ramp_temp.*(1+ramp_sigma*randn); %ramping input - 'ramp_train' is either 1 or 0
    
    % membrane currents variable x
    x = mean(des_out(:,1:10),2) + noise_sigma.*randn(N,1);
    x(N+1:N+3) = 0;
    %     x(N_dropout) = 0;
    
    % firing rates (of recurrent units)
    r_rec = f0./(1.0 + exp(-beta0.*(x(1:N)-theta0)));
        
    % biases (DC input)
    b = b_mean.*(1.0 + b_sigma.*randn(size(b_mean)));
    
    % training signal (from electrophts data)
    des_out_temp = [des_out;zeros(3,size(des_out,2))];
    
    % empty average weight change for this trial
    dW_temp_avg = [];
    
    ti = 0;
    
    % input currents
    x_in = zeros(N+3,T);
        
    for t = simtime
        ti = ti+1;
        
        x_in(:,ti) = W_temp*[r_rec;r_in(:,ti)] + b;
        
        x = x + (-x + x_in(:,ti))*eff_dt + [noise_sigma_eff*randn(N,1);zeros(3,1)];
        
        %         x(N_dropout) = 0;
                
        r_rec = f0./(1.0 + exp(-beta0.*(x(1:N)-theta0)));
        
        % learning
        if (mod(ti, learn_every(l_idx)) == 0)
            
            l_idx = l_idx + 1;
            
            % compute currents
            r = [r_rec;r_in(:,ti)];
            zt(:,l_idx-1) = W_temp*r + b;
         
            % update inverse correlation matrix and related quantities
            k = P * r;
            rPr = r'*k;
            c = 1.0./(1.0 + rPr);
            P = P - c*(k*k');
                        
            %             zt(N_dropout,l_idx-1) = 0;
            %             des_out_temp(N_dropout,ti/dt_ratio) = 0;
            
            % update the recurrent weight matrix using the error at each neuron
            dW_temp = c*(des_out_temp(:,ti/dt_ratio) - zt(:,l_idx-1))*k';
            W_temp = W_temp + dW_temp;
            
            % compute the average absolute change in weights
            dW_temp_avg(l_idx-1) = mean(abs(dW_temp(:)));
            
            %             M_temp(N_dropconnect) = 0;
        end
        
    end
    
    %     zt(N_dropout,:) = 0;
    
    % compute mean training error of trial i
    error_avg_train(i) = mean(sum(abs(zt(1:N,1:l_idx-1)...
        - des_out(:,learn_every(1:l_idx-1)/dt_ratio)))/size(des_out,2));
    
    %     error_avg_train(i) = mean(sum(abs(x_in(1:N,1:dt_ratio:end)...
    %         -des_out))/size(des_out,2));
    
    % update plot with average error and average absolute weight change
    if mod(i, 1) == 0
        figure(1)
        subplot(2,1,1)
        hold on
        plot(i,error_avg_train(i),'o','Color',[1-amp_list{trial_train_idx(i)},0,amp_list{trial_train_idx(i)}])
        xlabel('Test trial')
        ylabel('MSE (a.u.)')
        subplot(2,1,2)
        hold on
        plot(i,abs(mean(dW_temp_avg)),'o','Color',[0.5,0,0.5])
        xlabel('Test trial')
        ylabel('|\Delta W|')
        drawnow
        %                 pause(0.01);
    end
    
    %     dM_avg(i) = mean(dM_temp_avg);
    
    %     M(N_syn_nonzero) = M_temp(N_syn_nonzero);
    
    
    %% test trials (same structure as training trials)
    if mod(i,test_every) == 0
        
        i_test = i_test + 1;
        
        error_avg_test = zeros(N_trials_test,1);
        
        %% test
        disp(['Now testing... at i=',num2str(i), ' please wait.']);
        
        for ii=1:N_trials_test
            
            zt_t_vec = zeros(N+3,simtime_len);
            
            r_in_t(3,:) = ramp_train.*inp_ramp_temp.*(1+ramp_sigma*randn);
            
            rand_n = rand(1);
            
            x_t_vec = zeros(N+3,T);
            
            %%
            gg = 1;
            while gg<n_conds+1
                
                if  (prob(gg)< rand_n) && (rand_n <= prob(gg+1))
                    trial_test_idx(ii) = gg;
                    des_out = des_out_list{gg};
                    inp_stim_temp = zeros(simtime_len,1);
                    
                    if size(time_list{gg},2)==2
                        t_stim_interval_1 = time_list{gg}{1} + randi([-10 10],1);
                        t_stim_interval_2 = time_list{gg}{2} + randi([-10 10],1);
                        inp_stim_temp(t_stim_interval_1) = amp_list{gg}{1};
                        inp_stim_temp(t_stim_interval_2) = amp_list{gg}{2};
                        
                    else
                        t_stim_interval_1 = time_list{gg} + randi([-10 10],1);
                        inp_stim_temp(t_stim_interval_1) = amp_list{gg};
                    end
                    
                    inp_stim_temp = smooth(inp_stim_temp, fr_smooth);
                    
                    r_in_t(2,:) = inp_stim_temp.*(1+stim_sigma*randn);
                    
                    gg = gg+1;
                    
                    continue
                else
                    gg = gg+1;
                end
            end
            
            x_t = mean(des_out(:,1:10),2) + noise_sigma.*randn(N,1);
            x_t(N+1:N+3) = 0;
            
            r_t = f0./(1.0 + exp(-beta0.*(x_t(1:N)-theta0)));
            
            ti_t = 0;
            for tt = simtime
                ti_t = ti_t+1;
                
                x_t = x_t + (-x_t + W_temp*[r_t;r_in_t(:,ti_t)])*eff_dt + b + [noise_sigma_eff*randn(N,1); zeros(3,1)];
                
                x_t_vec(:,ti_t) = x_t;
                
                r_t = f0./(1.0 + exp(-beta0.*(x_t(1:N)-theta0)));
                
                zt_t = W_temp*[r_t;r_in_t(:,ti_t)] + b;
                
                zt_t_vec(:,ti_t) = zt_t;
                
            end
            
            error_avg_test(ii) = mean(sum(abs(zt_t_vec(1:N,1:dt_ratio:end)-des_out))./size(des_out,2));
             
%             disp(['Testing MAE: ' num2str(error_avg_test(ii))]);
            
            if mod(ii, 1) == 0
                figure(i_test + 2)
                title(['Testing after ', num2str(i), ' trials'])
                hold on
                plot(ii,error_avg_test(ii),'o','Color',[1-amp_list{trial_test_idx(ii)},0,amp_list{trial_test_idx(ii)}])
                xlabel('Test trial')
                ylabel('MSE (a.u.)')
                drawnow
                %                 pause(0.01);
            end
            
        end
        
        % plot test error 
        figure(2)
        hold on
        plot(i_test,mean(error_avg_test),'o','Color',[0.2,0.8,0.2])
        xlabel('Testing Session')
        ylabel('MSE (a.u.)')
        % save data after testing
        
        W = W_temp;
        
                
%         data_name = ['N=', num2str(i),'_g=', num2str(g),'_a=', num2str(alpha),...
%             '_sn=',num2str(noise_sigma), '_ss=',num2str(stim_sigma),  '_sr=',num2str(ramp_sigma),...
%             '_hin=', num2str(b_amp), '_Wc=', num2str(W_avg_chirp), '_Ws=', num2str(W_avg_stim), '_Wr=', num2str(W_avg_ramp),'_',...
%             datestr(datetime('now','Format','yyyy-MM-dd''T''HH:mmXXX'))];

        data_name = ['N=', num2str(i),'_g=', num2str(g),'_a=', num2str(alpha),...
            '_sn=',num2str(noise_sigma), '_ss=',num2str(stim_sigma),  '_sr=',num2str(ramp_sigma),...
            '_hin=', num2str(b_amp), '_Wc=', num2str(W_avg_chirp), '_Ws=', num2str(W_avg_stim), '_Wr=', num2str(W_avg_ramp)];
        
        %         close all
        save([output_file_path, data_name,'.mat'],'-v7.3')
        
    end    
    
    
%     disp(['Training MAE: ',num2str(error_avg_train(i))]);
    
    toc
end