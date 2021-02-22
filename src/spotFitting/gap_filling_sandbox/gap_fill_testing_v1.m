clear
close all

% set path to data status
data_root = 'D:\Nick\LivemRNA\Dropbox\ProcessedEnrichmentData\';
% data_status_root = 'D:\Nick\LivemRNA\Dropbox\LocalEnrichmentResults\';
% sheet_path = [data_status_root 'DataStatus.xlsx'];
% specify projects to use
project_list = {'Dl-Ven_snaBAC-mCh_v3','Dl-Ven_hbP2P-mCh_v2'};
% get list of sheets 
% [~,sheet_names]=xlsfinfo(sheet_path);
% sheet_indices = find(ismember(sheet_names,project_list));

% define some basic parameters
sim_res = 0.25;
max_jump = 5;
jump_vec_lr = -max_jump:1:max_jump;
jump_vec = -max_jump:sim_res:max_jump;

% load data
master_struct = struct;
iter = 1;
for p = 1:numel(project_list)
    d_path = [data_root project_list{p} '\nucleus_struct.mat'];
    load(d_path)
    fnames = fieldnames(nucleus_struct);
    for n = 1:numel(nucleus_struct)
        for f = 1:numel(fnames)
            master_struct(iter).(fnames{f}) = nucleus_struct(n).(fnames{f});
        end
        master_struct(iter).data_path = d_path;
        master_struct(iter).ProjectID = p;
        iter = iter + 1;
    end
    clear nucleus_struct;
end

% generate simple discretized jump transition matrix
os = ceil(numel(jump_vec_lr)/2);
jumpX_tr_counts = zeros(numel(jump_vec_lr),numel(jump_vec_lr));
jumpY_tr_counts = zeros(numel(jump_vec_lr),numel(jump_vec_lr));

xy_sigma = 1;

for m = 1:numel(master_struct)
    dx_vec = diff(round(double(master_struct(m).xPosParticle3D)));    
    dy_vec = diff(round(double(master_struct(m).yPosParticle3D)));
    for i = 1:numel(dx_vec)-1        
        dx_pre = dx_vec(i);
        dx_post = dx_vec(i+1);
        dy_pre = dy_vec(i);
        dy_post = dy_vec(i+1);
%         dx_post = find(ismembertol(jump_vec_lr,dx_vec(i+1),tol,'DataScale',1));
%         dx_pre = find(ismembertol(jump_vec_lr,dx_vec(i),tol,'DataScale',1));
%         dy_post = find(ismembertol(jump_vec_lr,dy_vec(i+1),tol,'DataScale',1));
%         dy_pre = find(ismembertol(jump_vec_lr,dy_vec(i),tol,'DataScale',1));
        if max([dx_pre dx_post dy_pre dy_post])<=jump_vec_lr(end) && min([dx_pre dx_post dy_pre dy_post]) >= jump_vec_lr(1) ...
                && ~any(isnan([dx_pre dx_post dy_pre dy_post]))                      
            jumpX_tr_counts(dx_post+os,dx_pre+os) = jumpX_tr_counts(dx_post+os,dx_pre+os) + 1;
            jumpY_tr_counts(dy_post+os,dy_pre+os) = jumpY_tr_counts(dy_post+os,dy_pre+os) + 1;
        end
    end
end

% generate testing sets 
% minDP = 20;
trace_len = 27;
lr_x_mat = [];
lr_y_mat = [];
for m = 1:numel(master_struct)
    x_vec = round(master_struct(m).xPosParticle);    
    y_vec = round(master_struct(m).yPosParticle);
    t_vec = master_struct(m).time;
    nn_ids = find(isnan(x_vec));
    nn_runs = diff(nn_ids);
    lr_ids = find(nn_runs>trace_len);
    for i = 1:numel(lr_ids)
        x_run = x_vec(nn_ids(lr_ids(i))+1:nn_ids(lr_ids(i))+trace_len);
        y_run = y_vec(nn_ids(lr_ids(i))+1:nn_ids(lr_ids(i))+trace_len);
        if sum(size(lr_x_mat)) == 0
            lr_x_mat = [x_run];
            lr_y_mat = [y_run];
        else
            lr_x_mat = [lr_x_mat ; x_run];
            lr_y_mat = [lr_y_mat ; y_run];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%% test performance of simple markov gap filler %%%%%%%%%%%%%%%%%

x_lr = repmat(jump_vec_lr,numel(jump_vec_lr),1);
y_lr = repmat(jump_vec_lr',1,numel(jump_vec_lr));
x_hr = repmat(jump_vec,numel(jump_vec),1);
y_hr = repmat(jump_vec',1,numel(jump_vec));

% interpolate count matrix to desired resolution
jumpX_tr_ct_hr = interp2(x_lr,y_lr,jumpX_tr_counts,x_hr,y_hr);
jumpY_tr_ct_hr = interp2(x_lr,y_lr,jumpY_tr_counts,x_hr,y_hr);
% focus on x for now
% make normalized transition matrices
jumpX_tr = jumpX_tr_ct_hr ./ sum(jumpX_tr_ct_hr);
jumpY_tr = jumpY_tr_ct_hr ./ sum(jumpY_tr_ct_hr);


%% calculate ss stats
x100 = jumpX_tr^100;
ssX = x100(:,1);
y100 = jumpY_tr^100;
ssY = y100(:,1);
% randomly generate gaps in run arrays
rng(123);
indices = 1:numel(lr_x_mat);
nan_indices = randsample(indices,round(numel(indices)/5),false);

x_mat_gaps = lr_x_mat;
x_mat_gaps(nan_indices) = NaN;
dx_mat_gaps = diff(x_mat_gaps')';
y_mat_gaps = lr_y_mat;
y_mat_gaps(nan_indices) = NaN;
dy_mat_gaps = diff(y_mat_gaps')';

% max aggregate movement (in pixels for now)
max_move = 125; % in pixels
% pre-generate useful quantities 
pos_option_vec = -max_move:sim_res:max_move;
pos_option_mat_fwd = repmat(pos_option_vec,numel(jump_vec),1) + repmat(jump_vec',1,numel(pos_option_vec));
pos_keep_indices_fwd = find(ismember(pos_option_mat_fwd(:),pos_option_vec));
pos_option_mat_bkd = repmat(pos_option_vec,numel(jump_vec),1) - repmat(jump_vec',1,numel(pos_option_vec));
pos_keep_indices_bkd = find(ismember(pos_option_mat_bkd(:),pos_option_vec));
bl_prob = 1 / numel(pos_option_vec);

% initialize arrays to store results
% markov chain model
x_guess_markov = NaN(size(x_mat_gaps));
x_error_markov = NaN(size(x_mat_gaps));
x_logL_markov  = NaN(size(x_mat_gaps));
x_prob_array_markov = NaN(numel(pos_option_vec),size(x_mat_gaps,2),size(x_mat_gaps,1));
% linear and spline interpolation
x_guess_lin = NaN(size(x_mat_gaps));
x_guess_spline = NaN(size(x_mat_gaps));

%% basic implementation of gap-filling algorithm 
tic
for i = 1:size(x_mat_gaps,1)
    % extract trace
    x_trace = x_mat_gaps(i,:);
    dx_trace = dx_mat_gaps(i,:);            
    
    % truncate jump outliers
    dx_trace(dx_trace>jump_vec(end)) = jump_vec(end);
    dx_trace(dx_trace<jump_vec(1)) = jump_vec(1);
    % identify gaps
    ngap_indices_x = find(isnan(x_trace));
    pgap_indices_x = find(~isnan(x_trace));
    ngap_indices_dx = find(isnan(dx_trace));
    pgap_indices_dx = find(~isnan(dx_trace));
    
    % first try simple linear and spline interpolation approaches
    % linear
    x_trace_lin = x_trace;
    x_trace_lin(ngap_indices_x) = interp1(pgap_indices_x,x_trace(pgap_indices_x),ngap_indices_x,'linear','extrap');
    x_guess_lin(i,:) = x_trace_lin;
    % spline
    x_trace_spline = x_trace;
    x_trace_spline(ngap_indices_x) = interp1(pgap_indices_x,x_trace(pgap_indices_x),ngap_indices_x,'spline','extrap');
    x_guess_spline(i,:) = x_trace_spline;
    
    % initialize fwd and bkd arrays
    fwd_array_dx = zeros(numel(jump_vec),numel(dx_trace)+1);
    fwd_array_dx(:,1) = ssX;
    bkd_array_dx = zeros(numel(jump_vec),numel(dx_trace)+1);
    bkd_array_dx(:,end) = ssX;
    
    %%%%%%%% calculate fwd probs %%%%%%%%%
    % fill in known values
    [~,dx_indices] = ismember(dx_trace(pgap_indices_dx),jump_vec);
    pind = sub2ind(size(fwd_array_dx),dx_indices,1+pgap_indices_dx);
    fwd_array_dx(pind) = 1;
    % fill in best guesses 
    for n = 1:numel(ngap_indices_dx)
        nind = ngap_indices_dx(n)+1;
        prev = fwd_array_dx(:,nind-1);
        fwd_array_dx(:,nind) = jumpX_tr*prev;
    end
    
    %%%%%%%% calculate bkd probs %%%%%%%%%
    % fill in known values
    pind = sub2ind(size(bkd_array_dx),dx_indices,pgap_indices_dx);
    bkd_array_dx(pind) = 1;
    % fill in best guesses 
    for n = fliplr(1:numel(ngap_indices_dx))
        nind = ngap_indices_dx(n);
        post = bkd_array_dx(:,nind+1)';
        bkd_array_dx(:,nind) = post*jumpX_tr;
    end
    
    %%%%%%%% apply to real space %%%%%%%%%  
    
    %%% fwd pass
    fwd_probs_x = zeros(numel(pos_option_vec),numel(x_trace));
    zero_pt = nanmean(x_trace);
    x_trace_norm = round((x_trace - zero_pt)/sim_res)*sim_res;
    
    % fill in what we know    
    [~,x_indices] = ismember(x_trace_norm(pgap_indices_x),pos_option_vec);
    fill_ind = sub2ind(size(fwd_probs_x),x_indices,pgap_indices_x);
    fwd_probs_x(fill_ind) = 1;          
    first_ind = pgap_indices_x(1);
    % uniform probs for obs prior to first spot detection 
    fwd_probs_x(:,1:first_ind-1) = bl_prob;
    
    % fill in best guesses for stuff AFTER to first obs    
    for n = ngap_indices_x(ngap_indices_x>first_ind)              
        % calculate forward weights
        x_weights = repmat(fwd_probs_x(:,n-1)',numel(jump_vec),1).*repmat(fwd_array_dx(:,n),1,numel(pos_option_vec));        
        % aggregate by position        
        fwd_probs_x(:,n) = accumarray((pos_option_mat_fwd(pos_keep_indices_fwd)+max_move)/sim_res+1,x_weights(pos_keep_indices_fwd));           
    end    
    fwd_probs_x = fwd_probs_x ./ sum(fwd_probs_x);
    
    %%% bkd pass
    bkd_probs_x = zeros(numel(pos_option_vec),numel(x_trace));    
    
    % fill in what we know        
    bkd_probs_x(fill_ind) = 1;          
    last_ind = pgap_indices_x(end);
    % uniform probs for obs after to last spot detection 
    bkd_probs_x(:,last_ind+1:end) = bl_prob;
    
    % fill in best guesses for stuff AFTER to first obs    
    for n = fliplr(ngap_indices_x(ngap_indices_x<last_ind))              
        % calculate backward weights
        x_weights = repmat(bkd_probs_x(:,n+1)',numel(jump_vec),1).*repmat(bkd_array_dx(:,n),1,numel(pos_option_vec));        
        % aggregate by position        
        bkd_probs_x(:,n) = accumarray((pos_option_mat_bkd(pos_keep_indices_bkd)+max_move)/sim_res+1,x_weights(pos_keep_indices_bkd));           
    end    
    bkd_probs_x = bkd_probs_x ./ sum(bkd_probs_x);
    
    %%%%%%%% calculate joint probabilities %%%%%%%%%%% 
    x_pos_log_raw = log(fwd_probs_x) + log(bkd_probs_x);
    x_pos_fwd_ref = logsumexp(log(fwd_probs_x) + log(fwd_probs_x),1);
    x_pos_bkd_ref = logsumexp(log(bkd_probs_x) + log(bkd_probs_x),1);
    x_denom = logsumexp(x_pos_log_raw,1);
    x_pos_log = x_pos_log_raw - x_denom;
    x_pos_prob_mat = exp(x_pos_log);
    
    % save
    x_prob_array_markov(:,:,i) = x_pos_prob_mat;
        
    % save best guess 
    x_guess_markov(i,:) = nansum(x_pos_prob_mat.*repmat(pos_option_vec',1,size(x_pos_prob_mat,2)))+ zero_pt;
    for j = 1:size(x_pos_prob_mat,2)
        x_error_markov(i,j) = nanstd(pos_option_vec+zero_pt,x_pos_prob_mat(:,j));
    end
    % likelihood "leakage" method
    x_logL_markov(i,:) = x_denom - .5*x_pos_fwd_ref - .5*x_pos_bkd_ref;
end    
toc

%% Estimate error for each gap-filling method
close all
x_vec_true = lr_x_mat(nan_indices);
x_vec_markov = x_guess_markov(nan_indices);
x_err_markov = x_error_markov(nan_indices);
x_loss_markov = x_logL_markovt(nan_indices);
x_vec_lin = x_guess_lin(nan_indices);
x_vec_spline = x_guess_spline(nan_indices);

markov_err = nanstd(x_vec_true-x_vec_markov)
markov_err_constrained = nanstd(x_vec_true(x_err_markov<2)-x_vec_markov(x_err_markov<2))
lin_err = nanstd(x_vec_true-x_vec_lin)
spline_err = nanstd(x_vec_true-x_vec_spline)
% makehistogram of errors
err_bins = linspace(-25,25);
error_fig = figure;
hold on
histogram(x_vec_true - x_vec_markov,err_bins,'Normalization','probability')
histogram(x_vec_true - x_vec_lin,err_bins,'Normalization','probability')
histogram(x_vec_true - x_vec_spline,err_bins,'Normalization','probability')

markov_only = figure;
hold on
histogram(x_vec_true - x_vec_markov,err_bins,'Normalization','probability')
histogram(x_vec_true(x_err_markov<2) - x_vec_markov(x_err_markov<2),err_bins,'Normalization','probability')

% check validity of estimated error
markov_loss_true = (x_vec_true - x_vec_markov).^2;
markov_e_fig = figure;
scatter(-x_loss_markov,  (x_vec_true - x_vec_markov).^2)
ylim([0 20])
xlim([0 2])