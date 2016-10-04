function xnim_out = fit_2Dfilters( xnim, Robs, Xstims, varargin )
% Usage: xnim = xnim.fit_2Dfilters( Robs, Xstims, <train_inds>, varargin )
%
% Estimates the specified 2-D filters of XNIM model
%
% INPUTS:
%   Robs: vector of response observations (e.g. spike counts)
%   Xstims: cell array of stimuli
%   <train_inds>: index values of data on which to fit the model [default to all indices in provided data]
%   optional flags:
%       ('subs',fit_subs): set of subunits whos filters we want to optimize [default is all]
%       ('gain_funs',gain_funs): matrix of multiplicative factors, one column for each subunit
%       ('fit_offsets',fit_offsets): vector of bools, (or single bool) specifying whether
%           to fit the additive offset terms associated with each subunit
%       ('optim_params',optim_params): struct of desired optimization parameters, can also
%           be used to override any of the default values for other optional inputs
%       ('silent',silent): boolean variable indicating whether to suppress the iterative optimization display
%       ('fit_spk_hist',fit_spk_hist): boolean indicating whether to hold the spk NL filter constant
%
% OUTPUTS:
%   xnim: new nim object with optimized subunit filters

Nsubs = length(xnim.twoD_subunits); % number of 2-D subunits

% Set defaults for optional inputs
defaults.subs = 1:Nsubs; % defualt to fitting all subunits (plus -1 for spkHist filter)
defaults.gain_funs = []; % default has no gain_funs
defaults.fit_spk_hist = xnim.spk_hist.spkhstlen > 0; % default is to fit the spkNL filter if it exists
defaults.fit_offsets = false(1,Nsubs); % default is NOT to fit the offset terms
defaults.silent = false; % default is to display optimization output
option_list = {'subs','gain_funs','silent','fit_spk_hist','fit_offsets'}; % list of possible option strings

% Over-ride any defaults with user-specified values
OP_loc = find(strcmp(varargin,'optim_params')); % find if optim_params is provided as input
if ~isempty(OP_loc)
	optim_params = varargin{OP_loc+1};
	varargin(OP_loc:(OP_loc+1)) = [];
	OP_fields = lower(fieldnames(optim_params));
	for ii = 1:length(OP_fields) % loop over fields of optim_params
		if ismember(OP_fields{ii},option_list); %if the field is a valid input option, over-ride the default
			eval(sprintf('%s = optim_params.(''%s'');',OP_fields{ii},OP_fields{ii}));
			optim_params = rmfield(optim_params,OP_fields{ii}); %and remove the field from optim_params, so that the remaining fields are options for the optimizer
		end
	end
else
	optim_params = [];
end

% Now parse explicit optional input args
[train_inds,parsed_options] = NIM.parse_varargin( varargin, [], defaults );
NIM.validate_parsed_options( parsed_options, option_list );

gain_funs = parsed_options.gain_funs;
fit_subs = parsed_options.subs;
assert(all(ismember(fit_subs,1:Nsubs)),'invalid target subunits specified');
silent = parsed_options.silent;
assert(ismember(silent,[0 1]),'silent must be 0 or 1');
fit_spk_hist = parsed_options.fit_spk_hist;
assert(ismember(fit_spk_hist,[0 1]),'fit_spk_hist must be 0 or 1');

fit_offsets = zeros(1,Nsubs); 
%mod_NL_types = {xnim.subunits.NLtype}; % NL types for each targeted subunit

% Validate inputs
if ~iscell(Xstims)
	tmp = Xstims;	clear Xstims
	Xstims{1} = tmp;
end
if size(Robs,2) > size(Robs,1); Robs = Robs'; end; % make Robs a column vector
%xnim.check_inputs( Robs, Xstims, train_inds, gain_funs ); % make sure input format is correct

Nfit_subs = length(fit_subs); % number of targeted subunits
non_fit_subs = setdiff( 1:Nsubs, fit_subs ); % elements of the model held constant
spkhstlen = xnim.spk_hist.spkhstlen; % length of spike history filter
if fit_spk_hist; assert(spkhstlen > 0,'no spike history term initialized!'); end;
if spkhstlen > 0 % create spike history Xmat IF NEEDED
	Xspkhst = xnim.create_spkhist_Xmat( Robs );
else
	Xspkhst = [];
end
if ~isnan(train_inds) %if specifying a subset of indices to train model params
	for nn = 1:length(Xstims)
		Xstims{nn} = Xstims{nn}(train_inds,:); %grab the subset of indices for each stimulus element
	end
	Robs = Robs(train_inds);
	if ~isempty(Xspkhst); Xspkhst = Xspkhst(train_inds,:); end;
	if ~isempty(gain_funs); gain_funs = gain_funs(train_inds,:); end;
end

%% PARSE INITIAL PARAMETERS and add new X-matrices
%xnimtmp = xnim;
%NXstims = length(Xstims); 
[init_params,lambda_L1,sign_con] = deal([]);
NKs = zeros(1,length(fit_subs)*2); % length of each kernel to be fit
for nn = 1:length(fit_subs)
	imod = fit_subs(nn);
	cur_kern = [xnim.twoD_subunits(imod).ks{1}; xnim.twoD_subunits(imod).ks{2};];
	NKs((nn-1)*2+(1:2)) = [length(xnim.twoD_subunits(imod).ks{1}) length(xnim.twoD_subunits(imod).ks{2})];
	if (xnim.twoD_subunits(imod).Ksign_con(1) ~= 0) % add sign constraints on the filters of this subunit if needed
		sign_con(length(init_params)+(1:NKs(2*nn-1))) = xnim.twoD_subunits(imod).Ksign_con;
	end
	if (xnim.twoD_subunits(imod).Ksign_con(2) ~= 0) % add sign constraints on the filters of this subunit if needed
		sign_con(length(init_params)+NKs(2*nn-1)+(1:NKs(2*nn))) = xnim.twoD_subunits(imod).Ksign_con;
	end
	lambda_L1(length(init_params) + (1:length(cur_kern))) = xnim.twoD_subunits(imod).reg_lambdas.l1;
	init_params = [init_params; cur_kern;]; % add coefs to initial param vector
end
lambda_L1 = lambda_L1'/sum(Robs); % since we are dealing with LL/spk

Nfit_filt_params = sum(NKs); % length(init_params); % number of filter coefficients in param vector

% Add in spike history coefs
if fit_spk_hist
	init_params = [init_params; xnim.spk_hist.coefs];
	lambda_L1 = [lambda_L1; zeros(size(xnim.spk_hist.coefs))];
end

init_params = [init_params; xnim.spkNL.theta]; % add constant offset
lambda_L1 = [lambda_L1; 0];

% Calculate output of non-targets
nontarg_g = xnim.process_stimulus( Xstims, 1:length(xnim.subunits), gain_funs ); % by definition not fitting regular subunits (although this could change....)
for imod = non_fit_subs
	nontarg_g = nontarg_g + xnim.twoD_subunits(imod).process_stim( Xstims );
end
if ~fit_spk_hist && (spkhstlen > 0) % add in spike history filter output, if we're not fitting it
	nontarg_g = nontarg_g + Xspkhst*xnim.spk_hist.coefs(:);
end

%% IDENTIFY ANY CONSTRAINTS 
use_con = 0;
LB = -Inf*ones(size(init_params));
UB = Inf*ones(size(init_params));
% Constrain any of the filters to be positive or negative
if any(sign_con ~= 0)
	LB(sign_con == 1) = 0;
	UB(sign_con == -1) = 0;
	use_con = 1;
end
if fit_spk_hist % if optimizing spk history term
	% negative constraint on spk history coefs
	if xnim.spk_hist.negCon
		spkhist_inds = Nfit_filt_params + (1:spkhstlen);
		UB(spkhist_inds) = 0;
		use_con = 1;
	end
end

%% GENERATE REGULARIZATION MATRICES
Tmats = xnim.make_Tikhonov_matrices();

fit_opts = struct( 'fit_spk_hist',fit_spk_hist, 'fit_subs',fit_subs, 'fit_offsets',fit_offsets ); % put any additional fitting options into this struct
% The function we want to optimize
opt_fun = @(K) internal_LL_filters_2D( xnim,K,Robs,Xstims,Xspkhst,nontarg_g,gain_funs,Tmats,fit_opts );

% Determine which optimizer were going to use
if max(lambda_L1) > 0
	assert(~use_con,'Can use L1 penalty with constraints');
	assert(exist('L1General2_PSSas','file') == 2,'Need Mark Schmidts optimization tools installed to use L1');
	optimizer = 'L1General_PSSas';
else
	if ~use_con %if there are no constraints
		if exist('minFunc','file') == 2
			optimizer = 'minFunc';
		else
			optimizer = 'fminunc';
		end
	else
		if exist('minConf_TMP','file')==2
			optimizer = 'minConf_TMP';
		else
			optimizer = 'fmincon';
		end
	end
end
optim_params = xnim.set_optim_params( optimizer, optim_params, silent );
if ~silent; fprintf('Running optimization using %s\n\n',optimizer); end;
switch optimizer %run optimization
	case 'L1General_PSSas'
		[params] = L1General2_PSSas(opt_fun,init_params,lambda_L1,optim_params);
	case 'minFunc'
		[params] = minFunc(opt_fun, init_params, optim_params);
	case 'fminunc'
		[params] = fminunc(opt_fun, init_params, optim_params);
	case 'minConf_TMP'
		[params] = minConf_TMP(opt_fun, init_params, LB, UB, optim_params);
	case 'fmincon'
		[params] = fmincon(opt_fun, init_params, [], [], [], [], LB, UB, [], optim_params);
end
[~,penGrad] = opt_fun(params);
first_order_optim = max(abs(penGrad));
% Warn if first-order opt is off range for typical spike data (but cancel if fit is different type -- see conditions)
if (first_order_optim > xnim.opt_check_FO) && ~use_con && (max(lambda_L1) == 0) && ~strcmp(xnim.noise_dist,'gaussian') 
	warning( 'First-order optimality: %.3f, fit might not be converged.', first_order_optim );
end

%% PARSE MODEL FIT
xnim_out = xnim;
xnim_out.spkNL.theta = params(end); % set new offset parameter
if fit_spk_hist
	xnim_out.spk_hist.coefs = params(Nfit_filt_params + (1:spkhstlen));
end
kOffset = 0; % position counter for indexing param vector
for ii = 1:Nfit_subs
	for jj = 1:2
		%filtLen = length(xnim.subunits(fit_subs(ii)).filtK);
		cur_kern = params((1:NKs(2*(ii-1)+jj)) + kOffset); % grab parameters corresponding to this subunit's filters
		xnim_out.twoD_subunits(fit_subs(ii)).ks{jj} = cur_kern(:); % assign new filter values
		kOffset = kOffset + NKs(2*(ii-1)+jj);
	end
end

[LL,~,mod_internals,LL_data] = xnim_out.eval_model(Robs,Xstims,'gain_funs',gain_funs);
xnim_out = xnim_out.set_subunit_scales( mod_internals.fgint ); % update filter scales
cur_fit_details = struct('fit_type','2Dfilter','LL',LL,'filt_pen',LL_data.filt_pen,...
    'NL_pen',LL_data.NL_pen,'FO_optim',first_order_optim,'fit_subs',fit_subs);
xnim_out.fit_props = cur_fit_details; % store details of this fit
xnim_out.fit_history = cat(1,xnim.fit_history,cur_fit_details);
end




%% ************************** INTERNAL FUNCTION ***************************

function [penLL, penLLgrad] = internal_LL_filters_2D( xnim, params, Robs, Xstims, Xspkhst, nontarg_g, gain_funs, Tmats, fit_opts )
% computes the penalized LL and its gradient wrt the filters for the given nim with parameter vector params

fit_subs = fit_opts.fit_subs;
Nfit_subs = length(fit_subs); % number of targeted subs
%fit_offsets = fit_opts.fit_offsets; % which filters are we fitting offset parameters for

% USEFUL VALUES
theta = params(end); % overall model offset
gint = nan(length(Robs),Nfit_subs*2); % initialize matrix for storing filter outputs
fgint = nan(length(Robs),Nfit_subs); % initialize matrix for storing NL filter outputs
filtLen = zeros(Nfit_subs*2,1); % store the length of each (target) sub's two filters
filtKs = cell(Nfit_subs*2,1); % store the filter coefs for all (target) subs)
param_inds = cell(Nfit_subs*2,1); % this will store the index values of each subunit's filter coefs within the parameter vector
Xtarg_set = [xnim.twoD_subunits(fit_subs).Xtargs]; % vector of Xfit_subs for set of subunits being optimized

G = theta + nontarg_g; % initialize overall generating function G with the offset term and the contribution from nontarget subs

NKtot = 0;  % init filter coef counter
for ii = 1:Nfit_subs % loop over subunits, get filter coefs and their indices within the parameter vector
	for jj = 0:1  % two filters
		filtLen(2*ii-1+jj) = length(xnim.twoD_subunits(fit_subs(ii)).ks{jj+1}); % length of filter
		param_inds{2*ii-1+jj} = NKtot + (1:filtLen(2*ii-1)); % set of param indices associated with this subunit's filters
		filtKs{2*ii-1+jj} = params(param_inds{2*ii-1+jj}); % store filter coefs
		NKtot = NKtot + filtLen(2*ii-1+jj); % inc counter
		gint(:,2*ii-1+jj) = Xstims{Xtarg_set(2*ii-1+jj)}*filtKs{2*ii-1+jj};
	end
end

for jj = 1:Nfit_subs
	fgint(:,jj) = xnim.twoD_subunits(fit_subs(jj)).apply_NL(gint(:,(jj-1)*2+(1:2)));
	if isempty(gain_funs)
		G = G + fgint(:,jj);
	else
		G = G + (fgint(:,jj).*gain_funs(:,fit_subs(jj)));  
	end
end

% Add contribution from spike history filter
if fit_opts.fit_spk_hist
	G = G + Xspkhst*params( NKtot+(1:xnim.spk_hist.spkhstlen) );
end

pred_rate = xnim.apply_spkNL(G);
[penLL,LL_norm] = xnim.internal_LL(pred_rate,Robs); % compute LL and its normalization

%residual = LL'[r].*F'[g]
residual = xnim.internal_LL_deriv(pred_rate,Robs) .* xnim.apply_spkNL_deriv(G,pred_rate <= xnim.min_pred_rate);

penLLgrad = zeros(length(params),1); % initialize LL gradient
penLLgrad(end) = sum(residual);      % calculate derivatives with respect to constant term (theta)

% Calculate derivative with respect to spk history filter
if fit_opts.fit_spk_hist
	penLLgrad( NKtot+length(offset_inds)+(1:xnim.spk_hist.spkhstlen)) = residual'*Xspkhst;
end

% Calculate gradient with respect to each 2-D filter
for ii = 1:Nfit_subs
	dxdy_out = xnim.twoD_subunits(fit_subs(ii)).apply_NL_deriv2D( gint(:,(ii-1)*2+(1:2)) );
	for jj = 1:2
		target_params = param_inds{2*(ii-1)+jj}; % indices of filter coefs for current set of targeted subunits
		if isempty(gain_funs)
			penLLgrad(target_params) = bsxfun(@times, dxdy_out(:,jj), residual )' * Xstims{Xtarg_set(2*(ii-1)+jj)};
		else
			penLLgrad(target_params) = bsxfun(@times, dxdy_out(:,jj), residual.*gain_funs(:,fit_subs(ii)) )' * Xstims{Xtarg_set(2*(ii-1)+jj)};
		end	
	end
end

net_penalties = zeros(size(Xtarg_set));
net_pen_grads = zeros(length(params),1);
for ii = 1:length(Tmats) % loop over the derivative regularization matrices
	filtXinds = find(Xtarg_set == Tmats(ii).Xtarg); % set of subunits acting on the stimulus given by this Tmat
	if ~isempty(filtXinds)
		penalties = sum((Tmats(ii).Tmat * cat(2,filtKs{filtXinds})).^2);
		pen_grads = 2*(Tmats(ii).Tmat' * Tmats(ii).Tmat * cat(2,filtKs{filtXinds}));
		cur_lambdas = xnim.get_reg_lambdas(Tmats(ii).type,'fit_filters',filtXinds); % current lambdas
		net_penalties(filtXinds) = net_penalties(filtXinds) + penalties.*cur_lambdas;
		net_pen_grads(cat(2,param_inds{filtXinds})) = net_pen_grads(cat(2,param_inds{filtXinds})) + reshape(bsxfun(@times,pen_grads,cur_lambdas),[],1);
	end
end

l2_lambdas = xnim.get_reg_lambdas('fit_filters',1:(length(fit_subs)*2),'l2');
if any(l2_lambdas > 0)
	net_penalties = net_penalties + l2_lambdas.*cellfun(@(x) sum(x.^2),filtKs)';
	for ii = 1:length(un_Xtargs)
		filtXinds = find(Xtarg_set == un_Xtargs(ii)); % set of targeted subunits that act on this Xtarg
		net_pen_grads(cat(2,param_inds{filtXinds})) = net_pen_grads(cat(2,param_inds{filtXinds})) + reshape(2*bsxfun(@times,l2_lambdas(filtXinds),cat(2,filtKs{filtXinds})),[],1);
	end
end

penLL = penLL - sum(net_penalties);
penLLgrad = penLLgrad - net_pen_grads;

% CONVERT TO NEGATIVE LLS AND NORMALIZE 
penLL = -penLL/LL_norm;
penLLgrad = -penLLgrad/LL_norm;

end

