classdef XNIM < NIM
% Class implementation of NIM that incorporates 2D non-parametric
% nonlinearities for combining two signals.

% TODO
% - reg should only be d2xt or l2

properties
    twoD_subunits;
    % Inherits from NIM:
    % - spkNL
    % - subunits
    % - stim_params
    % - noise_dist
    % - spk_hist
    % - fit_props
    % - fit_history
end

methods % defined in separate files
	xnim = fit_2Dfilters( xnim, Robs, Xstims, varargin );
	xnim_out = twoDreg_path( xnim, Robs, Xs, Uindx, XVindx, varargin );
end

%% *************************** constructor ********************************
methods

	function xnim = XNIM( nim )
	% Usage: xnim = XNIM( nim )

		if nargin == 0
			% handle the no-input-argument case by returning a null model. This is important when initializing arrays of objects  
			return 
		end
    
		nim_fields = fields(nim);
		for nn = 1:length(nim_fields)
			xnim.(nim_fields{nn}) = nim.(nim_fields{nn});
		end
		xnim.twoD_subunits = [];
    
		% Change superclass properties
		xnim.version = '0.0';
	end % method
    
end

%% *************************** setting methods ***************************
methods

	function xnim = set_2Dreg_params( xnim, varargin )
	% Usage: xnim = xnim.set_2Dreg_params( varargin )
	% Sets a desired set of regularization parameters to specified values, apply to specified set of 2D-subunits
	%
	% INPUTS:
	%   optional flags:
	%       ('subs',sub_inds): set of subunits to apply the new reg_params for (default = ALL)
	%       ('filts',filt_inds): set of filters within each 2-d subunit to apply the new reg_params for (default = both)
	%       ('lambda_type',lambda_val): first input is a string specifying the type of regularization,
	%            e.g. 'd2t' for temporal smoothness. This must be followed by a scalar or vector of
	%            length (Nsubs) giving the associated lambda values
	%
	% OUTPUTS:
	%   xnim: XNIM model object

		sub_inds = 1:length(xnim.twoD_subunits); %default is to apply the change to all subunits
		filt_inds = [1 2];
		
		% INPUT PARSING
		%[~, parsed_varargin] = NIM.parse_varargin(varargin, {'subs','filts','lambda_type'}, defaults);
	  %sub_inds = parsed_varargin.subs;
		%assert(all(ismember(sub_inds,1:length(xnim.twoD_subunits))),'invalid target 2-D subunits specified');
		%filt_inds = parsed_varargin.filts;
		%assert(all(ismember(sub_inds,1:2)),'invalid target filters specified');
		j = 1;
		reg_types = {}; reg_vals = [];
		while j <= length(varargin)
			switch lower(varargin{j})  
				case 'subs'
					sub_inds = varargin{j+1};
					assert(all(ismember(sub_inds,1:length(xnim.twoD_subunits))),'invalid target 2-D subunits specified');
				case 'filts'
					filt_inds = varargin{j+1};
					assert(all(ismember(filt_inds,1:2)),'invalid target filters specified');
				case xnim.allowed_reg_types
					reg_types = cat(1,reg_types,lower(varargin{j}));
					cur_vals = varargin{j+1};
					reg_vals = cat(2,reg_vals, cur_vals(:)); % build up a [KxP] matrix, where K is the number of subunits and P is the number of reg types
				otherwise
					error('Invalid input flag');
			end
			j = j + 2;
		end
		
		if size(reg_vals,1) == 1 % if reg_vals are specified as scalars, assume we want the same for all subuntis
			reg_vals = repmat(reg_vals,length(sub_inds),1);
		end
		
		if isempty(reg_vals)
			warning('No regularization values specified, no action taken');
		end
		
		for ii = 1:length(reg_types)
			assert(all(reg_vals(:,ii) >= 0),'regularization hyperparameters must be non-negative');
			for jj = 1:length(sub_inds)
				for kk = 1:length(filt_inds)
					xnim.twoD_subunits(sub_inds(jj)).reg_lambdas.(reg_types{ii})(filt_inds(kk)) = reg_vals(jj,ii);
				end
			end
		end	
	end
	
end % setting methods

%% *************************** fitting methods ****************************
methods 

	function xnim_out = fit_filters( xnim, Robs, Xstims, varargin )
	% Usage: xnim = xnim.fit_filters( Robs, Xstims, Uindx, varargin )
	% Overloaded fit_filters that fits all non-2D subunits. It just calls NIM.fit_filters, but with addition Xstim
	% corresponding to contributions of all 2D subunits (for which it fits the gains and translates this to the 
	% output of each NL2d

		Nstims = length(Xstims);
		N2d = length(xnim.twoD_subunits);

		% Convert to NIM with 2D-subunit output calculated
		[nimtmp,Xstims_mod] = xnim.convert2NIM( Xstims );

		% Fit filters of all regular subunits
		nimtmp = nimtmp.fit_filters( Robs, Xstims_mod, varargin{:} );

		xnim_out = XNIM( nimtmp ); % translate back to XNIM

		% Copy 2D-subunits, with potential gain changes
		if N2d > 0
			xnim_out.subunits = xnim_out.subunits(1:end-1);  % remove extra gain-subunit
			xnim_out.TwoDsubunits = xnim.TwoDsubunits;
			for nn = 1:N2d
				xnim_out.TwoDsubunits.NL2d = xnim_out.TwoDsubunits.NL2d * nimtemp.subunits(Nstims+nn).filtK;
			end
		end 
  
	end % method

% 	function xnim_out = fit_2Dfilters( xnim, Robs, Xstims, varargin )
% 	% Usage: xnim = xnim.fit_2Dfilters( Robs, Xstims, varargin )
% 
% 		N2Dsubs = length(xnim.twoD_subunits);
% 		if N2Dsubs == 0
% 			xnim_out = xnim;
% 			warning( 'No 2-D subunits to fit.' )
% 			return
% 		end
% 		if ~iscell(Xstims)
% 			tmp = Xstims; clear Xstims
% 			Xstims{1} = tmp;
% 		end
% 		
% 		% Make NIM with only 2-d subunits. All 1-D subunits become additional X-matrix (not processed by fit-filters)
% 		nimtmp = NIM();
% 		nim_fields = fields(xnim);
% 		for nn = 1:length(nim_fields)
% 			if ~(strcmp(nim_fields{nn},'twoD_subunits') || strcmp(nim_fields{nn},'subunits'))
% 				nimtmp.(nim_fields{nn}) = xnim.(nim_fields{nn});
% 			end
% 		end
% 		% Assign twoD_subunits to subunits argument
% 		nimtmp.subunits = xnim.twoD_subunits;
% 		
% 		% Make new Xstims based on twoD subunits
% 		for nn = 1:N2Dsubs
% 			Xs{nn} = [Xstims{xnim.twoD_subunits(nn).Xtargs{1}} Xstims{xnim.twoD_subunits(nn).Xtargs{2}}];
% 			stim_par_list(nn) = xnim.stim_params(xnim.twoD_subunits(nn).Xtargs{1});
% 			stim_par_list(nn).dims(1) = size(Xs{nn},2);
% 			nimtmp.subunits(nn).Xtarg = nn;
% 			nimtmp.subunits(nn).filtK = [xnim.twoD_subunits(nn).ks{1}; xnim.twoD_subunits(nn).ks{2};];
% 		end
% 		nimtmp.stim_params = stim_par_list;
% 		% Make single Xstim for rest of subunit output (and add subunit)
%     
% 	end  % method

	function xnim_out = fit_NL2d(xnim, Robs, Xstims, varargin)
	% Usage: xnim = xnim.fit_NL2d(Robs, Xstims, Uindx, varargin)
	% Estimates 2D nonlinearities of the XNIM model
	%
	% INPUTS:
	%   Robs:           vector of response observations (e.g. spike counts)
	%   Xstims:         cell array of stimuli
	%   <train_inds>:   index values of data on which to fit the model 
	%                   [default to all indices in provided data]
	%
	%   optional flags:
	%       ('subs',fit_subs): set of twoD_subunits whos coefficients we want to optimize [default is all]
	%       ('gain_funs',gain_funs): matrix of multiplicative factors, one column for each subunit
	%       ('fit_offsets',fit_offsets): vector of bools, (or single bool)  specifying whether to fit the additive 
	%           offset terms associated with each subunit
	%       ('optim_params',optim_params): struct of desired optimization parameters, can also be used to override 
	%           any of the default values for other optional inputs
	%       ('silent',silent): boolean variable indicating whether to suppress the iterative optimization display
	%       ('fit_spk_hist',fit_spk_hist): boolean indicating whether to hold the spk NL filter constant
	%
	% OUTPUTS:
	%   xnim: new xnim object with optimized 2d nonlinearities
    
		defaults.subs = 1:length(xnim.twoD_subunits);
		[~, parsed_varargin, modvarargin] = NIM.parse_varargin(varargin, {'subs'}, defaults);
		twoD_subs = parsed_varargin.subs; 
    
		if ~iscell(Xstims)
			tmp = Xstims; clear Xstims
			Xstims{1} = tmp;
		end
		
		% twoD_subunits will be turned into normal NIM subunits and appended to subunits array; update 'subs' field to reflect this
		modvarargin{end+1} = 'subs';
		modvarargin{end+1} = length(xnim.subunits) + parsed_varargin.subs;
    
		% Reformulate XNIM into an NIM for fitting; let fit_filters take care of calculating outputs of non-target subunits
		[nim, Xstims] = xnim.format4NIMfitting(Xstims);
  
		% Fit 2-D nonlinearity coefficients
		nim = nim.fit_filters(Robs, Xstims, modvarargin{:});
    
		% -------------- Translate everything back into XNIM ------------------
		xnim_out = xnim;
		for i = 1:length(twoD_subs)
			xnim_sub_indx = twoD_subs(i);
			nim_sub_indx = length(xnim.subunits) + parsed_varargin.subs(i);
			xnim_out.twoD_subunits(xnim_sub_indx).NL2d = reshape(...
					nim.subunits(nim_sub_indx).filtK, ...
					length(xnim_out.twoD_subunits(xnim_sub_indx).ticks{1}), ...
					length(xnim_out.twoD_subunits(xnim_sub_indx).ticks{2}));
		end
		xnim_out.spkNL = nim.spkNL;
		xnim_out.spk_hist = nim.spk_hist;
		% Modify fit history
		xnim_out.fit_props = nim.fit_props;
		xnim_out.fit_props.fit_type = 'NL2d';
		xnim_out.fit_history = cat(1, xnim_out.fit_history, nim.fit_history(end));
		xnim_out.fit_history(end).fit_type = 'NL2d';    
    
	end % method

	function xnim_out = reg_path_NL2d(xnim, Robs, Xstims, Uindx, XVindx, varargin)
	% Usage: xnim = xnim.reg_path_NL2d(Robs, Xstims, Uindx, XVindx, varargin)
	% Estimates 2D nonlinearities of the XNIM model while searching for optimal regularization parameters using XVindx.
	%
	% INPUTS:
	%   Robs:   vector of response observations (e.g. spike counts)
	%   Xstims: cell array of stimuli
	%   Uindx:  index values of data on which to fit the model
	%   XVindx: index values of data on which cross-validation is evaluated
	%
	%   optional flags:
	%       ('subs',fit_subs): set of twoD_subunits whos coefficients we want to optimize [default is all]
	%       ('gain_funs',gain_funs): matrix of multiplicative factors, one column for each subunit
	%       ('fit_offsets',fit_offsets): vector of bools, (or single bool) specifying whether to fit the additive 
	%           offset terms associated with each subunit
	%       ('optim_params',optim_params): struct of desired optimization parameters, can also be used to override
	%           any of the default values for other optional inputs
	%       ('silent',silent): boolean variable indicating whether to suppress the iterative optimization display
	%       ('fit_spk_hist',fit_spk_hist): boolean indicating whether to hold the spk NL filter constant
	%
	% OUTPUTS:
	%   xnim: new xnim object with optimized 2d nonlinearities
    
		defaults.subs = 1:length(xnim.twoD_subunits);
		[~, parsed_varargin, modvarargin] = NIM.parse_varargin(varargin, {'subs'}, defaults);
		twoD_subs = parsed_varargin.subs; % record for later
    
		% twoD_subunits will be turned into normal NIM subunits and appended to subunits array; update 'subs' field to reflect this
		modvarargin{end+1} = 'subs';
		modvarargin{end+1} = length(xnim.subunits) + parsed_varargin.subs;
    
		% Reformulate XNIM into an NIM for fitting; let fit_filters take care of calculating outputs of non-target subunits
		[nim, Xstims] = xnim.format4NIMfitting(Xstims);
    
		% Fit 2d nonlinearity coefficients
		nim = nim.reg_path( Robs, Xstims, Uindx, XVindx, modvarargin{:} );
    
		% -------------- Translate everything back into XNIM ------------------
		xnim_out = xnim;
		for i = 1:length(twoD_subs)
			xnim_sub_indx = twoD_subs(i);
			nim_sub_indx = length(xnim.subunits) + parsed_varargin.subs(i);
			xnim_out.twoD_subunits(xnim_sub_indx).NL2d = reshape( nim.subunits(nim_sub_indx).filtK, ...
					length(xnim_out.twoD_subunits(xnim_sub_indx).ticks{1}), ...
					length(xnim_out.twoD_subunits(xnim_sub_indx).ticks{2}) );
			% translate regularization back into nld2
			xnim_out.twoD_subunits(xnim_sub_indx).reg_lambdas.nld2 = nim.subunits(nim_sub_indx).reg_lambdas.d2xt;
			%xnim_out.twoD_subunits(xnim_sub_indx).reg_lambdas.d2xt = xnim.twoD_subunits(xnim_sub_indx).reg_lambdas.d2xt;
		end
		xnim_out.spkNL = nim.spkNL;
		xnim_out.spk_hist = nim.spk_hist;
		% Modify fit history
		xnim_out.fit_props = nim.fit_props;
		xnim_out.fit_props.fit_type = 'NL2d';
		xnim_out.fit_history = cat(1, xnim_out.fit_history, nim.fit_history(end));
		xnim_out.fit_history(end).fit_type = 'NL2d';    
    
	end % method
    
end
%% *************************** helper methods *****************************
methods
    
	function xnim = add_2d_subunit( xnim, subunit_1, subunit_2, Xstims, varargin )
	% Usage: xnim = xnim.add_2d_subunit( subunit_1, subunit_2, Xstims, varargin )
	% Takes two SUBUNIT inputs (probably from an NIM) and converts them into a TWODSUBUNIT
	%
	% INPUTS:
	%   subunit_1:  SUBUNIT object
	%   subunit_2:  SUBUNIT object
	%   Xstims:     cell array of stimulus matrices
	%
	%   optional arguments:
	%       'interaction':  string specifying initial type of interaction.
	%                       'add' or 'mult'. Default is 'add'.
	%       'Nbins':        2x1 cell array or scalar defining number of 
	%                       bins to represent each dimension.
	%
	% OUTPUTS:
	%   xnim:   updated XNIM object
    
		defaults.interaction = 'add';
		fields_to_remove = 'interaction';
		[~, parsed_options, modvarargin] = NIM.parse_varargin(varargin, fields_to_remove, defaults);
		ks{1} = subunit_1.filtK;
		ks{2} = subunit_2.filtK;
		Xtargs = [subunit_1.Xtarg subunit_2.Xtarg];
		if ~iscell(Xstims)
			tmp = Xstims; clear Xstims
			Xstims{1} = tmp;
		end
	
		% Initialize 2D subunit
		if ~isempty(xnim.twoD_subunits) % dumb kluge to initialize empty list
			xnim.twoD_subunits(end+1) = TWODSUBUNIT(ks, Xstims, 'Xtargs', Xtargs, modvarargin{:});
		else
			xnim.twoD_subunits = TWODSUBUNIT(ks, Xstims, 'Xtargs', Xtargs, modvarargin{:});
		end
	
		% Initialize nonlinearity
		xnim.twoD_subunits(end) = xnim.twoD_subunits(end).init_NL2d(parsed_options.interaction);
        
	end
    
	function xnim = convert_to_2d_subunit(xnim, sub_indx_1, sub_indx_2, Xstims, varargin)
	% Usage: xnim = xnim.add_2d_subunit(sub_indx_1, sub_indx_2, Xstims, varargin)
	% Takes two SUBUNITS from the current XNIM object and converts them into a TWODSUBUNIT, and deletes the original SUBUNITS
	%
    % INPUTS:
    %   sub_indx_1:  index into subunits property of xnim
    %   sub_indx_2:  index into subunits property of xnim
    %   Xstims:      cell array of stimulus matrices
    %
    %   optional arguments:
    %       'interaction':  string specifying initial type of interaction.
    %                       'add' or 'mult'. Default is 'add'.
    %       'Nticks':       2x1 cell array or scalar defining number of 
    %                       ticks to represent each dimension.
    %
    % OUTPUTS:
    %   xnim:   updated XNIM object
    
    % make sure there are enough subunits
    if length(xnim.subunits) < 2
        error('Not enough subunits to convert to TWODSUBNIT')
    end
    
    defaults.interaction = 'mult';
    defaults.Nticks = 21;
    fields_to_remove = {'interaction'};
    [~, parsed_options, modvarargin] = NIM.parse_varargin(varargin, fields_to_remove, defaults);
    
    ks{1} = xnim.subunits(sub_indx_1).filtK;
    ks{2} = xnim.subunits(sub_indx_2).filtK;
    Xtargs = [xnim.subunits(sub_indx_1).Xtarg xnim.subunits(sub_indx_2).Xtarg];
    
    % initialize 2D subunit
    xnim.twoD_subunits = cat(1, xnim.twoD_subunits, TWODSUBUNIT(ks, Xstims, 'Xtargs', Xtargs, modvarargin{:}));
    
    % initialize nonlinearity
    xnim.twoD_subunits(end) = xnim.twoD_subunits(end).init_NL2d(parsed_options.interaction);
    
    % get rid of converted SUBUNITS; make sure to do so in proper order
    if sub_indx_1 > sub_indx_2
        xnim.subunits(sub_indx_1) = [];
        xnim.subunits(sub_indx_2) = [];
    elseif sub_indx_1 < sub_indx_2
        xnim.subunits(sub_indx_2) = [];
        xnim.subunits(sub_indx_1) = [];
    else
        xnim.subunits(sub_indx_1) = []; % same subunit
    end
    
    % turn subunits property into a 0x0 instead of 1x0 array if empty
    if isempty(xnim.subunits)
        xnim.subunits = [];
		end
    
	end
    
	function [LL, pred_rate, mod_internals, LL_data] = eval_model( xnim, Robs, Xstims, varargin )
	% Usage: [LL, pred_rate, mod_internals, LL_data] = xnim.eval_model( Robs, Xstims, <eval_inds>, varargin ) 
	%
	% Overloaded for XNIM: evaluates model as NIM, but for mod_internals it will have the twoDsubunit outputs
	% as the last entries in fgint, leave the corresponding gint to be the same, and separately have gint2d
	
		[nim_tmp, Xs_mod] = xnim.convert2NIM( Xstims );
    
		[LL, pred_rate, mod_internals, LL_data] = nim_tmp.eval_model( Robs, Xs_mod, varargin{:} );

		if ~isempty(xnim.twoD_subunits)
			N2d = length(xnim.twoD_subunits);
			% mod_internals.fgint2d = mod_internals.fgint(:,end);
			% mod_internals.gint = mod_internals.gint(:,1:end-1);
			%mod_internals.fgint = mod_internals.fgint(:,1:end-1);
			[XVindx,~] = NIM.parse_varargin( varargin );
			if isempty(XVindx)
				XVindx = 1:size(mod_internals.gint,1);
			end
			mod_internals.gint2d = zeros(length(XVindx),N2d*2);
			if ~iscell(Xstims)
				tmp = Xstims; clear Xstims
				Xstims{1} = tmp;
			end
			for nn = 1:N2d
				for mm = 1:2
					mod_internals.gint2d(:,(nn-1)*2+mm) = Xstims{xnim.twoD_subunits(nn).Xtargs(mm)}(XVindx,:) * xnim.twoD_subunits(nn).ks{mm};
				end
			end
		else
			mod_internals.gint2d = [];
		end
			
	end % method
    
	function  fig_handle = display_model(xnim, varargin)
	% Usage: [] = xnim.display_model(varargin)
	% Displays all parts of model as single plot (with multiple subplots)
	%
	% INPUTS:  
	%  Optional flags:
	%    'Xstims': to enter Xstims for calculation of G-distributions
	%    'Robs': to enter Robs, as required for G-distributions where there is a spike-history term
	%    'mod_outs': output of eval_model that gives required internal parameters in place of Xstims and Robs
	%    'colormap': to specify colormap for any image plots in filter displays
	%    'time_rev': time-reverse temporal plots
	%    'dt': enter to plot actual time versus time in binned units
	
		[~, parsed_options, modvarargin] = NIM.parse_varargin( varargin, {'Xstims','Robs','mod_outs'} );
		valid_list = {'Xstims', 'Robs', 'mod_outs', 'colormap', 'single', ...
									'color', 'dt', 'time_rev', 'xt-separable', 'xt-spatial', ...
									'sign', 'y_offset', 'no_axes_space', 'no_axes_time'};
		NIM.validate_parsed_options(parsed_options, valid_list);
		Nmods = length(xnim.twoD_subunits);
		mod_outs = [];

		if isfield(parsed_options, 'Xstims')
			if iscell(parsed_options.Xstims)
				Xstims = parsed_options.Xstims;
			else
				Xstims{1} = parsed_options.Xstims;
			end			
			Robs = [];
			if isfield(parsed_options, 'Robs') 
				Robs = parsed_options.Robs;
			end
			[~,~,mod_outs] = xnim.eval_model(Robs, Xstims);
		elseif isfield(parsed_options, 'mod_outs')
			mod_outs = parsed_options.mod_outs;
		end

    % Will be spike-history or spkNL plot?
		extra_plots = [(xnim.spk_hist.spkhstlen > 0) ~isempty(mod_outs)]; 

		if sum(extra_plots) == 0
			Nrows = Nmods;
			Ncols = 3;
		else
			% then need extra column (and possibly extra row)
			Nrows = max([Nmods sum(extra_plots)]);
			Ncols = 4;
		end

		if nargout > 0
			fig_handle = figure;
		else
			figure;
		end
		
		% Plot subunit info
		for nn = 1:Nmods
			dims{1} = xnim.stim_params(xnim.twoD_subunits(nn).Xtargs(1)).dims;
			dims{2} = xnim.stim_params(xnim.twoD_subunits(nn).Xtargs(2)).dims;
			xnim.twoD_subunits(nn).display_filter( dims, [Nrows Ncols (nn-1)*Ncols+1], modvarargin{:} );
			subplot(Nrows, Ncols, (nn-1)*Ncols+3);
			if isempty(mod_outs)
				xnim.twoD_subunits(nn).display_NL();
			else
				xnim.twoD_subunits(nn).display_NL(mod_outs.gint(:,nn));
			end
		end

		% Plot spkNL
		if sum(extra_plots) == 0
			return
		end

    subplot(Nrows, Ncols, Ncols);
    if extra_plots(2) > 0
        xnim.display_spkNL(mod_outs.G);
        title( 'Spk NL' )
        if extra_plots(1) > 0
            subplot(Nrows, Ncols, 2*Ncols);
        end
    end
    if extra_plots(1) > 0
        xnim.display_spike_history();
		end
    
	end
    
end
%% *************************** hidden methods *****************************
methods (Hidden)

	function [nim, Xstims] = convert2NIM( xnim, Xstims )
	% Usage: nim = xnim.convert2NIM( <Xstims> )
	% Allows NIM function calls (that are not overloaded by XNIM). Sets
	% output of each twoD_subunit to be a univariate signal that a scale
	% factor can be fit to, if desired.

		if nargin < 2
			Xstims = [];
		elseif ~iscell(Xstims)
			tmp = Xstims; clear Xstims
			Xstims{1} = tmp;
		end

		nim = NIM();
		nim_fields = fields(xnim);
		for nn = 1:length(nim_fields)
			if ~strcmp(nim_fields{nn},'twoD_subunits')
				nim.(nim_fields{nn}) = xnim.(nim_fields{nn});
			end
		end

		% Add additional subunit to NIM to capture 2D subunit contributions
		N2d = length(xnim.twoD_subunits);
		if (N2d > 0) && ~isempty(Xstims)
        
			% Make X-matrix for 2D-subunits
			Nstims = length(Xstims);
			stim_params = NIM.create_stim_params([1 1 1]);

			for nn = 1:N2d
				Xstims{Nstims+nn} = xnim.twoD_subunits(nn).process_stim(Xstims);
				nim.stim_params(Nstims+nn) = stim_params;
				% convert to NIM subunits
				nim = nim.add_subunits({'lin'}, 1, 'xtargs', Nstims+nn, 'init_filts', {[1]});
			end
		end

	end

	function [nim, Xstims] = format4NIMfitting( xnim, Xstims )
	% Usage: [nim, Xstims] = xnim.format4NIMfitting( <Xstims> )
	% Turns twoD_subunits into linear subunits of an NIM so the 2D-NL can be fit using the standard
	% fit_filters routine in the NIM. Regularization should be in the nld2 reg_param

		if nargin < 2
			Xstims = [];
		elseif ~iscell(Xstims)
			tmp = Xstims; clear Xstims
			Xstims{1} = tmp;
		end

		nim = NIM();
		nim_fields = fields(xnim);
		for nn = 1:length(nim_fields)
			if ~strcmp(nim_fields{nn},'twoD_subunits')
				nim.(nim_fields{nn}) = xnim.(nim_fields{nn});
			end
		end

		% Add additional subunit to NIM for each 2D subunit
		N2d = length(xnim.twoD_subunits);
		if (N2d > 0) && ~isempty(Xstims)
			% Make X-matrix for 2D-subunits
			Nstims = length(Xstims);

			for nn = 1:N2d
				curr_sub = xnim.twoD_subunits(nn);
            
				% Create new Xstim matrix for linearized tent basis functions
				gs = curr_sub.apply_filters( Xstims );
				Xstims{Nstims+nn} = curr_sub.make_NL2d_Xmat( gs );

				% Augment NIM stim params
				stim_params = NIM.create_stim_params(...
							[length(curr_sub.ticks{1}), length(curr_sub.ticks{2}), 1], ...    
								'boundary_conds', [Inf, Inf, Inf]);
                          
						nim.stim_params(Nstims+nn) = stim_params;
            
				% Add subunit to NIM
				nim = nim.add_subunits({'lin'}, 1, ...
								'xtargs', Nstims+nn, ...
								'init_filts', {curr_sub.NL2d(:)} );
				% Convert regularization
				nim.subunits(end).reg_lambdas = curr_sub.reg_lambdas;
				nim.subunits(end).reg_lambdas.d2xt = curr_sub.reg_lambdas.nld2;
				% Zero out other xnim regularization
				nim.subunits(end).reg_lambdas.d2x = 0;
				nim.subunits(end).reg_lambdas.d2t = 0;  
				nim.subunits(end).reg_lambdas.nld2 = 0;
			end
		end

	end

% OVERLOADED FUNCTIONS FOR twoD_subunits
	function lambdas = get_reg_lambdas( xnim, varargin )
	% Usage: lambdas = xnim.get_reg_lambdas( varargin )
	% Overloaded for XNIM: Gets regularization lambda values of specified types (see varargin)
	% from a set of xnim 2D-subunits corresponding to each filter (i.e., 2 per subunit)
	%
	% INPUTS:
	%   optional flags:
	%      'filters_to_fit': vector specifying which filters to extract lambda values from (2 per subunit: number accordingly)
	%      'lambda_type': string specifying the regularization type
	% OUTPUTS:
	%   lambdas: [K,N] matrix of lambda values, K is the number of specified lambda_types and N is the number of subunits

		filters_to_fit = 1:2*length(xnim.twoD_subunits); % default is to grab reg values from all subunits

		jj = 1;  reg_types = {};
		while jj <= length(varargin)
			switch lower(varargin{jj})
				case 'fit_filters'
					filters_to_fit = varargin{jj+1};
					assert(all(ismember(filters_to_fit,1:2*length(xnim.twoD_subunits))),'invalid target filters specified');
					jj = jj + 2;
				case xnim.allowed_reg_types
					reg_types = cat(1,reg_types,lower(varargin{jj}));
					jj = jj + 1;
				otherwise
					error('Invalid input flag'); 
			end		
		end
		
		lambdas = nan(length(reg_types),length(filters_to_fit));
		if isempty(reg_types)
			warning( 'No regularization type specified, returning nothing' );
		end
		
		for ii = 1:length(reg_types)
			for jj = 1:length(filters_to_fit)
				lambda_vals = xnim.twoD_subunits(ceil(filters_to_fit(jj)/2)).reg_lambdas.(reg_types{ii});
				% Check to see if one or two lambda values specified
				if length(lambda_vals) == 1
					lambdas(ii,jj) = lambda_vals;
				else  % potentially different regularization for each filter
					lambdas(ii,jj) = lambda_vals( mod(filters_to_fit(jj)-1,2)+1 );
				end
			end
		end
	end
	
	function Tmats = make_Tikhonov_matrices( xnim )  % Overload for twoD-subunits
	% Usage: Tmats = xnim.make_Tikhonov_matrices()
	% Creates a struct containing the Tikhonov regularization matrices, given the stimulus and regularization 
	% parameters specified in the nim
    
		Nstims = length(xnim.stim_params); % number of unique stimuli 
		Xtargs = [xnim.twoD_subunits(:).Xtargs];

		deriv_reg_types = xnim.allowed_reg_types(strncmp(xnim.allowed_reg_types,'d',1)); % set of regularization types where we need a Tikhonov matrix
		cnt = 1;
		Tmats = [];
		for ii = 1:Nstims % for each stimulus
			cur_filters = find(Xtargs == ii); % get set of subunits acting on this stimulus (assume filters act on same Xtarg?
			for jj = 1:length(deriv_reg_types) % check each possible derivative regularization type
				cur_lambdas = xnim.get_reg_lambdas(deriv_reg_types{jj},'fit_filters',cur_filters);
				if any(cur_lambdas > 0)
					cur_Tmat = NIM.create_Tikhonov_matrix(xnim.stim_params(ii),deriv_reg_types{jj});
					Tmats(cnt).Tmat = cur_Tmat;
					Tmats(cnt).Xtarg = ii;
					Tmats(cnt).type = deriv_reg_types{jj};
					cnt = cnt + 1;
				end
			end
		end      
	end
	
end

end