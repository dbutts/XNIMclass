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

%% *************************** constructor ********************************
methods

    function xnim = XNIM(nim)
    % Usage: xnim = XNIM(nim)

    if nargin == 0
        % handle the no-input-argument case by returning a null model. This
        % is important when initializing arrays of objects  
        return 
    end
    
    nim_fields = fields(nim);
    for nn = 1:length(nim_fields)
        xnim.(nim_fields{nn}) = nim.(nim_fields{nn});
    end
    xnim.twoD_subunits = [];
    
    % change superclass properties
    xnim.version = '0.0';
    
    end % method
    
end

%% *************************** fitting methods ****************************
methods 

    function xnim_out = fit_filters(xnim, Robs, Xstims, varargin)
    % Usage: xnim = fit_filters( xnim, Robs, Xstims, Uindx, varargin )
    % Overloaded fit_filters that fits all non-2D subunits. It just calls NIM.fit_filters, but with addition Xstim
    % corresponding to contributions of all 2D subunits (for which it fits the gains and translates this to the 
    % output of each NL2d

    Nstims = length(Xstims);
    N2d = length(xnim.TwoDsubunits);

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

    function xnim = fit_2dfilters(xnim, Robs, Xstims, varargin)
    
    end  % method

    function xnim_out = fit_NL2d(xnim, Robs, Xstims, varargin)
    % Usage: xnim = xnim.fit_NL2d(Robs, Xstims, Uindx, varargin)
    % Estimates 2D nonlinearities of the XNIM model.
    %
    % INPUTS:
    %   Robs:           vector of response observations (e.g. spike counts)
    %   Xstims:         cell array of stimuli
    %   <train_inds>:   index values of data on which to fit the model 
    %                   [default to all indices in provided data]
    %
    %   optional flags:
    %       ('subs',fit_subs): set of twoD_subunits whos coefficients we 
    %           want to optimize [default is all]
    %       ('gain_funs',gain_funs): matrix of multiplicative factors, one 
    %           column for each subunit
    %       ('fit_offsets',fit_offsets): vector of bools, (or single bool) 
    %           specifying whether to fit the additive offset terms 
    %           associated with each subunit
    %       ('optim_params',optim_params): struct of desired optimization 
    %           parameters, can also be used to override any of the default
    %           values for other optional inputs
    %       ('silent',silent): boolean variable indicating whether to 
    %           suppress the iterative optimization display
    %       ('fit_spk_hist',fit_spk_hist): boolean indicating whether to 
    %           hold the spk NL filter constant
    %
    % OUTPUTS:
    %   xnim: new xnim object with optimized 2d nonlinearities
    
    defaults.subs = 1:length(xnim.twoD_subunits);
    [~, parsed_varargin, modvarargin] = NIM.parse_varargin(varargin, {'subs'}, defaults);
    twoD_subs = parsed_varargin.subs; % record for later
    
    % twoD_subunits will be turned into normal NIM subunits and appended to
    % subunits array; update 'subs' field to reflect this
    modvarargin{end+1} = 'subs';
    modvarargin{end+1} = length(xnim.subunits) + parsed_varargin.subs;
    
    % Reformulate XNIM into an NIM for fitting; let fit_filters take care
    % of calculating outputs of non-target subunits
    [nim, Xstims] = xnim.format4NIMfitting(Xstims);
    
    % fit 2d nonlinearity coefficients
    nim = nim.fit_filters(Robs, Xstims, modvarargin{:});
    
    % -------------- translate everything back into XNIM ------------------
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
    % modify fit history
    xnim_out.fit_props = nim.fit_props;
    xnim_out.fit_props.fit_type = 'NL2d';
    xnim_out.fit_history = cat(1, xnim_out.fit_history, nim.fit_history(end));
    xnim_out.fit_history(end).fit_type = 'NL2d';    
    
    end % method

    function xnim_out = reg_path_NL2d(xnim, Robs, Xstims, Uindx, XVindx, varargin)
    % Usage: xnim = xnim.reg_path_NL2d(Robs, Xstims, Uindx, XVindx, varargin)
    % Estimates 2D nonlinearities of the XNIM model while searching for 
    % optimal regularization parameters using XVindx.
    %
    % INPUTS:
    %   Robs:   vector of response observations (e.g. spike counts)
    %   Xstims: cell array of stimuli
    %   Uindx:  index values of data on which to fit the model
    %   XVindx: index values of data on which cross-validation is evaluated
    %
    %   optional flags:
    %       ('subs',fit_subs): set of twoD_subunits whos coefficients we 
    %           want to optimize [default is all]
    %       ('gain_funs',gain_funs): matrix of multiplicative factors, one 
    %           column for each subunit
    %       ('fit_offsets',fit_offsets): vector of bools, (or single bool) 
    %           specifying whether to fit the additive offset terms 
    %           associated with each subunit
    %       ('optim_params',optim_params): struct of desired optimization 
    %           parameters, can also be used to override any of the default
    %           values for other optional inputs
    %       ('silent',silent): boolean variable indicating whether to 
    %           suppress the iterative optimization display
    %       ('fit_spk_hist',fit_spk_hist): boolean indicating whether to 
    %           hold the spk NL filter constant
    %
    % OUTPUTS:
    %   xnim: new xnim object with optimized 2d nonlinearities
    
    defaults.subs = 1:length(xnim.twoD_subunits);
    [~, parsed_varargin, modvarargin] = NIM.parse_varargin(varargin, {'subs'}, defaults);
    twoD_subs = parsed_varargin.subs; % record for later
    
    % twoD_subunits will be turned into normal NIM subunits and appended to
    % subunits array; update 'subs' field to reflect this
    modvarargin{end+1} = 'subs';
    modvarargin{end+1} = length(xnim.subunits) + parsed_varargin.subs;
    
    % Reformulate XNIM into an NIM for fitting; let fit_filters take care
    % of calculating outputs of non-target subunits
    [nim, Xstims] = xnim.format4NIMfitting(Xstims);
    
    % fit 2d nonlinearity coefficients
    nim = nim.reg_path(Robs, Xstims, Uindx, XVindx, modvarargin{:});
    
    % -------------- translate everything back into XNIM ------------------
    xnim_out = xnim;
    for i = 1:length(twoD_subs)
        xnim_sub_indx = twoD_subs(i);
        nim_sub_indx = length(xnim.subunits) + parsed_varargin.subs(i);
        xnim_out.twoD_subunits(xnim_sub_indx).NL2d = reshape(...
            nim.subunits(nim_sub_indx).filtK, ...
            length(xnim_out.twoD_subunits(xnim_sub_indx).ticks{1}), ...
            length(xnim_out.twoD_subunits(xnim_sub_indx).ticks{2}));
        xnim_out.twoD_subunits(xnim_sub_indx).reg_lambdas = ...
            nim.subunits(nim_sub_indx).reg_lambdas;
    end
    xnim_out.spkNL = nim.spkNL;
    xnim_out.spk_hist = nim.spk_hist;
    % modify fit history
    xnim_out.fit_props = nim.fit_props;
    xnim_out.fit_props.fit_type = 'NL2d';
    xnim_out.fit_history = cat(1, xnim_out.fit_history, nim.fit_history(end));
    xnim_out.fit_history(end).fit_type = 'NL2d';    
    
    end % method
    
end
%% *************************** helper methods *****************************
methods
    
    function xnim = add_2d_subunit(xnim, subunit_1, subunit_2, Xstims, varargin)
    % Usage: xnim = xnim.add_2d_subunit(subunit_1, subunit_2, Xstims, varargin)
    % Takes two SUBUNIT inputs (probably from an NIM) and converts them
    % into a TWODSUBUNIT
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
    Xtargs{1} = subunit_1.Xtarg;
    Xtargs{2} = subunit_2.Xtarg;
    
    % initialize 2D subunit
    xnim.twoD_subunits(end+1) = TWODSUBUNIT(ks, Xstims, 'Xtargs', Xtargs, modvarargin{:});
    
    % initialize nonlinearity
    xnim.twoD_subunits(end) = xnim.twoD_subunits(end).init_NL2d(parsed_options.interaction);
        
    end
    
    function xnim = convert_to_2d_subunit(xnim, sub_indx_1, sub_indx_2, Xstims, varargin)
    % Usage: xnim = xnim.add_2d_subunit(sub_indx_1, sub_indx_2, Xstims, varargin)
    % Takes two SUBUNITS from the current XNIM object and converts them
    % into a TWODSUBUNIT, and deletes the original SUBUNITS
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
    Xtargs{1} = xnim.subunits(sub_indx_1).Xtarg;
    Xtargs{2} = xnim.subunits(sub_indx_2).Xtarg;
    
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
	
    [nim_tmp, Xs_mod] = xnim.convert2NIM(Xstims);
    
    [LL, pred_rate, mod_internals, LL_data] = ...
                            nim_tmp.eval_model(Robs, Xs_mod, varargin{:});

    if ~isempty(xnim.twoD_subunits)
        mod_internals.fgint2d = mod_internals.fgint(:,end);
        mod_internals.gint = mod_internals.gint(:,1:end-1);
        mod_internals.fgint = mod_internals.fgint(:,1:end-1);
    else
        mod_internals.fgint2d = [];
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

    % will be spike-history or spkNL plot?
    extra_plots = [(xnim.spk_hist.spkhstlen > 0) ~isempty(mod_outs)]; 

    if sum(extra_plots) == 0
        Nrows = Nmods;
        Ncols = 3;
    else
        % Then need extra column (and possibly extra row)
        Nrows = max([Nmods sum(extra_plots)]);
        Ncols = 4;
    end

    if nargout > 0
        fig_handle = figure;
    else
        figure;
    end
    % Plot Subunit info
    for nn = 1:Nmods
        dims{1} = xnim.stim_params(xnim.twoD_subunits(nn).Xtargs{1}).dims;
        dims{2} = xnim.stim_params(xnim.twoD_subunits(nn).Xtargs{2}).dims;
        xnim.twoD_subunits(nn).display_filter(...
            dims, ...
            [Nrows Ncols (nn-1)*Ncols+1], ...
            modvarargin{:});
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

    function [nim, Xstims] = convert2NIM(xnim, Xstims)
    % Usage: nim = xnim.convert2NIM(<Xstims>)
    % Allows NIM function calls (that are not overloaded by XNIM). Sets
    % output of each twoD_subunit to be a univariate signal that a scale
    % factor can be fit to, if desired.

    if nargin < 2
        Xstims = [];
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

    function [nim, Xstims] = format4NIMfitting(xnim, Xstims)
    % Usage: [nim, Xstims] = xnim.format4NIMfitting(<Xstims>)
    % Turns twoD_subunits into linear subunits of an NIM that can be fit
    % using the standard fit_filters routine in the NIM.

    if nargin < 2
        Xstims = [];
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
            
            % create new Xstim matrix for linearized tent basis functions
            gs = curr_sub.apply_filters(Xstims);
            Xstims{Nstims+nn} = curr_sub.make_NL2d_Xmat(gs);
            
            % augment NIM stim params
            stim_params = NIM.create_stim_params(...
                              [length(curr_sub.ticks{1}), length(curr_sub.ticks{2}), 1], ...
                              'boundary_conds', [Inf, Inf, Inf]);
                          
            nim.stim_params(Nstims+nn) = stim_params;
            
            % Add subunit to NIM
            nim = nim.add_subunits({'lin'}, 1, ...
                                   'xtargs', Nstims+nn, ...
                                   'init_filts', {curr_sub.NL2d(:)});
            nim.subunits(end).reg_lambdas = curr_sub.reg_lambdas;
        end
    end

    end
    
end

end