classdef XNIM < NIM
% Class implementation of NIM that incorporates 2D non-parametric
% nonlinearities for combining two signals.

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

    nim_fields = fields(nim);
    for nn = 1:length(nim_fields)
        xnim.(nim_fields{nn}) = nim.(nim_fields{nn});
    end
    xnim.twoD_subunits = [];
    
    % change superclass properties
    xnim.allowed_reg_types = {'d2xy', 'l2'};
    xnim.version = '0.0';
    
    end % method
    
end

%% *************************** fitting methods ****************************
methods 

    function xnim_out = fit_filters( xnim, Robs, Xstims, varargin )
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

    function xnim = fit_2dfilters( xnim, Robs, Xstims, varargin )
    
    end  % method

    function xnim = fit_2dNL( xnim, Robs, Xstims, varargin )
    % Usage: xnim = xnim.fit_2dNL( Robs, Xstims, Uindx, varargin )

    defaults.subs = 1:length(xnim.TwoDsubunits);
    [Uindx,parsed_options] = NIM.parse_varargin( varargin, [], defaults );

    % Calculate rest of model (regular subunits and non-target 2d's)
    [~,~,mod_internals] = xnim.eval_model( [], Xstims, Uindx );
    Gregsubs = mod_internals.G;
    
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
    xnim.twoD_subunits(end) = xnim.twoD_subunits(end).init_NL2D(parsed_options.interaction);
        
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
    %       'Nbins':        2x1 cell array or scalar defining number of 
    %                       bins to represent each dimension.
    %
    % OUTPUTS:
    %   xnim:   updated XNIM object
    
    % make sure there are enough subunits
    if length(xnim.subunits) < 2
        error('Not enough subunits to convert to TWODSUBNIT')
    end
    
    defaults.interaction = 'add';
    fields_to_remove = 'interaction';
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
    
    end
    
    function [LL, pred_rate, mod_internals, LL_data] = eval_model( xnim, Robs, Xstims, varargin )
	% Usage: [LL, pred_rate, mod_internals, LL_data] = xnim.eval_model( Robs, Xstims, <eval_inds>, varargin ) 
	
    nimtmp = xnim.convert2NIM( Xstims );
    [LL,pred_rate,mod_internals,LL_data] = nimtmp.eval_model( Robs, Xstims, varargin{:} );

    if length(xnim.TwoDsubunits) > 0
        mod_internals.fgint2d = mod_internals.fgint(end);
        mod_internals.gint = mod_internals.gint(:,1:end-1);
        mod_internals.fgint = mod_internals.fgint(:,1:end-1);
    else
        mod_internals.fgint2d = [];
    end
			
    end % method
    
end
%% *************************** hidden methods *****************************
methods (Hidden)

    function [nim,Xstims] = convert2NIM( xnim, Xstims )
    % Usage: nim = xnim.convert2NIM( <Xstims> )
    % Allows NIM function calls (that are not overloaded by XNIM)

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

    % Add additional subunit to NIM to capture 2Dsubunit contributions
    N2d = length(xnim.twoD_subunits);
    if (N2d > 0) && ~isempty(Xstims)
        % Make X-matrix for 2D-subunits
        Nstims = length(Xstims);
        NT = size(Xstims{1},1);
        gain_stim_param = nim.stim_params(1);
        gain_stim_params.dims = [1 1 1];

        for nn = 1:N2d
            Xstims{Nstims+nn} = xnim.twoD_subunits(nn).process_stim( Xstims );
            nim.stim_params(Nstims+nn) = gain_stim_params;
            % Add gain-modules
            nim = nim.add_subunits({'lin'}, 1, 'xtargs', Nstims+nn, 'init_filts', {[1]});
        end
    end

    end
    
end

end