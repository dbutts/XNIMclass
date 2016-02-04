classdef XNIM < NIM
	
properties
	TwoDsubunits;
	
	
end

% Methods defined in other files
methods (Static)
	flag = PointInTriangle( P, A,B,C ); 
	lambda = BaryCentric( P, A,B,C ); % calculates Barycenteric Coordinate of point P with respect to triangle ABC
end
	
methods
	
	function xnim = XNIM( nim )
	% Usage: xnim = XNIM( nim )
	
		nimfields = fields(nim);
		for nn = 1:length(nimfields)
			xnim.(nimfields{nn}) = nim.(nimfields{nn});
		end
		xnim.TwoDsubunits = [];
	end

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
	end
	
	function xnim = fit_2dfilters( xnim, Robs, Xstims, varargin )
	end
	
	function xnim = fit_2dNL( xnim, Robs, Xstims, varargin )
	% Usage: xnim = xnim.fit_2dNL( Robs, Xstims, Uindx, varargin )
	
		defaults.subs = 1:length(xnim.TwoDsubunits);
		[Uindx,parsed_options] = NIM.parse_varargin( varargin, [], defaults );
		
		% Calculate rest of model (regular subunits and non-target 2d's)
		[~,~,mod_internals] = xnim.eval_model( [], Xstims, Uindx );
		Gregsubs = mod_internals.G;
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
			
	end

methods (Hidden)
	
	function [nim,Xstims] = convert2NIM( xnim, Xstims )
	% Usage: nim = xnim.convert2NIM( <Xstims> )
	% Allows NIM function calls (that are not overloaded by XNIM
	
		if nargin < 2
			Xstims = [];
		end
		
		nim = NIM();
		nimfields = fields(xnim);
		for nn = 1:length(nimfields)
			if ~strcmp(nimfields{nn},'TwoDsubunits')
				nim.(nimfields{nn}) = xnim.(nimfields{nn});
			end
		end
		
		% Add additional subunit to NIM to capture 2Dsubunit contributions
		N2d = length(xnim.TwoDsubunits);
		if (N2d > 0) && ~isempty(Xstims)
			% Make X-matrix for 2D-subunits
			Nstims = length(Xstims);
			NT = size(Xstims{1},1);
			gain_stim_param = nimtmp.stim_params(1);
			gain_stim_params.dims = [1 1 1];

			for nn = 1:N2d
				Xstims{Nstims+nn} = xnim.Two2subunits(nn).process_stim( Xstims );
				nim.stim_params(Nstims+nn) = gain_stim_params;
				% Add gain-modules
				nim = nim.add_subunits( {'lin'}, 1, 'xtargs',Nstims+nn, 'init_filts', {[1]} );
			end
		end
		
	end
end
end