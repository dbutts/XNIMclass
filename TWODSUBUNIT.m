classdef TWODSUBUNIT

properties
	Xtarg;
	ks; 
	xticks;
	yticks;
	NL2d;
end

properties (Hidden)
	TB;  % tent-basis object, defined on initialization
end

%% --------------------- Constructor --------------------- %%
methods
	
	function sub2d = TWODSUBUNIT( ks, stims, varargin )
	% Usage: sub2d = TWODSUBUNIT( ks, stims, varargin )
	%
	% INPUTS:
	%   ks: cell list with 2 elements -- each filter. If has one element, the delay by 3 the first filter?
	%   stims: Xstim or stims -- whatever is being used
	%   optional arguments:
	%     'Xtarg': X-target for subunit. Default will be 1
	%     'Nbins': number of bins to represent given dimention (NL will be NBxNB)
		
		defaults.Xtarg = 1;
		defaults.Nbins = 20;
		[~,parsed_options] = NIM.parse_varargin( varargin, [], defaults );
		
		sub2d.Xtarg = parsed_options.Xtarg;
		sub2d.ks = ks;
		if length(ks) < 2
			sub2d.ks{2} = sub2d.ks{1};
			sub2d.ks{2}(3:end) = sub2d.ks{1}(1:end-2);
			sub2d.ks{2}(1:2) = 0;
		end
		
		[sub2d,mean_gs] = sub2d.normalize_filters( stims );
		
		% Establish scale of filters (and NL-resolution)
		% will go to plus/minus 3-stds
		dx = 6.0/parsed_options.Nbins;
		sub2d.xticks = -3:dx:3;
		sub2d.yticks = -3:dx:3;
		
		sub2d.NL2d = zeros( parsed_options.Nbins, parsed_options.Nbins );
		sub2d.TB = TentBasis2D( sub2d.xticks, sub2d.yticks );
	end
end

%% ------------- CALCULATIONS ------------------
methods
	function gs = apply_filters( sub2d, stims )
	% Usage: gs = sub2d.apply_filters( stims )
	
		NT = size(stims{1},1);
		gs = zeros(NT,2);
		gs(:,1) = stims{sub2d.Xtarg}*sub2d.ks{1};
		gs(:,2) = stims{sub2d.Xtarg}*sub2d.ks{2};
	end
		
	function sub_out = apply_NL( sub2d, gs )
	% Usage: sub_out = sub2d.apply_NL( gen_signals )
	
		Xmat = sub2d.TB.InputNL2D( gs ); % produce X-matrix
		%sub_out = Xmat*sub2d.TB.tentcoeff(:);
		sub_out = Xmat*sub2d.NL2d(:);  % same thing but not private
	end
	
	function sub_out = process_stim( sub2d, stims )
	% Usage: sub_out = sub2d.process_stim( stims )
	
		sub_out = sub2d.apply_NL( sub2d.apply_filter( stims ));
	end
	
	function [] = display_filter( sub2d, dims, varargin )
	% Usage: [] = display_filter( sub2d, dims, varargin )
	
		assert((nargin > 1) && ~isempty(dims), 'Must enter filter dimensions.' )

		[plotloc,parsed_options,modvarargin] = NIM.parse_varargin( varargin, {'notitle','xt-spatial','xt-separable'} );
		if isempty(plotloc)
			plotloc = [1 2 1];
		end
		assert( plotloc(3) < prod(plotloc(1:2)), 'Invalid plot location.' )
		
		subplot( plotloc(1),plotloc(2),plotloc(3) );
		imagesc(reshape(sub2d.ks{1},dims(1:2)));
		
		subplot( plotloc(1),plotloc(2),plotloc(3)+1 );
		imagesc(reshape(sub2d.ks{2},dims(1:2)));
		
	end
	
	function [] = display_NL( sub2d, varargin )
	% Usage: [] = display_NL( sub2d, varargin )
		imagesc(sub2d.NL2d)
	end
	
	
end
methods (Hidden)
	
	function [sub2d_out,mean_gs] = normalize_filters( sub2d, stims );
	% Usage: [sub2d_out,mean_gs] = sub2d.normalize_filters( stims );
	
		sub2d_out = sub2d;
		gs = sub2d.apply_filters( stims );
		mean_gs = mean(gs);
		for nn = 1:2
			nrm = std(gs(:,nn));
			if nrm > 0
				sub2d_out.ks{nn} = sub2d.ks{nn}/nrm;
			end
		end
	end
	
end
end % (classdef)