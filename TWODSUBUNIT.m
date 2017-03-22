classdef TWODSUBUNIT

properties
    NL2d;           % 2D nonlinearity coefficients
    ks;             % 2x1 cell array of filter coefficients (x-dir/y-dir)
    Xtargs;         % 2x1 array of indices of stimulus the filters act on
    ticks;          % 2x1 cell array of tick values along x/y-dir
    reg_lambdas;    % struct of regularization hyperparameters
    Ksign_con;      % 2x1 cell array defining constraint matrices A and B <- overly complicated to define here. Why....
                    % such that A*TWODSUBUNIT.NL2d(:) < B
    pre_scale;      % 2x2 scale matrix
end

%% *************************** constructor ********************************
methods

	function sub2d = TWODSUBUNIT( ks, Xstims, varargin )
	% Usage: sub2d = TWODSUBUNIT(ks, Xstims, varargin)
	%
	% INPUTS:
	%   ks:     2x1 cell array of filter coefficients. If has one element, 
	%           then second filter is delayed (3 lags) version of the first
	%   Xstims: cell array of stimulus matrices
	%
	%   optional arguments:
	%     'Xtargs': 2x1 array of X-target indices for subunit. Default will be 1 for both
	%     'nld2': regularization for 2-d nonlinearity
	%     'Nticks': 2x1 cell array or scalar defining number of ticks/coeffs to represent each dimension.
	%     'ticks':  2x1 cell array of tick values for each dimension

		% Parse input options
		defaults.Xtargs = [1 1];
		defaults.Nticks = {21,21}; % 20 bins
		defaults.nld2 = 0;
		[~, parsed_options] = NIM.parse_varargin(varargin, [], defaults);
    
		% Get proper number of ticks
		if ~iscell(parsed_options.Nticks) && isscalar(parsed_options.Nticks)
			% duplicate scalar values if supplied
			Nticks{1} = parsed_options.Nticks;
			Nticks{2} = parsed_options.Nticks;
		else
			assert(iscell(parsed_options.Nticks), 'Nticks must be a 2x1 cell array or scalar')
			Nticks = parsed_options.Nticks;
		end
    
		% Set XNIM properties
		sub2d.NL2d = zeros(Nticks{1}, Nticks{2});
    
		sub2d.ks = ks;
		if length(ks) < 2
			sub2d.ks{2} = sub2d.ks{1};
			sub2d.ks{2}(3:end) = sub2d.ks{1}(1:end-2);
			sub2d.ks{2}(1:2) = 0;
		end
    
		if length(parsed_options.Xtargs) == 1
			% duplicate scalar values if supplied
			sub2d.Xtargs = parsed_options.Xtargs*[1 1];
		else
			assert(length(parsed_options.Xtargs) == 2, 'Xtargs must be a 2x1 array or scalar')
			sub2d.Xtargs = parsed_options.Xtargs;
		end
    
		% Establish scale of filters and NL-resolution
		sub2d = sub2d.normalize_filters( Xstims );
		sub2d = sub2d.scale_NLaxes( Xstims );
    
		sub2d.reg_lambdas = TWODSUBUNIT.init_reg_lambdas(); % set all to zero
		sub2d.reg_lambdas.nld2 = parsed_options.nld2;
		%subunit.Ksign_con = {[], []};
		sub2d.Ksign_con = [0 0];
		sub2d.pre_scale = eye(2);
    
	end % method
    
end
%% *************************** output methods *****************************
methods
    
	function gs = apply_filters(subunit, Xstims)
	% Usage: gs = subunit.apply_filters(Xstims)

		if ~iscell(Xstims)
			tmp = Xstims; clear Xstims
			Xstims{1} = tmp;
		end
		
		NT = size(Xstims{1}, 1);
		gs = zeros(NT,2);
		gs(:,1) = Xstims{subunit.Xtargs(1)} * subunit.ks{1};
		gs(:,2) = Xstims{subunit.Xtargs(2)} * subunit.ks{2};
    
	end % method

	function sub_output = apply_NL(subunit, gen_signals)
	% Usage: sub_output = subunit.apply_NL(gen_signals)

		Xmat = subunit.make_NL2d_Xmat(gen_signals); % produce X-matrix
		sub_output = Xmat*subunit.NL2d(:);
	end

	function [dxdy_out, grads] = apply_NL_deriv( twoD_subunit, gen_signals )   % overload
		[dxdy_out, grad_x, grad_y] = apply_NL_deriv2D( twoD_subunit, gen_signals );
		grads = [grad_x; grad_y;];
	end
	
	function [dxdy_out, grad_x, grad_y] = apply_NL_deriv2D( subunit, gen_signals )
	% Usage: [dx_out, dy_out, grad] = subunit.apply_NL_deriv( gen_signals )
	% Calculates the piecewise-constant gradient of the 2D nonlinearity, and evaluates the gradient for each term in gen_signals
	% This is Yuwei's function, and not properly overloading regular subunit function (see above)
	%
	% INPUTS:
	%   gen_signals: Tx2 matrix of filtered stimuli
	%
	% OUTPUTS:
	%   dxdy_out: Tx2 vector containing [df/dx, df/dy] for each time point
	%   grad_x:   Nbins_x x Nbins_y x 2 tensor containing gradient with
	%             respect to x for lower (:,:,1) and upper triangles
	%   grad_y:   Nbins_x x Nbins_y x 2 tensor containing gradient with
	%             respect to y for lower (:,:,1) and upper triangles
    
		% calculate grad
		NLx = subunit.ticks{1};
		NLy = subunit.ticks{2};
		Nx = length(NLx);
		Ny = length(NLy);
		grad_x = zeros(Nx-1,Ny-1,2);
		grad_y = zeros(Nx-1,Ny-1,2);
    
		for x=1:Nx-1
			for y=1:Ny-1
				% vertices of lower triangle
				LT = [NLx(x) NLy(y); NLx(x+1) NLy(y); NLx(x) NLy(y+1)];
				LTf = [subunit.NL2d(x,y) subunit.NL2d(x+1,y) subunit.NL2d(x,y+1)];
				% vertices of upper triangle
				UT = [NLx(x+1) NLy(y); NLx(x+1) NLy(y+1); NLx(x) NLy(y+1)];
				UTf = [subunit.NL2d(x+1,y) subunit.NL2d(x+1,y+1) subunit.NL2d(x,y+1)];
				% calculate gradients
				[grad_x(x,y,1), grad_y(x,y,1)] = TWODSUBUNIT.triangle_grad(LT, LTf);
				[grad_x(x,y,2), grad_y(x,y,2)] = TWODSUBUNIT.triangle_grad(UT, UTf);
			end
		end
    
		% Calculate value of grad for each time point
		[NT, Ndim] = size(gen_signals);
		assert(Ndim == 2, 'gen_signals must be a 2D signal')

		% gen_signals = gen_signals*subunit.prescale;
		gen_signals = subunit.truncate_sig(gen_signals);

		% bin centers
		[~, binx] = histc(gen_signals(:,1), NLx);
		[~, biny] = histc(gen_signals(:,2), NLy);

		dxdy_out = zeros(NT,2);

		% go through all bins
		for x = 1:Nx-1
			for y = 1:Ny-1                    

				% data points in the current bin
				inds = find(binx(:) == x & biny(:) == y);
				if isempty(inds)
					continue;
				end

				% vertices of lower triangle
				LT = [NLx(x) NLy(y); NLx(x+1) NLy(y); NLx(x) NLy(y+1)];
				% vertices of upper triangle
				% UT = [NLx(x+1) NLy(y); NLx(x+1) NLy(y+1); NLx(x) NLy(y+1)];

				% separate points in this bin into two different triangles
				flag = TWODSUBUNIT.point_in_triangle(gen_signals(inds,:), ...
								LT(1,:), LT(2,:), LT(3,:));
				indLT = inds(flag==1);
				indUT = inds(flag==0);

				% assign derivatives to points in lower triangle
				if ~isempty(indLT)
					dxdy_out(indLT,1) = grad_x(x,y,1);
					dxdy_out(indLT,2) = grad_y(x,y,1);
				end
				% assign derivatives to points in upper triangle
				if ~isempty(indUT)
					dxdy_out(indUT,1) = grad_x(x,y,2);
					dxdy_out(indUT,2) = grad_y(x,y,2);
				end
			end
		end
		% dxdy_out = dxdy_out*(TB.prescale');

	end
    
	function sub_output = process_stim(subunit, Xstims)
	% Usage: sub_output = subunit.process_stim( Xstims )
	% Generate output of twoD-subunit

		sub_output = subunit.apply_NL(subunit.apply_filters(Xstims));
    
	end % method

	function [subunit, mean_gs] = normalize_filters(subunit, stims)
	% Usage: [subunit, mean_gs] = subunit.normalize_filters(stims);

		gs = subunit.apply_filters( stims );
		mean_gs = mean(gs);
		for nn = 1:2
			nrm = std(gs(:,nn));
			if nrm > 0
				subunit.ks{nn} = subunit.ks{nn}/nrm;
			end
		end
    
	end % method

end

%% *************************** display methods ****************************
methods
    
	function [] = display_filter(subunit, dims, varargin)
	% Usage: [] = subunit.display_filter(dims, varargin)
    
		assert((nargin > 1) && ~isempty(dims), 'Must enter filter dimensions.' )

		[plotloc, parsed_options, modvarargin] = NIM.parse_varargin(varargin, {'notitle','xt-spatial','xt-separable'});
		if isempty(plotloc)
			plotloc = [1 2 1];
		end
		assert( plotloc(3) < prod(plotloc(1:2)), 'Invalid plot location.')

		subplot(plotloc(1), plotloc(2), plotloc(3));
		if (dims{1}(1) > 1) && (dims{1}(2) > 1)
			filt = reshape(subunit.ks{1}, dims{1}(1:2));
			M = max(abs(filt(:)));
			colormap gray
			imagesc(filt, [-M, M]);
		else
			plot(subunit.ks{1});
		end

		subplot(plotloc(1), plotloc(2), plotloc(3)+1);
		if (dims{2}(1) > 1) && (dims{2}(2) > 1)
			imagesc(reshape(subunit.ks{2}, dims{2}(1:2)));
			colormap gray
		else
			plot(subunit.ks{2});
		end
			
	end % method

	function [] = display_NL(subunit, varargin)
	% Usage: [] = display_NL(subunit, varargin)
    
		M = max(abs(subunit.NL2d(:)));
		imagesc(subunit.ticks{1}, subunit.ticks{2}, subunit.NL2d', [-M, M])
		xlabel('X-dir')
		ylabel('Y-dir')
		set(gca, 'Ydir', 'normal')
		colorbar;
		%colormap(jet);
    
	end % method

end
%% *************************** hidden methods *****************************
methods (Hidden)

	function subunit = init_NL2d(subunit, interaction)
	% Usage: subunit = subunit.init_NL2d(interaction)
	% Initializes values for NL2d
	%
	% INPUTS:
	%   interaction:    'add' | 'mult'
	%
	% OUTPUTS:
	%   subunit:    updated subunit object
    
		if strcmp(interaction, 'add')
			subunit.NL2d = bsxfun(@plus, subunit.ticks{1}(:), subunit.ticks{2}(:)');
		elseif strcmp(interaction, 'mult')
			subunit.NL2d = subunit.ticks{1}(:) * subunit.ticks{2}(:)';
		else
			error('Invalid interaction string specified')
		end
    
	end % method
    
	function twoDsub_out = scale_NLaxes( twoDsub, Xstims )
	% Usage: twoDsub_out = twoDsub.scale_NLaxes( Xstims )
	
		% will go to plus/minus 3-stds
		%subunit = subunit.normalize_filters( Xstims );
		%dx = 6.0/(Nticks{1}-1);
		%dy = 6.0/(Nticks{2}-1);
		%subunit.ticks{1} = -3:dx:3;
		%subunit.ticks{2} = -3:dy:3;

		% will go between 1st and 99th percentile - more appropriate for 2p data
		twoDsub_out = twoDsub.normalize_filters( Xstims );
		gs = twoDsub_out.apply_filters( Xstims );
		prctiles = prctile( gs(:,1), [1,99] );
		twoDsub_out.ticks{1} = linspace( prctiles(1), prctiles(2), size(twoDsub.NL2d,1) );
		prctiles = prctile( gs(:,2), [1,99] );
		twoDsub_out.ticks{2} = linspace( prctiles(1), prctiles(2), size(twoDsub.NL2d,2) );
		
	end % method
	
	function [NLoutput, Counts, binx, biny] = make_NL2d_Xmat(subunit, gen_signals)
	% calculates output of 2D nonlinearity with respect to subunit.NL2d coefficients.
	%
	% Input:
	%       gen_signals:    NTx2 matrix of generating signals
	%
	% Output:
	%       NLoutput:       NTxNpar matrix, where Npar = Nbins_x*Nbins_y
	%       Counts:         ticks{1} x ticks{2} matrix, number of samples
	%                       per vertex
	%
	% Yuwei Cui, Created by Oct 20, 2012
	% Oct 28, 2012 Extend Pixel filter to 2D tent filter

		[NT, Ndim] = size(gen_signals);
		assert(Ndim == 2, 'gen_signals must be a 2D signal')

		% gen_signals = gen_signals*subnit.prescale;
		gen_signals = subunit.truncate_sig(gen_signals);
    
		NLx = subunit.ticks{1};
		NLy = subunit.ticks{2};
    
		% number of sample pts along x/y dimensions
		NNx = length(NLx);
		NNy = length(NLy);

		% bin centers
		[~, binx] = histc(gen_signals(:,1), NLx);
		[~, biny] = histc(gen_signals(:,2), NLy);

		NLoutput = zeros(NT, NNx, NNy);
		if nargout > 1
			Counts = zeros(NNx,NNy);
		end

		% Go through all bins
		for x=1:NNx-1
			for y=1:NNy-1           
            
				% data points in the current bin
				inds = find(binx(:) == x & biny(:) == y);
				if isempty(inds)
					continue;
				end

				% vertices of lower triangle
				LT = [NLx(x) NLy(y); NLx(x+1) NLy(y); NLx(x) NLy(y+1)];
				% vertices of upper triangle
				UT = [NLx(x+1) NLy(y); NLx(x+1) NLy(y+1); NLx(x) NLy(y+1)];
            
				% separate points in this bin into two different triangles
				flag = TWODSUBUNIT.point_in_triangle(gen_signals(inds,:), LT(1,:), LT(2,:), LT(3,:));
				indLT = inds(flag == 1);
				indUT = inds(flag == 0);

				% calculate barycentric coordinates in lower triangle
				lambdasLT = TWODSUBUNIT.barycentric(gen_signals(indLT,:), LT(1,:), LT(2,:), LT(3,:));
				% convert to coefficient loadings
				NLoutput(indLT,x,y)   = lambdasLT(:,1);
				NLoutput(indLT,x+1,y) = lambdasLT(:,2);
				NLoutput(indLT,x,y+1) = lambdasLT(:,3);

				% calculate barycentric coordinates in lower triangle
				lambdasUT = TWODSUBUNIT.barycentric(gen_signals(indUT,:), UT(1,:), UT(2,:), UT(3,:));
				% convert to coefficient loadings
				NLoutput(indUT,x+1,y)   = lambdasUT(:,1);
				NLoutput(indUT,x+1,y+1) = lambdasUT(:,2);
				NLoutput(indUT,x,y+1)   = lambdasUT(:,3);

				if nargout>1
					% count sample per vertex
					NLT = length(indLT); 
					Counts(x,y)   = Counts(x,y) + NLT;
					Counts(x+1,y) = Counts(x+1,y) + NLT;
					Counts(x,y+1) = Counts(x,y+1) + NLT;

					NUT = length(indUT);
					Counts(x+1,y)   = Counts(x+1,y) + NUT;
					Counts(x+1,y+1) = Counts(x+1,y+1) + NUT;
					Counts(x,y+1)   = Counts(x,y+1) + NUT;
				end
			end
		end
    
		NLoutput = reshape(NLoutput, NT, NNx*NNy);
    
	end % method
        
	function [sig, trucated_indxs] = truncate_sig(subunit, sig)
	% truncate stimulus at boundaries of nonlinearity            

		NLx = subunit.ticks{1};
		NLy = subunit.ticks{2};
    
		if nargout >= 2
			trucated_indxs = find(sig(:,1) <= NLx(1) | sig(:,2) <= NLy(1) | sig(:,1) >= NLx(end) | sig(:,2) >= NLy(end));
		end

		sig(sig(:,1) <= NLx(1),1) = NLx(1);
		sig(sig(:,2) <= NLy(1),2) = NLy(1);
		sig(sig(:,1) >= NLx(end),1) = NLx(end)-1e-7; % move off boundary
		sig(sig(:,2) >= NLy(end),2) = NLy(end)-1e-7; % move off boundary

	end % method

end

%% *************************** static methods *****************************
methods (Static)
    
	function reg_lambdas = init_reg_lambdas()
	% Usage: reg_lambdas = init_reg_lambdas()
	% Creates reg_lambdas struct and sets default values to 0; this method
	% comes from the SUBUNIT class and should be kept up to date with it,
	% else fit_NL2d won't work.
     
		reg_lambdas.nld2 = 0; %second derivative of tent basis coefs
		reg_lambdas.d2xt = 0; %spatiotemporal laplacian
		reg_lambdas.d2x = 0; %2nd spatial deriv
		reg_lambdas.d2t = 0; %2nd temporal deriv
		reg_lambdas.l2 = 0; %L2 on filter coefs
		reg_lambdas.l1 = 0; %L1 on filter coefs
	end
    
	function lambda = barycentric(P, A, B, C)
	% Usage: lambda = TWODSUBUNIT.barycentric(P, A, B, C)
	%
	% Calculate barycenteric coordinate of point P with respect to triangle ABC
	% Yuwei Cui, Oct 30, 2012

		Npt = size(P,1);
		lambda = zeros(Npt,3);

		x = [A(1) B(1) C(1)];
		y = [A(2) B(2) C(2)];

		M = ((y(2)-y(3))*(x(1)-x(3))+(x(3)-x(2))*(y(1)-y(3)));
		lambda(:,1) = ((y(2)-y(3))*(P(:,1)-x(3))+(x(3)-x(2))*(P(:,2)-y(3)))/M;
		lambda(:,2) = ((y(3)-y(1))*(P(:,1)-x(3))+(x(1)-x(3))*(P(:,2)-y(3)))/M;

		% lambda(:,1) = ((B(2)-C(2))*(P(:,1)-C(1))+(C(1)-B(1))*(P(:,2)-C(2)))./...
		%     ((B(2)-C(2))*(A(1)-C(1))+(C(1)-B(1))*(A(2)-C(2)));
		% lambda(:,2) = ((C(2)-A(2))*(P(:,1)-C(1))+(A(1)-C(1))*(P(:,2)-C(2)))./...
		%     ((B(2)-C(2))*(A(1)-C(1))+(C(1)-B(1))*(A(2)-C(2)));

		lambda(:,3) = ones(Npt,1)-lambda(:,1)-lambda(:,2);
    
	end % method
    
	function flag = point_in_triangle(P, A, B, C) 
	% Usage: flag = TWODSUBUNIT.point_in_triangle(P, A, B, C) 
	%
	% This function will decide whether point P is inside the triangle 
	% defined by the points A, B, C or not (edge included). The Input P can
	% also be a Nx2 vector, each row is a point
	%
	% This method is sometimes referred as the "Barycentric Technique"
	% For more info, see the following link:
	% http://en.wikipedia.org/wiki/Barycentric_coordinates_%28mathematics%29
	% or the one I used:
	% http://www.blackpawn.com/texts/pointinpoly/default.html
	% Implemented by Yuwei Cui, Oct 30. 2012

		A = A(:)';
		B = B(:)';
		C = C(:)';
		Npt = size(P,1);

		% Compute vectors        
		v0 = C - A;
		v1 = B - A;
		v2 = P - ones(Npt,1)*A;

		% Compute dot products
		dot00 = v0*v0';
		dot01 = v0*v1';
		dot11 = v1*v1';
		dot02 = v2*v0';
		dot12 = v2*v1';

		% Compute barycentric coordinates
		invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
		u = (dot11 * dot02 - dot01 * dot12) * invDenom;
		v = (dot00 * dot12 - dot01 * dot02) * invDenom;

		% Check if point is in triangle
		flag = (u >= 0) & (v >= 0) & (u + v <= 1);
	end
    
	function [dx, dy] = triangle_grad(xy, z)
	% Usage: [dX, dY] = TWODSUBUNIT.triangle_grad(T, f)
	% Calculates the derivate along the x and y directions wrt z of an 
	% arbitrarily-oriented triangle in 3D space. The formulas here can be
	% derived by considering the 3 (x,y,z) coordinates as defining a plane,
	% and calculating the equation of that plane as z = a*x + b*y + c; then
	% dz/dx = a and dz/dy = b (or dx = a and dy = b, as per the notation of
	% this function).
	%
	% INPUTS:
	%   xy: 3x2 matrix of x-y coordinates of triangle vertices
	%   z:  3x1 vector of z-coordinates of triangle vertices
	%
	% OUTPUTS:
	%   dx: derivative in x-direction
	%   dy: derivative in y-direction
    
		x = xy(:,1);
		y = xy(:,2);

		M = ((y(2)-y(3))*(x(1)-x(3))+(x(3)-x(2))*(y(1)-y(3)));

		dx = (z(1)*(y(2)-y(3)) + z(2)*(y(3)-y(1)) + z(3)*(-(y(2)-y(3))-(y(3)-y(1))) ) / M;

		dy = (z(1)*(x(3)-x(2)) + z(2)*(x(1)-x(3))+ z(3)*(-(x(3)-x(2))-(x(1)-x(3))) ) / M;
    
	end
    
end

end