function xnim_out = twoDreg_path( xnim, Robs, Xs, Uindx, XVindx, varargin )
% Usage: xnim_out = xnim.twoDreg_path( Robs, Xs, Uindx, XVindx, <varargin> )
%
% Rudimentary function to fit filters or non-parametric nonlinearities for range of reg values. This function
% is not necessarily robust or tested outside of particular situations.
%
% optional inputs to varargin:
%   'subs': which 2-D subunits to fit (default all)
%   'filts': which filters within each 2-D subunit to fit (default both)
%   'lambda-ID' refers to which lamdba to do regularization path: default = 'd2t'
%   'silent' = 1 will suppress reg_path output, and 2 will suppress except for reporting maximum

if nargin < 5
	error('Require Uindx and XVindx as arguments.')
end
Nsubs = length(xnim.twoD_subunits);

% Set default options
defaults.subs = 1:Nsubs;
defaults.filts = [1 2];
defaults.silent = 0;
defaults.lambdaID = 'd2t';
defaults.L2s = [];

[~,parsed_options,fit_options] = NIM.parse_varargin( varargin, ...
    {'silent','lambdaID','lambdaid','L2s','l2s','subs','filts','optim_params'}, defaults );
targets = parsed_options.subs;

assert(all(ismember(targets,1:Nsubs)),'invalid target subunits specified');
silent = parsed_options.silent;
L2s = parsed_options.L2s;
lambdaID = parsed_options.lambdaID;
if ~strcmp( lambdaID, 'nld2' )
	filts_to_fit = parsed_options.filts;
else
	filts_to_fit = 1; % only one thing to fit if fitting nonlinearity
end

modvarargin{1} = Uindx;
modvarargin{2} = 'silent'; % to pass into subfunctions
modvarargin{3} = 1;
for nn = 1:length(fit_options)
	modvarargin{3+nn} = fit_options{nn};
end
modvarargin{4+length(fit_options)} = 'subs';
modcounter = length(modvarargin)+1;

if strcmpi(lambdaID,'l1')
	regchar = '1';
else
	regchar = '2';
end

% Add targets to modvarargin
%modvarargin{modcounter} = 'subs';
%modcounter = modcounter + 1;

for tar = targets
	for filtN = filts_to_fit
		if isempty(L2s)  
			%% Do order-of-mag reg first
			%L2s = [0 0.1 1.0 10 100 1000 10000 1e5];
			Niter = 0;  notdone = 1; 
			lambda = 0;
			if ~silent 
				fprintf( 'Order-of-magnitude L%c reg path: target %d, filt %d\n', regchar, tar, filtN )
			end		
			LLregs = []; L2s = [];
			while ((Niter < 10) && notdone)
				Niter = Niter + 1;
				regfit = xnim.set_2Dreg_params( 'subs', tar, 'filts', filtN, lambdaID, lambda );
			
				% Only refit specific target
				modvarargin{modcounter} = tar;
			
				if strcmp( lambdaID, 'nld2' )
					regfit = regfit.fit_NL2d( Robs, Xs, modvarargin{:} );
				else
					regfit = regfit.fit_2Dfilters( Robs, Xs, modvarargin{:} );
				end
			
				fitsaveM{Niter} = regfit;
				[LL,~,~,LLdata] = regfit.eval_model( Robs, Xs, XVindx, fit_options{:} );
				L2s(Niter) = lambda;
				LLregs(Niter) = LL-LLdata.nullLL;
				if ~silent 
					fprintf( '  %8.2f: %f\n', lambda, LLregs(Niter) )
				end
			
%			if Niter > 2
%				if (LLregs(Niter) < LLregs(Niter-1)) && (LLregs(Niter) < LLregs(Niter-2))
%					notdone = 0;
%				end
%			end
			% Check that magnitude of LL is not too big
				if lambda == 0
					lambda = 0.1;
				else
					if (LLdata.filt_pen < (10*(LL-LLdata.nullLL))) || (Niter < 3)
						lambda = lambda*10;
					elseif ((Niter == 2) && (LLdata.filt_pen > abs(LL-LLdata.nullLL)))
						Niter = Niter-1;
						lambda = lambda/10;
					else
						notdone = 0;
					end
				end
			end
	
			[~,bestord] = sort(LLregs);
			% Check to see if contiguous
			if abs(bestord(end)-bestord(end-1)) == 1
				loweredge = min(bestord(end-1:end));
			else  % otherwise go the right direction on the best
				if bestord(end) > bestord(end-1)
					loweredge = bestord(end)-1;
				else
					loweredge = bestord(end);
				end
			end
			mag = L2s(loweredge);
			LLbounds = LLregs(loweredge+[0 1]); 

			% Zoom in on best regularization
			if mag == 0
			 L2s = [0 0.01 0.02 0.04 0.1];
			else
				L2s = mag*[1 2 4 6 8 10];
			end
			Nreg = length(L2s);
			LLregs = zeros(Nreg,1);
			LLregs(1) = LLbounds(1);    LLregs(end) = LLbounds(2);
			fitsave{1} = fitsaveM{loweredge};
			fitsave{Nreg} = fitsaveM{loweredge+1};

			if ~silent
				fprintf( 'Zooming in on L%c reg path (%0.1f-%0.1f):\n', regchar, mag, mag*10 )
			end

			for nn = 2:(Nreg-1)
				regfit = xnim.set_2Dreg_params( 'subs', tar, 'filts', filtN, lambdaID, L2s(nn) );
			
				if strcmp( lambdaID, 'nld2' )
					regfit = regfit.fit_NL2d( Robs, Xs, modvarargin{:} );
				else
					regfit = regfit.fit_2Dfilters( Robs, Xs, modvarargin{:} );
				end
				fitsave{nn} = regfit;

				[LL,~,~,LLdata] = regfit.eval_model( Robs, Xs, XVindx, fit_options{:} );
				LLregs(nn) = LL-LLdata.nullLL;
				if ~silent
					fprintf( '  %8.2f: %f\n', L2s(nn), LLregs(nn) )
				end
				if nn > 2
				 if (LLregs(nn) < LLregs(nn-1)) && (LLregs(nn) < LLregs(nn-2))
				   nn = length(L2s)+1;
				 end
				end
			end
  
		else
		
			%% Use L2 list estalished in function call
			if ~silent
				fprintf( 'L%c reg path (%d): target %d, filt %d', regchar, Nreg, tar, filtN )
			end
			LLregs = zeros(Nreg,1);

			for nn = 1:length(L2s)
				regfit = xnim.set_2Dreg_params( 'sub_inds', tar, 'filts', filtN, lambdaID, L2s(nn) );
			
				if strcmp( lambdaID, 'nld2' )
					regfit = regfit.fit_NL2d( Robs, Xs, modvarargin{:} );
				else
					regfit = regfit.fit_2Dfilters( Robs, Xs, modvarargin{:} );
				end
				fitsave{nn} = regfit;
				[LL,~,~,LLdata] = regfit.eval_model( Robs, Xs, XVindx );
				LLregs(nn) = LL-LLdata.nullLL;
				if ~silent
					fprintf( '  %8.2f: %f\n', L2s(nn), LLregs(nn) )
				end
			end
		end

		[~,bestnn] = max(LLregs);
		%L2best = L2s(bestnn);
	
		xnim_out = fitsave{bestnn};
		if silent ~= 1
			fprintf( '    Subunit %d, filt %d: Best reg = %0.2f\n', tar, filtN, L2s(bestnn) )
		end
		L2s = [];
	end
end

end % function
