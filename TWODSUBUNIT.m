classdef TWODSUBUNIT

properties
    NL2d;           % 2D nonlinearity coefficients
    k_x;            % filter coefficients along x-dir
    k_y;            % filter coefficients along y-dir
    Xtarg_x;        % index of stimulus the filter along x-dir acts on
    Xtarg_y;        % index of stimulus the filter along x-dir acts on
    xticks;         % tick values along x-dir (range of k_x output)
    yticks;         % tick values along y-dir (range of k_y output)
    reg_lambdas;    % struct of regularization hyperparameters
    Ksign_con;      % 2x1 cell array defining constraint matrices A and B
                    % such that A*TWODSUBUNIT.NL2d(:) < B
    pre_scale;      % 2x2 scale matrix
end

%% *************************** constructor ********************************
methods

    function subunit = TWODSUBUNIT(ks, stims, varargin)
    % Usage: sub2d = TWODSUBUNIT(ks, stims, varargin)
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

    subunit.Xtarg = parsed_options.Xtarg;
    subunit.ks = ks;
    if length(ks) < 2
        subunit.ks{2} = subunit.ks{1};
        subunit.ks{2}(3:end) = subunit.ks{1}(1:end-2);
        subunit.ks{2}(1:2) = 0;
    end

    [subunit,mean_gs] = subunit.normalize_filters( stims );

    % Establish scale of filters (and NL-resolution)
    % will go to plus/minus 3-stds
    dx = 6.0/parsed_options.Nbins;
    subunit.xticks = -3:dx:3;
    subunit.yticks = -3:dx:3;

    subunit.NL2d = zeros( parsed_options.Nbins, parsed_options.Nbins );
    subunit.TB = TentBasis2D( subunit.xticks, subunit.yticks );
    
    end % method
    
end
%% *************************** output methods *****************************
methods
    
    function gs = apply_filters( sub2d, stims )
    % Usage: gs = sub2d.apply_filters( stims )

    NT = size(stims{1},1);
    gs = zeros(NT,2);
    gs(:,1) = stims{sub2d.Xtarg}*sub2d.ks{1};
    gs(:,2) = stims{sub2d.Xtarg}*sub2d.ks{2};
    
    end % method

    function sub_out = apply_NL( sub2d, gs )
    % Usage: sub_out = sub2d.apply_NL( gen_signals )

    Xmat = sub2d.TB.InputNL2D( gs ); % produce X-matrix
    %sub_out = Xmat*sub2d.TB.tentcoeff(:);
    sub_out = Xmat*sub2d.NL2d(:);  % same thing but not private
    
    end

    function sub_out = apply_NL_deriv(subnit, gs)
    % Usage: sub_out = sub2d.apply_NL( gen_signals )
    
    end
    
    function sub_out = process_stim( sub2d, stims )
    % Usage: sub_out = sub2d.process_stim( stims )

    sub_out = sub2d.apply_NL( sub2d.apply_filter( stims ));
    
    end % method

end
%% *************************** display methods ****************************
methods
    
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

    end % method

    function [] = display_NL( sub2d, varargin )
    % Usage: [] = display_NL( sub2d, varargin )
    imagesc(sub2d.NL2d)
    
    end % method

end
%% *************************** hidden methods *****************************
methods (Hidden)

    function [sub2d_out,mean_gs] = normalize_filters( sub2d, stims )
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
    
    end % method

end
%% *************************** static methods *****************************
methods (Static)
    
	function reg_lambdas = init_reg_lambdas()
    % Usage: reg_lambdas = init_reg_lambdas()
	% Creates reg_lambdas struct and sets default values to 0
     
    reg_lambdas.nld2 = 0; %second derivative of tent basis coefs
    reg_lambdas.d2xt = 0; %spatiotemporal laplacian
    reg_lambdas.d2x = 0; %2nd spatial deriv
    reg_lambdas.d2t = 0; %2nd temporal deriv
    reg_lambdas.l2 = 0; %L2 on filter coefs
    reg_lambdas.l1 = 0; %L1 on filter coefs
    
    end % methods
    
end

end