function [stim, output, dfdx, dfdy, vis_info] = create_nonlinearity(func_type, num_obs, num_ticks)
% creates random 2d signals and combines them using a specified
% nonlinearity.
%
% INPUTS
%   func_type   'add' | 'mult' | 'poly' | 'trig'
%   num_obs     number of observations
%   num_ticks   number of ticks for discretized version of nonlinearity
%
% OUTPUTS
%   output      Tx1 vector of nonlinearity output
%   dfdx        Tx1 vector of partial derivatives wrt x for each obs
%   dfdy        Tx1 vector of partial derivatives wrt y for each obs
%   nl          num_ticks x num_ticks matrix of nonlinearity values
%   dnl         num_ticks x num_ticks x 2 tensor of derivatives wrt x
%               (first entry in 3rd dim) and derivatives wrt y (second
%               entry in 3rd dim)

% create random 2D signal
x = rand(num_obs,1)-0.5;
y = rand(num_obs,1)-0.5;
stim = [x, y];

% create meshgrid for nonlinearity val/grad visualization
a = min(x);
b = max(x);
xticks = a:((b-a)/(num_ticks-1)):b;
yticks = a:((b-a)/(num_ticks-1)):b;
[X, Y] = meshgrid(xticks, yticks);
X = X'; % match orientation of NL2d
Y = Y'; % match orientation of NL2d
dnl = zeros(length(xticks), length(yticks), 2);

% function and gradient evals
if strcmp(func_type, 'add')
    % signal
    output = x + y;
    dfdx = ones(size(x));
    dfdy = ones(size(y));
    % visualization
    nl = X + Y;
    dnl(:,:,1) = ones(size(X));
    dnl(:,:,2) = ones(size(Y));
elseif strcmp(func_type, 'mult')
    % signal
    output = x .* y;
    dfdx = y;
    dfdy = x;
    % visualization
    nl = X .* Y;
    dnl(:,:,1) = Y;
    dnl(:,:,2) = X;
elseif strcmp(func_type, 'poly')
    % signal
    output = x.^2 + 0.5*y.^3 + 5*x.*y;
    dfdx = 2*x + 5*y;
    dfdy = 1.5*y.^2 + 5*x;
    % visualization
    nl = X.^2 + 0.5*Y.^3 + 5*X.*Y;
    dnl(:,:,1) = 2*X + 5*Y;
    dnl(:,:,2) = 1.5*Y.^2 + 5*X;
elseif strcmp(func_type, 'trig')
    % signal
    output = cos(10*x).*sin(5*y + 2) + cos(2*x.*y + 1);
    dfdx = -10*sin(10*x).*sin(5*y + 2) - 2*y.*sin(2*x.*y + 1);
    dfdy = 5*cos(10*x).*cos(5*y + 2) - 2*x.*sin(2*x.*y + 1);
    % visualization
    nl = cos(10*X).*sin(5*Y + 2) + cos(2*X.*Y + 1);
    dnl(:,:,1) = -10*sin(10*X).*sin(5*Y + 2) - 2*Y.*sin(2*X.*Y + 1);
    dnl(:,:,2) = 5*cos(10*X).*cos(5*Y + 2) -2*X.*sin(2*X.*Y + 1);
end

vis_info.X = X;
vis_info.Y = Y;
vis_info.nl = nl;
vis_info.dnl = dnl;