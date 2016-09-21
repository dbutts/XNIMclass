%% create data

func_type = 'custom';   % 'add' | 'mult' | 'custom'

% create signals
num_obs = 500;
stim = (rand(num_obs,2)-0.5); % random 2D signal
x = stim(:,1);
y = stim(:,2);

% for discretization
num_ticks = 21;
a = min(stim(:,1));
b = max(stim(:,1));
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
elseif strcmp(func_type, 'custom')
    % signal
    output = x.^2 + 0.5*y.^3 + 5*x.*y;
    dfdx = 2*x + 5*y;
    dfdy = 1.5*y.^2 + 5*x;
    % visualization
    nl = X.^2 + 0.5*Y.^3 + 5*X.*Y;
    dnl(:,:,1) = 2*X + 5*Y;
    dnl(:,:,2) = 1.5*Y.^2 + 5*X;
end

%% view 2d-nonlinearity, normalized like XNIM
figure;
subplot(131)
myimagesc(nl);
title('true func')
subplot(132)
myimagesc(dnl(:,:,1));
title('true x-deriv')
subplot(133)
myimagesc(dnl(:,:,2));
title('true y-deriv')

%% fit XNIM

% create NIM with 2 subunits
stim_params = NIM.create_stim_params([1, 1, 1]);
stim_params(2) = NIM.create_stim_params([1, 1, 1]);
Xstims{1} = NIM.create_time_embedding(stim(:,1), stim_params(1));
Xstims{2} = NIM.create_time_embedding(stim(:,2), stim_params(2));
nim = NIM(stim_params, {'lin', 'lin'}, [1, 1], ...
          'init_filts', {1, 1}, ...
          'Xtargets', [1, 2], ...
          'spkNL', 'lin', ...
          'noise_dist', 'gaussian');

% convert to an XNIM
xnim = XNIM(nim);
xnim = xnim.convert_to_2d_subunit(1, 2, Xstims);
% change default behavior of convert_to_2d_subunit for uniform signal
a = min(Xstims{1});
b = max(Xstims{1});
xnim.twoD_subunits.ticks{1} = a:((b-a)/(num_ticks-1)):b;
a = min(Xstims{2});
b = max(Xstims{2});
xnim.twoD_subunits.ticks{2} = a:((b-a)/(num_ticks-1)):b;
xnim.twoD_subunits.ks{1} = 1;
xnim.twoD_subunits.ks{2} = 1;

xnim.twoD_subunits.reg_lambdas.d2xt = 1; % looks great with only 500 obs
% xnim.twoD_subunits.reg_lambdas.l2 = 1; % looks much worse

% fit model
xnim = xnim.fit_NL2d(output, Xstims);

% plot results
dxAdy = xnim.twoD_subunits(1).apply_NL_deriv([X(:), Y(:)]);
figure;
subplot(131)
myimagesc(xnim.twoD_subunits(1).NL2d);
title('est. func')
subplot(132)
myimagesc(reshape(dxdy(:,1), num_ticks, num_ticks));
title('est. x-deriv')
subplot(133)
myimagesc(reshape(dxdy(:,2), num_ticks, num_ticks));
title('est. y-deriv')