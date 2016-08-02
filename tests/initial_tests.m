%% create data

func_type = 'custom';   % 'add' | 'mult' | 'custom'

% signals
stim = 2*(rand(5000,2)-0.5); % random 2D signal between [-1 1]
x = stim(:,1);
y = stim(:,2);

% function and gradient evals
if strcmp(func_type, 'add')
    output = x + y;
    dfdx = ones(size(x));
    dfdy = ones(size(y));
elseif strcmp(func_type, 'mult')
    output = x .* y;
    dfdx = y;
    dfdy = x;
elseif strcmp(func_type, 'custom')
    output = x.^2 + 0.5*y.^3 + 5*x.*y;
    dfdx = 2*x + 5*y;
    dfdy = 1.5*y.^2 + 5*x;
end


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

% fit model
