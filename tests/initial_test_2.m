%% create data

num_2d_subs = 3;

func_type = {'poly', 'trig', 'mult'};
num_obs = 25000;  
% for visualization
num_ticks = {21, 19, 17}; 

for i = 1:num_2d_subs
    [stim{i}, output{i}, ~, ~, vis_info{i}] = ...
                  create_nonlinearity(func_type{i}, num_obs, num_ticks{i});
end

%% view 2d-nonlinearities
figure;

for i = 1:num_2d_subs
    subplot(num_2d_subs, 3, (i-1)*3 + 1)
    myimagesc(vis_info{i}.nl);
    title('Function')
    subplot(num_2d_subs, 3, (i-1)*3 + 2)
    myimagesc(vis_info{i}.dnl(:,:,1));
    title('x derivative')
    subplot(num_2d_subs, 3, (i-1)*3 + 3)
    myimagesc(vis_info{i}.dnl(:,:,2));
    title('y derivative')
end

%% fit XNIM

% create NIM with 4 subunits
for i = 1:num_2d_subs
    stim_params(2*i-1) = NIM.create_stim_params([1, 1, 1]);
    stim_params(2*i) = NIM.create_stim_params([1, 1, 1]);
    Xstims{2*i-1} = NIM.create_time_embedding(stim{i}(:,1), stim_params(2*i-1));
    Xstims{2*i} = NIM.create_time_embedding(stim{i}(:,2), stim_params(2*i));
end

nim = NIM(stim_params, repmat({'lin'}, 1, 2*num_2d_subs), ones(1, 2*num_2d_subs), ...
          'init_filts', repmat({1}, 1, 2*num_2d_subs), ...
          'Xtargets', 1:(2*num_2d_subs), ...
          'spkNL', 'lin', ...
          'noise_dist', 'gaussian');

% convert to an XNIM
xnim = XNIM(nim);
for i = 1:num_2d_subs
    xnim = xnim.convert_to_2d_subunit(1, 2, Xstims, 'Nticks', num_ticks{i});
end

% change default behavior of convert_to_2d_subunit for uniform signal
total_output = zeros(size(output{1}));
for i = 1:num_2d_subs
    
    a = min(Xstims{2*i-1});
    b = max(Xstims{2*i});
    xnim.twoD_subunits(i).ticks{1} = a:((b-a)/(num_ticks{i}-1)):b;
    
    a = min(Xstims{2*i});
    b = max(Xstims{2*i});
    xnim.twoD_subunits(i).ticks{2} = a:((b-a)/(num_ticks{i}-1)):b;
    
    xnim.twoD_subunits(i).ks{1} = 1;
    xnim.twoD_subunits(i).ks{2} = 1;
    xnim.twoD_subunits(i).reg_lambdas.l2 = 1;
    
    total_output = total_output + output{i};
end

% fit model
xnim = xnim.fit_NL2d(total_output, Xstims);

%% plot results

figure;
for i = 1:num_2d_subs
    dxdy = xnim.twoD_subunits(i).apply_NL_deriv([vis_info{i}.X(:), vis_info{i}.Y(:)]);

    subplot(num_2d_subs, 3, (i-1)*3 + 1)
    myimagesc(xnim.twoD_subunits(i).NL2d);
    subplot(num_2d_subs, 3, (i-1)*3 + 2)
    myimagesc(reshape(dxdy(:,1), num_ticks{i}, num_ticks{i}));
    subplot(num_2d_subs, 3, (i-1)*3 + 3)
    myimagesc(reshape(dxdy(:,2), num_ticks{i}, num_ticks{i}));

end
