classdef Autoencoder
    
% Object-oriented implementation of an autoencoder with a single layer 
%   (Notes)
%
% Reference: Lots of good information on basic autoencoder implementations
% (like this one) can be found in Stanford's UFLDL tutorial wiki
% http://deeplearning.stanford.edu/wiki/index.php/Neural_Networks
%
% Author: Matt whiteway
%   03/03/16
%
% TODO
%   example script

properties
    data_params;         % struct defining the data parameters
        % num_ex
        % num_io_nodes
    net_params;          % struct defining the network parameters
        % num_hid_nodes
        % act_func_hid  
        % act_func_out  
        % weight_tie
    reg_params;          % struct defining the regularization parameters
        % l2_weight
        % l2_bias1 
        % l2_bias2 
        % sp_hid_act
        % kl_param
    optim_params;        % struct defining the optimization parameters
        % opt_routine
        % batch_size
        % max_iters
        % Display
        % monitor
    fit_history         % struct with output of optimization routine
    init_params
        % init_weights
        % rng_state
        
    w1;                  % matrix of weights from input to hidden layer
    w2;                  % matrix of weights from hidden to output layer
    b1;                  % vector of bias weights from input to hidden layer
    b2;                  % vector of bias weights from hidden to output layer
end


properties (Hidden)
    version = '0.2';     % version number
    date = date;         % date of fit
    allowed_act = {'lin','relu','sigmoid'};		% allowed activation functions
    allowed_init_opts = {'gaussian','uniform'}; % allowed initializations
    allowed_opt_routines = {'sd','mat_sd','cg','lbfgs','mat_bfgs','con','sgd','opti'};
end


    
%% ********************  constructor ********************************
methods
    
    function net = Autoencoder(init_params,varargin)
    % net = Autoencoder(init_params,'key1_handle','key1_value',...)
    %
    % constructor function for an ExtendedAutoencoder object; sets
    % properties, the defaults of which can be found below. Optional flags 
    % allow modification of defaults upon initializing an instance of the
    % object.
    % INPUTS:
    %   init_params:    1x2 vector [num_ex, num_io_nodes]
    %   optional_flags:
    %
    %       ('num_io_nodes',num_io_nodes): number of input/output layer nodes
    %       ('num_hid_nodes',num_hid_nodes): number of hidden layer nodes
    %       ('act_func_hid','act_func_hid'): string specifying hidden layer
    %           activation functions; 'linear','relu' or 'sigmoid' supported
    %       ('act_func_out','act_func_out'): string specifying output layer
    %           activation functions; 'linear','relu' or 'sigmoid' supported
    %       ('weight_tie',value): 1 to force encoding and decoding weights
    %           to be the same, 0 otherwise
    %
    %       ('l2_weight',value): scalar specifying the penalty used for
    %           L2 norm of weights
    %       ('l2_bias1',value): scalar specifying the penalty used for
    %           L2 norm of biases into hidden layer
    %       ('l2_bias2',value): scalar specifying the penalty used for
    %           L2 norm of biases into output layer
    %       ('sp_hid_act',value): scalar specifying the penalty used for the
    %           activation values of the hidden layer
    %       ('kl_param',value): scalar in [0 1] specifying the KL sparsity 
    %           penalty value
    %
    %       ('opt_routine',opt_routine): string specifying the optimization
    %           routine used for learning the weights and biases.
    %           'sd','cg','lbfgs' and 'sgd' supported
    %       ('batch_size',batch_size): number of examples to use for each
    %           gradient descent step if using sgd
    %       ('Display',Display): 'off' to suppress output, 'iter' for output at each
    %           iteration, 'batch' for output after each pass through data (in sgd)
    %       ('maxIter',maxIter): maximum number of iterations of opt
    %           routine
    %       ('monitor',monitor): 'off' to suppress saving output, 'iter' to save output 
    %           at each iteration, 'batch' to save output after each pass through 
    %           data (in sgd), 'both' to save both
    % 
    %       ('init_weights',input): input is either a string specifying
    %           type of random initialization for weights ('gaussian' or 
    %           'uniform') or a weight vector of appropriate length
    
    if nargin == 0
        return % handle the no-input-argument case by returning a null model. This is important when initializing arrays of objects
    end

    % handle init_params
    assert(~isempty(init_params.num_ex),'init_params must contain number of examples');
    assert(~isempty(init_params.num_io_nodes),'init_params must contain number of predictors');
    num_ex = init_params.num_ex;
    num_io_nodes = init_params.num_io_nodes;
    
    % DEFINE DEFAULTS
    % net_params
    num_hid_nodes = 1;          % number of hidden layer nodes
    act_func_hid  = 'relu';     % type of activation function in hidden layer
    act_func_out  = 'lin';      % type of activation function in output layer
    weight_tie    = 0;          % tie encoding and decoding weights
    % reg_params
    l2_weight  = 0.01;          % L2 penalty shared by weights in first and second layers
    l2_bias1   = 0.01;          % L2 penalty on biases in first layer
    l2_bias2   = 0.01;          % L2 penalty on biases in second layer
    sp_hid_act = 0.1;           % sparsity penalty on hidden layer outputs (L1 for lin/rect_lin, KL for sigmoid)
    kl_param   = 0.05;          % target probability of activation for sigmoidal hidden layer nodes
    % optim_params defaults
    opt_routine = 'lbfgs';      % optimization routine to use
    batch_size  = 10;           % batch size for use in sgd routine
    maxIter     = 1000;         % maximum number of iterations in optimization routine
    Display     = 'off';        % control optimization routine output
    monitor     = 'iter';       % save optimization routine output
    % init_params   
    init_weights = 'uniform';   % what weights to use for optimization routine
    
    % PARSE INPUTS
    i = 1;
    while i <= length(varargin)
        switch lower(varargin{i})
            % net_params
            case 'num_hid_nodes'
                assert(varargin{i+1} > 0,...
                    'must provide positive number of hidden nodes')
                num_hid_nodes = varargin{i+1};
            case 'num_io_nodes'
                assert(varargin{i+1}==net.data_params.num_io_nodes,...
                    'number of input/output layer nodes must match input dimension of data')
                num_io_nodes = varargin{i+1};
            case 'act_func_hid'
                assert(ismember(varargin{i+1},net.allowed_act),...
                    'invalid activation function');
                act_func_hid = varargin{i+1};
            case 'act_func_out'
                assert(ismember(varargin{i+1},net.allowed_act),...
                    'invalid activation function');
                act_func_out = varargin{i+1};
            case 'weight_tie'
                assert(ismember(varargin{i+1},[0,1]),...
                    '''weight_tie'' input must be 0 or 1');
                weight_tie = varargin{i+1};
            % reg_params
            case 'l2_weight'
                assert(varargin{i+1} > 0,...
                    'reg value must be positive')
                l2_weight = varargin{i+1};
            case 'l2_bias1'
                assert(varargin{i+1} > 0,...
                    'reg value must be positive')
                l2_bias1 = varargin{i+1};
            case 'l2_bias2'
                assert(varargin{i+1} > 0,...
                    'reg value must be positive')
                l2_bias2 = varargin{i+1};
            case 'sp_hid_act'
                assert(varargin{i+1} > 0,...
                    'reg value must be positive')
                sp_hid_act = varargin{i+1};
            case 'kl_param'
                assert(varargin{i+1} >= 0 && varargin{i+1} <= 1,...
                    'KL param must lie in interval [0,1]')
                kl_param = varargin{i+1};
            % optim_params
            case 'opt_routine'
                assert(ismember(varargin{i+1},net.allowed_opt_routines),...
                    'invalid optimization routine')
                opt_routine = varargin{i+1};
            case 'batch_size'
                assert(mod(net.data_params.num_ex,varargin{i+1}) == 0,...
                    'batch size should evenly divide number of examples')
                batch_size = varargin{i+1};
            case 'display'
                assert(ismember(varargin{i+1},{'off','iter','batch'}),...
                    'invalid parameter for ''Display''')
                Display = varargin{i+1};
            case 'maxiter'
                assert(varargin{i+1} > 0,...
                    'maximum number of iterations must be greater than zero')
                maxIter = varargin{i+1};
            case 'monitor'
                assert(ismember(varargin{i+1},{'off','iter','batch','both'}),...
                    'invalid parameter for ''monitor''')
                monitor = varargin{i+1};
            % init_params
            case 'init_weights'
                if ischar(varargin{i+1})
                    assert(ismember(varargin{i+1},{'gaussian','uniform'}),'invalid random init option')
                    init_weights = varargin{i+1};
                elseif isvector(varargin{i+1})
                    % do size checking later
                    init_weights = varargin{i+1}; 
                else
                    error('''init_weights'' must be a string or a vector')
                end
            otherwise
                error('Invalid input flag');
        end
        i = i+2;
    end
    
    % INITIALIZE wEIGHTS
    [w1_,w2_,b1_,b2_,init_params] = Autoencoder.set_init_weights_stat(...
                        init_weights,num_io_nodes,num_hid_nodes,weight_tie);
    
    % SET PROPERTIES
    % data_params
    net.data_params.num_ex       = num_ex;
    net.data_params.num_io_nodes = num_io_nodes;
    % net_params
    net.net_params.num_hid_nodes = num_hid_nodes; 
    net.net_params.num_io_nodes  = num_io_nodes;
    net.net_params.act_func_hid  = act_func_hid;    
    net.net_params.act_func_out  = act_func_out;     
    net.net_params.weight_tie    = weight_tie;        
    % reg_params
    net.reg_params.l2_weight     = l2_weight;        
    net.reg_params.l2_bias1      = l2_bias1;        
    net.reg_params.l2_bias2      = l2_bias2;           
    net.reg_params.sp_hid_act    = sp_hid_act;           
    net.reg_params.kl_param      = kl_param;            
    % optim_params defaults
	net.optim_params 			 = Autoencoder.set_init_optim_params();
    net.optim_params.opt_routine = opt_routine;    
    net.optim_params.batch_size  = batch_size;           
    net.optim_params.maxIter     = maxIter;
    net.optim_params.maxFunEvals = 2*maxIter;         
	net.optim_params.Display 	 = Display;
    net.optim_params.monitor     = monitor;       
    % fit_history
    net.fit_history = struct([]);
	% init_params   
    net.init_params              = init_params;           
    % weights
    net.w1 = w1_;
    net.w2 = w2_;
    net.b1 = b1_;
    net.b2 = b2_;
    
    end

end
%% ********************  setting methods ********************************
methods  

    function net = set_net_params(net,varargin)
    % net = net.set_net_params('key1',key1_value,'key2',key2_value,...)
    %
    % Takes a sequence of key-value pairs to set network parameters for an 
    % Autoencoder object
    % INPUTS:
    %   ('num_io_nodes',num_io_nodes): number of input/output layer nodes
    %   ('num_hid_nodes',num_hid_nodes): number of hidden layer nodes
    %   ('act_func_hid','act_func_hid'): string specifying hidden layer
    %       activation functions; 'linear','relu' or 'sigmoid' supported
    %   ('act_func_out','act_func_out'): string specifying output layer
    %       activation functions; 'linear','relu' or 'sigmoid' supported
    %   ('weight_tie',value): 1 to force encoding and decoding weights
    %       to be the same, 0 otherwise
    %   
    % OUTPUTS:
    %   net: updated Autoencoder object
    
    % check for appropriate number of inputs
    assert(mod(length(varargin),2)==0,'Input should be a list of key-value pairs')
    
    i = 1;
    while i <= length(varargin)
        switch lower(varargin{i})
            case 'num_hid_nodes'
                assert(varargin{i+1} > 0,...
                    'must provide positive number of hidden nodes')
                if varargin{i+1} ~= net.net_params.num_hid_nodes
                    %warning('changing number of hidden nodes; randomly reinitializing weights!')
                    net.net_params.num_hid_nodes = varargin{i+1};
                    net = net.set_init_weights('uniform');
                end
            case 'num_io_nodes'
                assert(varargin{i+1}==net.data_params.num_io_nodes,...
                    'number of input/output layer nodes must match input dimension of data')
                net.net_params.num_io_nodes = varargin{i+1};
            case 'act_func_hid'
                assert(ismember(varargin{i+1},net.allowed_act),...
                    'invalid activation function');
                net.net_params.act_func_hid = varargin{i+1};
            case 'act_func_out'
                assert(ismember(varargin{i+1},net.allowed_act),...
                    'invalid activation function');
                net.net_params.act_func_out = varargin{i+1};
            case 'weight_tie'
                assert(ismember(varargin{i+1},[0,1]),...
                    '''weight_tie'' input must be 0 or 1');
                if varargin{i+1} < net.net_params.weight_tie
                    % introducing new weights
                    net.net_params.weight_tie = varargin{i+1};
                    fprintf('Removing weight-tying constraint! Enter:\n')
                    fprintf('\t0 to reinitialize ALL weights randomly\n')
                    fprintf('\t1 to copy current weights\n')
                    resp = input('(0/1): ');
                    if resp == 1
                        net.w1 = net.w2';
                    elseif resp == 0
                        net = net.set_init_weights('uniform');
                    else
                        warning('invalid input; reinitializing all weights randomly')
                        net = net.set_init_weights('uniform');
                    end
                elseif varargin{i+1} > net.net_params.weight_tie
                    % introducing weight-tying
                    net.net_params.weight_tie = varargin{i+1};
                    fprintf('Introducing weight-tying constraint! Enter:\n')
                    fprintf('\t0 to reinitialize ALL weights randomly\n')
                    fprintf('\t1 to use current weights\n')
                    resp = input('(0/1): ');
                    if resp == 1
                        % no changes
                    elseif resp == 0
                        net = net.set_init_weights('uniform');
                    else
                        warning('invalid input; reinitializing all weights randomly')
                        net = net.set_init_weights('uniform');
                    end
                else
                    warning('No change in weight_tie; ignoring input')
                end
            otherwise
                error('Invalid input flag')
        end
        i = i + 2;
    end
    end
    
    function net = set_reg_params(net,varargin)
    % net = net.set_reg_params('key1',key1_value,...)
    %
    % Takes a sequence of key-value pairs to set regularization parameters 
    % for an Autoencoder object
    % INPUTS:
    %   ('l2_weight',value): scalar specifying the penalty used for
    %       L2 norm of weights
    %   ('l2_bias1',value): scalar specifying the penalty used for
    %       L2 norm of biases into hidden layer
    %   ('l2_bias2',value): scalar specifying the penalty used for
    %       L2 norm of biases into output layer
    %   ('sp_hid_act',value): scalar specifying the penalty used for the
    %       activation values of the hidden layer
    %   ('kl_param',value): scalar in [0 1] specifying the KL sparsity 
    %       penalty value
    %
    % OUTPUTS:
    %   net: updated Autoencoder object
    
    % check for appropriate number of inputs
    assert(mod(length(varargin),2)==0,'Input should be a list of key-value pairs')
    
    i = 1;
    while i <= length(varargin)
        switch lower(varargin{i})
            case 'l2_weight'
                assert(varargin{i+1} >= 0,...
                    'reg value must be nonnegative')
                net.reg_params.l2_weight = varargin{i+1};
            case 'l2_bias1'
                assert(varargin{i+1} >= 0,...
                    'reg value must be nonnegative')
                net.reg_params.l2_bias1 = varargin{i+1};
            case 'l2_bias2'
                assert(varargin{i+1} >= 0,...
                    'reg value must be nonnegative')
                net.reg_params.l2_bias2 = varargin{i+1};
            case 'sp_hid_act'
                assert(varargin{i+1} >= 0,...
                    'reg value must be nonnegative')
                net.reg_params.sp_hid_act = varargin{i+1};
            case 'kl_param'
                assert(varargin{i+1} >= 0 && varargin{i+1} <= 1,...
                    'KL param must lie in interval [0,1]')
                net.reg_params.kl_param = varargin{i+1};
            otherwise
                error('Invalid input flag');
        end
        i = i + 2;
    end
    end
        
    function net = set_optim_params(net,varargin)
    % net = net.set_optim_params('key1',key1_value,...)
    %
    % Takes a sequence of key-value pairs to set optimization parameters 
    % for an Autoencoder object
    % INPUTS:
    %   ('opt_routine',opt_routine): string specifying the optimization
    %       routine used for learning the weights and biases.
    %       'sd','cg','lbfgs','sgd' and 'con' supported
    %   ('batch_size',batch_size): number of examples to use for each
    %       gradient descent step if using sgd
    %   ('Display',Display): 'off' to suppress output, 'iter' for output at each
    %       iteration, 'batch' for output after each pass through data (in sgd)
    %   ('maxIter',maxIter): maximum number of iterations of opt
    %       routine
    %   ('monitor',monitor): 'off' to suppress saving output, 'iter' to save output 
    %       at each iteration, 'batch' to save output after each pass through 
    %       data (in sgd), 'both' to save both
    %
    % OUTPUTS:
    %   net: updated Autoencoder object
    
    % check for appropriate number of inputs
    assert(mod(length(varargin),2)==0,'Input should be a list of key-value pairs')
    
    % allowed optimization routines
    i = 1;
    while i <= length(varargin)
        switch lower(varargin{i})
            case 'opt_routine'
                assert(ismember(varargin{i+1},net.allowed_opt_routines),...
                    'invalid optimization routine')
                net.optim_params.opt_routine = varargin{i+1};
				% set necessary fields of optim_params for matlab routines
				if strcmp(varargin{i+1},'mat_sd')
					net.optim_params.Algorithm = 'quasi-newton';
					net.optim_params.HessUpdate = 'steepdesc';
					net.optim_params.GradObj = 'on';
				elseif strcmp(varargin{i+1},'mat_bfgs')
					net.optim_params.Algorithm = 'quasi-newton';
					net.optim_params.HessUpdate = 'bfgs';
					net.optim_params.GradObj = 'on';
				end
            case 'batch_size'
                assert(mod(net.data_params.num_ex,varargin{i+1}) == 0,...
                    'batch size should evenly divide number of examples')
                net.optim_params.batch_size = varargin{i+1};
            case 'display'
                assert(ismember(varargin{i+1},{'off','iter','batch'}),...
                    'invalid parameter for ''Display''')
                net.optim_params.Display = varargin{i+1};
            case 'maxiter'
                assert(varargin{i+1} > 0,...
                    'maximum number of iterations must be greater than zero')
                net.optim_params.maxIter = varargin{i+1};
                net.optim_params.maxFunEvals = 2*varargin{i+1};
            case 'monitor'
                assert(ismember(varargin{i+1},{'off','iter','batch','both'}),...
                    'invalid parameter for ''monitor''')
                net.optim_params.monitor = varargin{i+1};
            otherwise
                error('Invalid input flag');
        end
        i = i + 2;
    end
    end
    
    function net = set_hidden_order(net,data)
    % net = net.set_hidden_order(data)
    %
    % reorders the hidden nodes to represent highest to lowest stddev of 
    % output
    % INPUTS:
    %   data:           num_ex x num_io_nodes matrix of data
    %   
    % OUTPUTS:
    %   net:            updated Autoencoder object

    % ensure weights exist
    assert(~isempty(net.w2),'Autoencoder object must have weights to reorder')
    
    % get stddev of hidden units
    fgint = net.get_model_internals(data);
    activity = fgint{1};
    activity = std(activity,[],1);
    % sort by stddev (and flip to go from highest to lowest)
    [~,sorted_index] = sort(activity);
    sorted_index = fliplr(sorted_index);
    % rearrange weights and biases
    net.b1 = net.b1(sorted_index);
    net.w1 = net.w1(:,sorted_index); % only reorder w1 if it exists
    net.w2 = net.w2(sorted_index,:);
    
    end
    
    function net = set_init_weights(net,init_weights)
    % net = net.set_init_weights(init_weights)
    %
    % function that sets w1, w2, b1 and b2 properties of Autoencoder
    % object
    % INPUTS:
    %   init_weights:   either a string specifying type of random 
    %                   initialization for weights ('gaussian' or 
    %                   'uniform') or a weight vector of appropriate length  
    %
    % OUTPUT:
    %   net:            updated Autoencoder object
    
    % call static set_init_weights used in Constructor
    [w1_,w2_,b1_,b2_,init_params_struct] = ...
    Autoencoder.set_init_weights_stat(init_weights,net.net_params.num_io_nodes,...
                    net.net_params.num_hid_nodes,net.net_params.weight_tie);
    % set properties
    net.w1 = w1_;
    net.w2 = w2_;
    net.b1 = b1_;
    net.b2 = b2_;
    net.init_params = init_params_struct;
    end
    
end
%% ********************  getting methods ********************************
methods
      
    function [fgint, gint] = get_model_internals(net,data,varargin)
    % net = net.get_model_internals(data,'key1','key1_value',...)
    %
    % function that retrieves informaton from autoencoder model
    % INPUTS:
    %   data:                 num_ex x num_io_nodes matrix of data
    %   optional_flags:
    %       ('indx_tr',indx): subset of 1:num_io_nodes that specifies 
    %           portion of data used for evaluation, or set to NaN to use
    %           all data
    %
    % OUTPUTS:
    %   fgint       2x1 cell array, each cell containing a 
    %               num_ex x num_io_nodes matrix of the signal after being
    %               passed through the activation function
    %   gint        2x1 cell array, each cell containing a 
    %               num_ex x num_io_nodes matrix of the signal before being
    %               passed through the activation function

    % check inputs
    net.check_inputs(data);
    
    % define defaults
    indx_tr = NaN; % NaN means we use all available data
    
    % parse inputs
    i = 1;
    while i <= length(varargin)
        switch lower(varargin{i})
            case 'indx_tr'
                assert(all(ismember(varargin{i+1},1:size(data,1))) ...
                    || isnan(varargin{i+1}),'invalid fitting indices')
                indx_tr = varargin{i+1};
            otherwise
                error('invalid input flag');
        end
        i = i+2;
    end
    
    % use indx_tr
    if ~isnan(indx_tr)
        data = data(indx_tr,:);
    end
    
    % set aside constants to keep things clean
    w1_ = net.w1;
    w2_ = net.w2;
    b1_ = net.b1;
    b2_ = net.b2;
    
    % evaluate model
	gint{1}  = bsxfun(@plus,data*w1_,b1_');
	fgint{1} = net.apply_act_func(gint{1},1);
	gint{2}  = bsxfun(@plus,fgint{1}*w2_,b2_');
	fgint{2} = net.apply_act_func(gint{2},2);
    
	end
 
    function reg_pen = get_reg_pen(net,data,varargin)
    % net = net.get_reg_pen(data,'key1','key1_value',...)
    %
    % function that retrieves informaton from autoencoder model
    % INPUTS:
    %   data:                 num_ex x num_io_nodes matrix of data
    %   optional_flags:
    %       ('indx_tr',indx): subset of 1:num_io_nodes that specifies 
    %           portion of data used for evaluation, or set to NaN to use
    %           all data
    %
    % OUTPUTS:
    %   reg_pen:        struct with fields for penalties associated with
    %                   each regularization type

    % check inputs
    net.check_inputs(data);
    
    % define defaults
    indx_tr = NaN; % NaN means we use all available data
    
    % parse inputs
    i = 1;
    while i <= length(varargin)
        switch lower(varargin{i})
            case 'indx_tr'
                assert(all(ismember(varargin{i+1},1:size(data,1))) ...
                    || isnan(varargin{i+1}),'invalid fitting indices')
                indx_tr = varargin{i+1};
            otherwise
                error('invalid input flag');
        end
        i = i+2;
    end
    
    % use indx_tr
    if ~isnan(indx_tr)
        data = data(indx_tr,:);
    end

	% set aside constants to keep things clean
    w1_ = net.w1;
    w2_ = net.w2;
    b1_ = net.b1;
    b2_ = net.b2;
    lambda_w  = net.reg_params.l2_weight;
    lambda_b1 = net.reg_params.l2_bias1;
    lambda_b2 = net.reg_params.l2_bias2;
    beta      = net.reg_params.sp_hid_act;
    
    % sparsity penalty eval on hidden layer
    fgint = net.get_model_internals(data);
    sparse_func = net.get_sparse_penalty(fgint{1});
    
    % get penalty terms
	Z = numel(data);
	reg_pen.l2_weight1 = 0.5*lambda_w*sum(sum(w1_.^2))/Z; % L2 penalty on layer 1 weights
    reg_pen.l2_weight2 = 0.5*lambda_w*sum(sum(w2_.^2))/Z; % L2 penalty on layer 2 weights
    reg_pen.l2_bias1 = 0.5*lambda_b1*sum(b1_.^2)/Z;       % L2 penalty on layer 1 biases
    reg_pen.l2_bias2 = 0.5*lambda_b2*sum(b2_.^2)/Z;       % L2 penalty on layer 2 biases
    reg_pen.sp_act_hid = beta*sparse_func/Z;              % sparsity penalty   
    end
    
    function [r2s,cost_func,mod_internals,reg_pen] = get_model_eval(net,data,varargin)
    % net = net.get_model_internals(data,'key1','key1_value',...)
    %
    % function that retrieves informaton from autoencoder model
    % INPUTS:
    %   data:                 num_ex x num_io_nodes matrix of data
    %   optional_flags:
    %       ('indx_tr',indx): subset of 1:num_io_nodes that specifies 
    %           portion of data used for evaluation, or set to NaN to use
    %           all data
    %
    % OUTPUTS:
    %   r2s:            r squared values for all output nodes
    %   cost_func:      value of cost function given current weights and input
    %                   data
    %   mod_internals:  struct with fields gint and fgint holding output of
    %                   get_model_internals method
    %   reg_pen:        struct with fields for penalties associated with
    %                   each regularization type

    % check inputs
    net.check_inputs(data);
    
    % define defaults
    indx_tr = NaN; % NaN means we use all available data
    
    % parse inputs
    i = 1;
    while i <= length(varargin)
        switch lower(varargin{i})
            case 'indx_tr'
                assert(all(ismember(varargin{i+1},1:size(data,1))) ...
                    || isnan(varargin{i+1}),'invalid fitting indices')
                indx_tr = varargin{i+1};
            otherwise
                error('invalid input flag');
        end
        i = i+2;
    end
    
    % use indx_tr
    if ~isnan(indx_tr)
        data = data(indx_tr,:);
    end
    
    [fgint, gint] = net.get_model_internals(data);
    
    % evaluate cost_func and r2
    cost_func = 0.5*sum(sum((data-fgint{2}).^2))/numel(data);
    r2s = (1-sum((data-fgint{2}).^2,1)./...
         sum((bsxfun(@minus,data,mean(data,1))).^2,1))';

    % put internal generating functions in a struct if desired in output
    if nargout > 2
        mod_internals = struct('gint',gint,'fgint',fgint);
    end
    % put reg penalties in a struct if desired in output
    if nargout > 3
        reg_pen = net.get_reg_pen(data);
    end
    end
   
    function [w1, w2, b1, b2] = get_weights(net,weight_vec)
    % net = net.get_weights(weight_vec)
    %
    % function that takes a column vector of the network parameters
    % (weights and biases) and returns the separated parameters in matrix
    % and vector form, respectively
    % INPUTS:
    %     weight_vec:   column vector of network parameters in the form
    %                   [-----w1-----|------w2------|--b1--|--b2--]'
    %
    % OUTPUTS:
    %     w1:           num_io_nodes x num_hid_nodes matrix of weights
    %     w2:           num_hid_nodes x num_io_nodes matrix of weights
    %     b1:           num_hid_nodes x 1 vector of hidden layer biases
    %     b2:           out_layer_nodes x 1 vector of output layer biases

	% should do some error checking but this is mostly used as an internal 
	% function during optimization
    if net.net_params.weight_tie % no w1 matrix; just [---w1----|-b1-|-b2-]
        % get cut points
        cut1 = net.net_params.num_hid_nodes*net.net_params.num_io_nodes;
        cut2 = cut1+net.net_params.num_hid_nodes;
        % get weights
        w2 = weight_vec(1:cut1);
        w2 = reshape(w2,net.net_params.num_hid_nodes,net.net_params.num_io_nodes);
        w1 = w2';
        % get biases
        b1 = weight_vec(cut1+1:cut2);
        b2 = weight_vec((cut2+1):end);
    else
        % get cut points
        cut1 = net.net_params.num_hid_nodes*net.net_params.num_io_nodes;
        cut2 = 2*cut1;
        cut3 = cut2+net.net_params.num_hid_nodes;
        % get weights
        w1 = weight_vec(1:cut1);
        w1 = reshape(w1,net.net_params.num_io_nodes,net.net_params.num_hid_nodes);
        w2 = weight_vec(cut1+1:cut2);
        w2 = reshape(w2,net.net_params.num_hid_nodes,net.net_params.num_io_nodes);
        % get biases
        b1 = weight_vec((cut2+1):cut3);
        b2 = weight_vec((cut3+1):end);
    end
    end

end
%% ********************  fitting methods ********************************
methods
    
    function net = fit_weights(net,data,varargin)
    % net = net.fit_weights(data)
	%    
	% INPUTS:
    %   data:   num_ex x num_io_nodes matrix of data
	%   optional_flags:
    %       ('indx_tr',indx): subset of 1:num_ex that specifies portion of
    %           data used for fitting
    %       ('optim_params',params): struct of desired optimization
    %           parameters that are used to override defaults
    %
	% OUTPUTS:
    %   net:    updated Autoencoder object

	% check inputs
	net.check_inputs(data);

    % define defaults
    indx_tr = NaN;  % use all indices for fitting
    
    % parse optional inputs
    i = 1;
    while i <= length(varargin)
        switch lower(varargin{i})
            case 'indx_tr'
                assert(all(ismember(varargin{i+1},1:size(data,1))) ...
                    && ~any(isnan(varargin{i+1})),'Invalid fitting indices')
                indx_tr = varargin{i+1};
            case 'optim_params'
                net.optim_params = varargin{i+1}; % TODO
            otherwise
                error('invalid input flag');
        end
        i = i+2;
    end
   
    % use indx_tr
    if ~isnan(indx_tr)
        data = data(indx_tr,:);
    end
    % update data_params
    net.data_params.num_ex = size(data,1);
    
	% define constants to be used in nested function
    lambda_w  = net.reg_params.l2_weight;
    lambda_b1 = net.reg_params.l2_bias1;
    lambda_b2 = net.reg_params.l2_bias2;
    beta      = net.reg_params.sp_hid_act;
	Z = numel(data); % will normalize by # examples and # dims

	% define init params
	if net.net_params.weight_tie
		init_weights = [net.w2(:); net.b1; net.b2];
	else
		init_weights = [net.w1(:); net.w2(:); net.b1; net.b2];
	end
	
	% if checking derivatives
	if 0
		net.optim_params.opt_routine = 'mat_sd';
		net.optim_params.Algorithm = 'quasi-newton';
		net.optim_params.HessUpdate = 'steepdesc';
		net.optim_params.GradObj = 'on';
		net.optim_params.DerivativeCheck = 'on';
		net.optim_params.FinDiffType = 'central';
	end

	% define function handle to pass to optimizer
	obj_fun = @(x) objective_fun(x);    
	
    % net is passed in automatically, we need to specify an initial point
    % and data
    switch net.optim_params.opt_routine
        case 'mat_sd'
            [weights,f,~,output] = fminunc(obj_fun,init_weights,net.optim_params);
        case 'mat_bfgs'
            [weights,f,~,output] = fminunc(obj_fun,init_weights,net.optim_params);
        case 'sd'
            [weights,f,~,output] = minFunc(obj_fun,init_weights,net.optim_params);
        case 'cg'
            [weights,f,~,output] = minFunc(obj_fun,init_weights,net.optim_params);
        case 'lbfgs'
            [weights,f,~,output] = minFunc(obj_fun,init_weights,net.optim_params);
        case 'con'
%             [weights,f,~,output] = minConf_SPG(obj_fun,init_weights,@(t,b)max(t,0),net.optim_params);
            [weights,f,output] = minConf_PQN(obj_fun,init_weights,@(t,b)max(t,0),net.optim_params);
            %[weights,f,~,output] = minConf_BBST(obj_fun,init_weights,@(t,b) max(t,0),net.optim_params);
        case 'sgd'
            obj_fun = @(x,data) objective_fun_sgd(x,data);
            [weights,f,output] = sgd(obj_fun,init_weights,net.optim_params,data);
%         case 'opti'
%             % Options
%             opts = optiset('display','iter');
%             
%             % construct A and b for inequality
%             
%             nW1 = numel(net.w1);
%             nW2 = numel(net.w2);
%             nB1 = numel(net.b1);
%             nB2 = numel(net.b2);
%             if net.net_params.weight_tie
%                 A = diag([-ones(nW1,1) ; zeros(nB1+nB2,1)]);
%                 b = zeros(nW1+nB1+nB2,1);
%             else
%                 A = diag([-ones(nW1+nW2,1) ; zeros(nB1+nB2,1)]);
%                 b = zeros(nW1+nW2+nB1+nB2,1);
%             end
%             
%             % Build OPTI Problem
%             Opt = opti('fun',obj_fun,'ineq',A,b,'x0',init_weights,'options',opts);
%             
%             % Solve NLP
%             [weights,f,exitflag,output] = solve(Opt);
%             
    end

	% see if routine has converged
	[~,grad_pen] = objective_fun(weights);
	first_order_optim = max(abs(grad_pen));
	if first_order_optim > 1e-2
		warning('first-order optimality: %.3f, fit might not be converged!',first_order_optim);
	end

	% parse outputs of opt routine
    [net.w1,net.w2,net.b1,net.b2] = net.get_weights(weights);
    net.fit_history = cat(1,net.fit_history,struct('cost_func',f,'output',output));
	% end of fit_weights
 
		function [func, grad] = objective_fun(x)
        % function for calculating the mean square cost function and
        % the gradient with respect to the weights and biases for the 2
        % layer autoencoder network.
        % INPUTS:
        % OUTPUTS:

        [w1_, w2_, b1_, b2_] = net.get_weights(x);

        % Note: tradeoff b/t speed and memory here. Saving gint1 and
        % gint2 adds a larger memory requirement, but each of these
        % quantities is used twice; rough estimate on speedup is ~20%
        % when using gint1/2
        
        % FUNCTION
        % get activation of hidden and output layers
        gint1 = bsxfun(@plus,data*w1_,b1_');
        a2 = net.apply_act_func(gint1,1);
        gint2 = bsxfun(@plus,a2*w2_,b2_');
        a3 = net.apply_act_func(gint2,2);
        % cost function eval
        func = 0.5*sum(sum((data-a3).^2));    % cost function
          
        % GRADIENT
        % sparsity penalty eval on hidden layer
        [sparse_func, sparse_grad] = net.get_sparse_penalty(a2);
        % grad wrt b2 & b1
        gb2 = -net.apply_act_deriv(gint2,2).*(data-a3);
        gb1 = net.apply_act_deriv(gint1,1).*(bsxfun(@plus,gb2*w2_',beta*sparse_grad'));
        % grad wrt w2 & w1
        if net.net_params.weight_tie
            gw1 = [];
            gw2 = a2'*gb2+gb1'*data;
        else
            gw1 = data'*gb1;
            gw2 = a2'*gb2;
        end
        gb1 = sum(gb1,1)';
        gb2 = sum(gb2,1)';
        grad = [gw1(:); gw2(:); gb1; gb2];

        % REG PENALTIES
        reg_pen_func = 0.5*lambda_w*sum(sum(w1_.^2))+...
                       0.5*lambda_w*sum(sum(w2_.^2))+...  % L2 penalty on weights
                       0.5*lambda_b1*sum(b1_.^2)+...
                       0.5*lambda_b2*sum(b2_.^2)+...      % L2 penalty on biases
                       beta*sparse_func;                 % sparsity penalty
        if net.net_params.weight_tie
            reg_pen_grad = [2*lambda_w*w2_(:); lambda_b1*b1_; lambda_b2*b2_];
        else
            reg_pen_grad = [lambda_w*w1_(:); lambda_w*w2_(:); lambda_b1*b1_; lambda_b2*b2_];
        end

        func = (func + reg_pen_func)/Z;
        grad = (grad + reg_pen_grad)/Z;
    
        end % internal function
		
		
		function [func, grad] = objective_fun_sgd(x,data)
        % function for calculating the mean square cost function and
        % the gradient with respect to the weights and biases for the 2
        % layer autoencoder network. for use with the sgd routine, which
        % for generality requires a data input variable
        % INPUTS:
        % OUTPUTS:

        [w1_, w2_, b1_, b2_] = net.get_weights(x);

        % Note: tradeoff b/t speed and memory here. Saving gint1 and
        % gint2 adds a larger memory requirement, but each of these
        % quantities is used twice; rough estimate on speedup is ~20%
        % when using gint1/2
        
        % FUNCTION
        % get activation of hidden and output layers
        gint1 = bsxfun(@plus,data*w1_,b1_');
        a2 = net.apply_act_func(gint1,1);
        gint2 = bsxfun(@plus,a2*w2_,b2_');
        a3 = net.apply_act_func(gint2,2);
        % cost function eval
        func = 0.5*sum(sum((data-a3).^2));    % cost function
          
        % GRADIENT
        % sparsity penalty eval on hidden layer
        [sparse_func, sparse_grad] = net.get_sparse_penalty(a2);
        % grad wrt b2 & b1
        gb2 = -net.apply_act_deriv(gint2,2).*(data-a3);
        gb1 = net.apply_act_deriv(gint1,1).*(bsxfun(@plus,gb2*w2_',beta*sparse_grad'));
        % grad wrt w2 & w1
        if net.net_params.weight_tie
            gw1 = [];
            gw2 = a2'*gb2+gb1'*data;
        else
            gw1 = data'*gb1;
            gw2 = a2'*gb2;
        end
        gb1 = sum(gb1,1)';
        gb2 = sum(gb2,1)';
        grad = [gw1(:); gw2(:); gb1; gb2];

        % REG PENALTIES
        reg_pen_func = 0.5*lambda_w*sum(sum(w1_.^2))+...
                       0.5*lambda_w*sum(sum(w2_.^2))+...  % L2 penalty on weights
                       0.5*lambda_b1*sum(b1_.^2)+...
                       0.5*lambda_b2*sum(b2_.^2)+...      % L2 penalty on biases
                       beta*sparse_func;                 % sparsity penalty
        if net.net_params.weight_tie
            reg_pen_grad = [2*lambda_w*w2_(:); lambda_b1*b1_; lambda_b2*b2_];
        else
            reg_pen_grad = [lambda_w*w1_(:); lambda_w*w2_(:); lambda_b1*b1_; lambda_b2*b2_];
        end

        func = (func + reg_pen_func)/Z;
        grad = (grad + reg_pen_grad)/Z;
    
        end % internal function sgd
		
    end
    
    function fig_handle = disp_weights(net)
    % net.disp_weight()
    %
    % function that plots weights between hidden layer and output layer
    % (decoding weights)
    % INPUTS:
    %   none
    %
    % OUTPUTS:
    %   fig_handle: figure handle of output
    
    % check inputs
    
    fig_handle = figure;
    weights = net.w2';
    imagesc(weights,[-max(abs(weights(:))),max(abs(weights(:)))]);
    colormap(jet);
    
    end
    
    function fig_handle = disp_model(net,data)
    % net.disp_model()
    %
    % function that displays the hidden unit activation functions on top of
    % the distribution of generating signals (input)
    % INPUTS:
    %   data:       num_ex x num_io_nodes matrix of data
    %
    % OUTPUTS:
    %   fig_handle: figure handle of output
    % TODO:
    % add indx_tr
    
    % check inputs
    net.check_inputs(data);
    
    fig_handle = figure;
    if net.net_params.num_hid_nodes > 1
        num_cols = 2;
    else 
        num_cols = 1;
    end
    num_rows = ceil(net.net_params.num_hid_nodes/num_cols);
    
    % get input to hidden units
    [~,gint] = net.get_model_internals(data);
    m = min(gint{1}(:));        % get minimum value
    M = max(gint{1}(:));        % get maximum value
    edge = max([abs(m),abs(M)]);
    edges = linspace(-edge,edge,21);
    counts = histc(gint{1},edges,1);
    % scale counts by largest value of bins
    counts = bsxfun(@rdivide,counts,max(counts,[],1))*edges(end);
    
    % get hidden unit activation function plot points
    switch net.net_params.act_func_hid
        case 'lin'
            y = edges;
        case 'relu'
            y = edges;
            y(y<=0) = 0;
        case 'sigmoid'
            y = net.sigmoid(edges);
    end
    counter = 0;
    for n = 1:net.net_params.num_hid_nodes
        counter = counter+1;
        subplot(num_rows,num_cols,counter)
        bar(edges,counts(:,n),'histc')
        hold on;
        plot(edges,y,'LineWidth',2,'Color','k')
        xlim([-edge,edge]);
    end
    end
    
end
%% ********************  hidden methods ********************************
methods (Hidden)

    function sig = apply_act_func(net,sig,layer)
    % sig = net.apply_act_func(sig,layer)
    % 
	% internal function that applies activation function of the nodes
    % in a specified layer to a given input
    % INPUTS:
    %   sig:    in_layer_nodes x 1 vector
    %   layer:  1 or 2 for hidden or output layer, respectively
    %
	% OUTPUTS:
    %   out:    linearly weighted combo of input passed through
    %           activation function

    if layer == 1
        switch net.net_params.act_func_hid
            case 'lin'
            case 'relu'
                sig = max(0,sig);
            case 'sigmoid'
                sig = Autoencoder.sigmoid(sig);
        end
    elseif layer == 2
        switch net.net_params.act_func_out
            case 'lin'
            case 'relu'
                sig = max(0,sig);
            case 'sigmoid'
                sig = Autoencoder.sigmoid(sig);
        end
    end

    end
    
    function sig = apply_act_deriv(net,sig,layer)
    % sig = net.apply_act_deriv(sig,layer)
	%
    % internal function that calculates the gradient of the  activation 
    % function of the nodes in a specified layer to a given input
    % INPUTS:
    %     sig:      in_layer_nodes x 1 vector
    %     layer:    1 or 2 for hidden or output layer, respectively
	%    
	% OUTPUTS:
    %     sig:      linearly weighted combo of input passed through
    %               activation function

    if layer == 1
        switch net.net_params.act_func_hid
            case 'lin'
                sig = ones(size(sig));
            case 'relu'
				sig(sig<=0) = 0;
				sig(sig>0)  = 1;
            case 'sigmoid'
                sig = net.sigmoid_grad(sig);
        end
    elseif layer == 2
        switch net.net_params.act_func_out
            case 'lin'
                sig = ones(size(sig));
            case 'relu'
				sig(sig<=0) = 0;
				sig(sig>0)  = 1;
            case 'sigmoid'
                sig = net.sigmoid_grad(sig);
        end
    end

    end
    
    function [func, grad] = get_sparse_penalty(net,a2)
    % [func, grad] = net.get_sparse_penalty(a2)
	%
    % internal function that evaluates sparsity term in the cost function
    % and its gradient; enforces sparsity in activity of hidden layer
    % INPUTS:
    %   a2:     activation of the second (hidden) layer
	%    
	% OUTPUTS:
    %   func:   function value (scalar)
    %   grad:   gradient (column vector)
    
    switch net.net_params.act_func_hid
        case 'lin'
            % optimized
            func = sum(a2(:));
            % old
            % func = sum(mean(a2,1));
            grad = ones(size(a2,2),1);
        case 'relu'
            % optimized
            func = sum(a2(:));
            grad = sign(sum(abs(a2),1))';
            % old
            % func = sum(mean(abs(a2),1));
            % grad = sign(mean(abs(a2),1))';
        case 'sigmoid'
            p = net.reg_params.kl_param;
            pj = sum(a2,1);
            func = sum(p*log(p./pj)+(1-p)*log((1-p)./(1-pj)));
            grad = (-p./pj+(1-p)./(1-pj))';
    end
    end
    
	function check_inputs(net,data)
	% net.check_inputs(data)
	%
	% internal function that checks if the input data dims are
	% consistent with the specified network structure
	% INPUTS:
	% 	data:	data vector
	%
	% OUTPUTS:
	%	none; throws error flag if parameters are not consistent

	assert(size(data,2)==net.data_params.num_io_nodes,...
		'input data does not have proper number of parameters')
	end
		
    function check_weight_dims(net,weight_vec)
    % net.check_weight_dims(weights)
	%    
	% internal function that checks if the weight dimensions are
    % consistent with the specified network structure
    % INPUTS:
    %   weight_vec: 	vector of weights
	%    
	% OUTPUTS:
    %   none; throws error flag if parameters are not consistent
    
    if net.net_params.weight_tie
        d = net.net_params.num_io_nodes*net.net_params.num_hid_nodes + ...
            net.net_params.num_io_nodes+net.net_params.num_hid_nodes;
        assert(length(weight_vec) == d,'weight vector inconsistent with network structure')
    else
        d = net.net_params.num_io_nodes*net.net_params.num_hid_nodes + ...
            net.net_params.num_io_nodes+net.net_params.num_hid_nodes;
        assert(length(weight_vec) == d,'weight vector inconsistent with network structure')
    end
    end

end
%% ********************  static methods ********************************
methods (Static)
    
    function init_params = set_init_params(data)
    % init_params = Autoencoder.set_init_params(data)
	%    
	% takes data matrix and extracts relevant information for initializing
    % a network model. 
    % INPUTS:
    %   data:               number of examples x number of predictors data
    %                       matrix
    % OUTPUTS:
    %   init_params
    %       num_ex:         number of data points
    %       num_io_nodes:   dimensionality of input
    
    [init_params.num_ex, init_params.num_io_nodes] = size(data); 
        
    end

	function optim_params = set_init_optim_params()
	% optim_params = Autoencoder.set_init_optim_params();
	%
	% function that sets default optimization parameters for the various
	% optimization routines

	% both matlab and mark schmidt options
	optim_params.maxIter = 1000;
    optim_params.optTol = 1e-10;      				% termination tolerance on first order optimality (max(abs(grad))
    optim_params.progTol = 1e-16;     				% termination tolerance on function/parameter values
	optim_params.TolX = 1e-10;						% termination tolerance on function values
	optim_params.TolFun = 1e-7;						% termination tolerance on function values
    optim_params.maxFunEvals = 2*optim_params.maxIter;% mostly for sd in minFunc

	% just mark schmidt options
	optim_params.Method = 'lbfgs';

	% just matlab options
	optim_params.Algorithm = 'quasi-newton';
	optim_params.HessUpdate = 'steepdesc';			% default is bfgs, which is incredibly slow
    optim_params.GradObj = 'on';
    optim_params.DerivativeCheck = 'off';
    optim_params.numDiff = 0;

    end

    function [w1,w2,b1,b2,init_params] = ...
    set_init_weights_stat(init_weights,num_io_nodes,num_hid_nodes,weight_tie)
    % [w1,w2,b1,b2,init_params] = Autoencoder.set_init_weights(
    %                   init_weights,num_io_nodes,num_hid_nodes,weight_tie)
    % 
    % static function that initializes weights and sets init_params
    % structure based on input. Called from the constructor (which is why
    % its a static method) and from the non-static method set_init_weights
    % INPUTS:
    %   init_weights:   either a string specifying type of random 
    %                   initialization for weights ('gaussian' or 
    %                   'uniform') or a weight vector of appropriate length
    %   num_io_nodes:   number of input/output nodes in network
    %   num_hid_nodes:  number of hidden nodes in network
    %   weight_tie:     0 or 1 specifying whether or not encoding and
    %                   decoding weights are the same (1) or different (0)
    % 
    % OUTPUTS:
    %   w1:             num_io_nodes x num_hid_nodes weight matrix
    %   w2:             num_hid_nodes x num_io_nodes weight matrix
    %   b1:             num_hid_nodes x 1 bias vector for hidden layer
    %   b2:             num_io_nodes x 1 bias vector for output layer
    %   init_params:    struct specifying init_weights and rng_state, if
    %                   applicable
    
    if ischar(init_weights)
        % randomly initialize weights; start biases off at 0
        init_params.init_weights = lower(init_weights);
        init_params.rng_state = rng();
        % create filter
        s = 0.5;
        switch lower(init_weights)
            case 'gaussian'
                w1 = s*randn(num_io_nodes,num_hid_nodes);
                if weight_tie
                    w2 = w1';
                else
                    w2 = s*randn(num_hid_nodes,num_io_nodes);
                end
            case 'uniform'
                r = 4*sqrt(6) / sqrt(num_hid_nodes+num_io_nodes+1);   % we'll choose weights uniformly from the interval [-r, r]
                w1 = rand(num_io_nodes,num_hid_nodes)*2*r-r;
                if weight_tie
                    w2 = w1';
                else
                    w2 = rand(num_hid_nodes,num_io_nodes)*2*r-r; 
                end
        end
        b1 = zeros(num_hid_nodes,1);
        b2 = zeros(num_io_nodes,1);        
    elseif isvector(init_weights)
        % use 'init_weights' to initialize weights
        % split up initial vector into appropriate matrices and vectors
        if weight_tie % no w1 matrix; just [---w1----|-b1-|-b2-]
            assert(length(init_weights)==(num_hid_nodes*num_io_nodes+num_hid_nodes+num_io_nodes),...
            'init_weights vector has improper size')
            % get cut points
            cut1 = num_hid_nodes*num_io_nodes;
            cut2 = cut1+num_hid_nodes;
            % get weights
            w1 = init_weights(1:cut1);
            w1 = reshape(w1,num_io_nodes,num_hid_nodes);
            w2 = w1';
            % get biases
            b1 = init_weights(cut1+1:cut2);
            b2 = init_weights(cut2+1:end);
        else
            assert(length(init_weights)==(2*num_hid_nodes*num_io_nodes+num_hid_nodes+num_io_nodes),...
            'init_weights vector has improper size')
            % get cut points
            cut1 = num_hid_nodes*num_io_nodes;
            cut2 = 2*cut1;
            cut3 = cut2+num_hid_nodes;
            % get weights
            w1 = init_weights(1:cut1);
            w1 = reshape(w1,num_io_nodes,num_hid_nodes);
            w2 = init_weights(cut1+1:cut2);
            w2 = reshape(w2,num_hid_nodes,num_io_nodes);
            % get biases
            b1 = init_weights(cut2+1:cut3);
            b2 = init_weights(cut3+1:end);
        end
        init_params.rng_state = NaN;
    else
        warning('init_weights must be a string or a vector')
    end
    end
    
end
%% ******************** static hidden methods *****************************
methods (Static, Hidden)
    
    function y = sigmoid(x)
    % evaluate sigmoid units
        y = 1 ./ (1 + exp(-x));
    end

    function y = sigmoid_grad(x)
    % evaluates gradient of sigmoid units
        y = exp(-x) ./ (1 + exp(-x)).^2;
    end
    
end 
    
end



