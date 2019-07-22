function [x,f,output] = sgd(obj_fun,x,options,data,varargin)
% [x,f,output] = minFunc(funObj,x0,options)
%
% Implementation of stochastic gradient descent method
% INPUTS:
%   obj_fun:        function handle that returns a function value
%                   and gradient for given initial parameter values x and
%                   given data d
%   x:              initial parameter values
%   options:        struct containing optimizer options
%       batch_size  number of data points to use on each step
%       maxIter     maximum number of passes through data
%       Display     'off' to suppress output
%                   'iter' to display info for each batch of data in every pass through data
%                   'batch' to display info after each pass through data
%       monitor     'off' to suppress saving output
%                   'iter' to save info for each batch of data in every pass through data
%                   'batch' to save info after each pass through the data
%					'both' is also supported
%   data:           num_examples x dim data matrix
%   varargin:       additional input parameters to obj_fun
%
% OUTPUTS:
%   x:              minimizer of objective function
%   f:              minimized function value
%   output:         struct containing optimization info
%       msg         message specifying exit condition        
%       trace.fval  function values across iterations
%       trace.norm_grad
%       trace.firstorderopt
%
% Reference:
% Author: Matt Whiteway
%   02/09/16

% get data info
N = size(data,1);

% define default params
batch_size = 10;
maxIter = 500;
Display = 'batch';
monitor = 'batch';

% parse input options
if isfield(options,'batch_size')
    batch_size = options.batch_size;
end
if isfield(options,'maxIter')
    maxIter = options.maxIter;
end
if isfield(options,'Display')
	Display = options.Display;
end
if isfield(options,'monitor')
	monitor = options.monitor;
end
    
% monitor progress
if strcmp(monitor,'iter') || strcmp(monitor,'both')
    trace.fval			= zeros(0,1);		% function value 
    trace.dirderiv      = zeros(0,1);		% directional derivative
    trace.firstorderopt = zeros(0,1);		% first order optimality
end
if strcmp(monitor,'batch') || strcmp(monitor,'both')
    trace.fval_batch          = zeros(0,1); % function value 
    trace.dirderiv_batch      = zeros(0,1); % directional derivative
    trace.firstorderopt_batch = zeros(0,1); % first order optimality
end

% monitor progress even if not saving for exit conditions
temp_fval_old = Inf;

% print output log
if strcmp(Display,'iter')
    fprintf('%10s %15s %15s %15s\n','Iteration','Function Val','Grad Norm','1st Order Opt');
elseif strcmp(Display,'batch')
    fprintf('%10s %15s %15s %15s\n','Batch Iter','Function Val','Grad Norm','1st Order Opt');
end
    
% define parameters; these definitely need to be tuned for particular
% appliations; step size is defined as a/(sqrt(b+j)), where j indexes the
% iteration; guarantees that step size will decrease
a = 0.01;
b = 0;

eta = 0;				% for momentum term on the parameter vector; eta > 1 is equivalent to "short memory" averaging 
dw = 0;					% for momentum term on the gradient
% G = zeros(size(x));	% for adaptive gradient	
rng(0);					% for reproducibility

% maxIter is the number of passes through the entire dataset; within each
% pass there will be 'batch_iters' iterations of steepest descent
batch_iters = floor(N/batch_size);
for i = 1:maxIter

    % randomly permute entries in data matrix
    w = randperm(N);
    data = data(w,:);

    % monitor progress even if not saving for exit conditions
    temp_fval = zeros(batch_iters,1);	% function value
    temp_dd   = zeros(batch_iters,1);	% directional derivative
    temp_foo  = zeros(batch_iters,1);	% first-order optimality
    
	for j = 1:batch_iters
        
        % pick subset of data for each iteration
        d = data((j-1)*batch_size+(1:batch_size),:);

        % compute descent direction
        [~,g] = obj_fun(x,d,varargin{:});

		% below are different options to play with. Again, different
		% strategies will work better for different problems. This is
		% definitely in need of cleaning up
		
        % 1. momentum on gradient
        dw = 0.95*g+0.05*dw;
        x = x - a/sqrt(b+j)*dw;

        % 2. adaptive gradient
%         G = G+g.*g;
%         x = x - a/(sqrt(G)*(b+j))*g;

        % 3. take step using averaging and decreasing step size
%         if eta > 0
%             x0 = x - a/sqrt(b+j)*g;
%             x = j/(j+eta)*x + eta/(j+eta)*x0;
%         else
%             x = x - a/sqrt(b+j)*g;
%         end

        % update diagnostics
        [f,g] = obj_fun(x,d,varargin{:});
        temp_fval(j) = f;				% calculate function value
        temp_dd(j) = g'*g;				% calculate directional derivative
        temp_foo(j) = max(abs(g));		% calculate first-order optimality
        
		% update iter trace
        if strcmp(monitor,'iter') || strcmp(monitor,'both')
            trace.fval(end+1,1)			= temp_fval(j);
            trace.dirderiv(end+1,1)		= temp_dd(j);
            trace.firstorderopt(end+1,1)= temp_foo(j);
        end
        
        % output iter info
        if strcmp(Display,'iter')
            fprintf('%10d %15.5e %15.5e %15.5e\n',...
            (i-1)*batch_iters+j,temp_fval(j),temp_dd(j),temp_foo(j));
        end

	end % end batch loop
    
	% evaluate objective function on all data
	[f,g] = obj_fun(x,data);

    % update batch trace
    if strcmp(monitor,'batch') || strcmp(monitor,'both')
		trace.fval_batch(end+1,1)			= f;
		trace.dirderiv_batch(end+1,1)		= g'*g;
		trace.firstorderopt_batch(end+1,1)	= max(abs(g));
    end
    
    % output batch info
    if strcmp(Display,'batch')
		fprintf('%10d %15.5e %15.5e %15.5e\n',i,f,g'*g,max(abs(g)));
	end
    
	% Need to figure out better exit conditions
    % Check for lack of progress
%     if i > 1
%         rel_chg = temp_fval_old-mean(temp_fval);
%         if rel_chg < 0 && abs(rel_chg) < 1e-5
%             msg = 'Function value increase from previous pass';
%             x = x_old;
%             break
%         elseif abs(rel_chg) < 1e-5
%             msg = 'Function Value changing by less than progTol';
%             break
%         end
%     end

    % Check for going over iteration/evaluation limit
    if i == maxIter
        msg = 'Reached Maximum Number of Batch Iterations';
        break;
    end

    % save for next iteration
    temp_fval_old = mean(temp_fval);
    x_old = x;
    
end % end iters

if strcmp(Display,'iter') || strcmp(Display,'batch')
    fprintf('%s\n',msg);
end

output.msg = msg;
output.trace = trace;

