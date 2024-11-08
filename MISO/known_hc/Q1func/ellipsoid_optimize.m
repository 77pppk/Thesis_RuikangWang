function [obj_best, x_best] = ellipsoid_optimize(obj, obj_diff, cst, cst_diff, x0, max_abs, opt_err)
%% Optimization tool based on ellipsoid method (Y)
% Especially for problems with complex form so that cvx cannot recognize.

% 1. Explanations:
% 1.1 INPUT:
%        obj            --  function handle of objective function (to be minimized)
%        obj_diff       --  function handle of derivative of objective function
%        cst, cst_diff  --  cell structure including function handles for constraints
%        x0             --  initial (feasible) input
%        max_abs        --  maximum distance in each dimension from optimal point to initial point x0
%        opt_err        --  maximum tolenrant error to the optimal
% 1.2 OUTPUT:
%        obj_best       --  optimal value (minimized)
%        x_best         --  optimal point x

% default output for error occurance
obj_best = [];
x_best = [];

if size(obj,1) ~= 1 || size(obj,2) ~= 1 || size(obj_diff,1) ~= 1 || size(obj_diff,2) ~= 1 ...
        || size(cst,1) ~= 1 || size(cst_diff,1) ~= 1 || size(cst,2) ~= size(cst_diff,2)
    disp('Wrong objective format or wrong constraint format.')
    return;
end
max_abs = reshape( max_abs, size(max_abs,1) * size(max_abs,2), 1);
x0 = reshape( x0, size(x0,1) * size(x0,2), 1);
if length(max_abs) ~= length(x0)
    disp('Wrong initial input or wrong range.');
    return;
elseif opt_err <= 0
    disp('Negative threshold error.');
    return;
end

opt_cst_nr = size(cst, 2);

opt_n = length(x0);
opt_P = diag(opt_n * max_abs .^2);
opt_loopi = 0;
opt_objbest = inf;
opt_x = x0;
%opt_err=1e-6;

for opt_i = 1 : opt_cst_nr
    if cst{opt_i}(real(opt_x)) >= 0
        disp('Initial point is not strictly feasible.');
        return;
    end
end

while 1 %opt_loopi<5000 %%% note that for a large N a large number of iterations may be required.
    
    % mark the iteration -- objective iteration or constraint iteration
    opt_mark = 1;
    
    for opt_i = 1 : opt_cst_nr
        if cst{opt_i}(real(opt_x)) > 0
            opt_g = cst_diff{opt_i}(real(opt_x));
            opt_g = reshape(opt_g, opt_n, 1);
            opt_gn = opt_g / (opt_g' * opt_P * opt_g) ^0.5;
            opt_a = cst{opt_i}(real(opt_x)) / (opt_g' * opt_P * opt_g) ^0.5;
            opt_mark = 0;
            break;
        end
    end
    
    % objective iteration
    if opt_mark == 1
%         opt_x
        opt_obj = obj(real(opt_x));
        opt_g = obj_diff(real(opt_x));
        opt_g = reshape(opt_g, opt_n, 1);
        if opt_obj < opt_objbest
            opt_objbest = opt_obj;
            x_best = real(opt_x);
        end
        opt_gn = opt_g / (opt_g' * opt_P * opt_g) ^0.5;
        opt_a = (opt_obj - opt_objbest) / (opt_g' * opt_P * opt_g) ^0.5;
    end
    
    % stopping iteration checking
    if (opt_g' * opt_P * opt_g) ^0.5 <= opt_err && opt_mark == 1
        break;
    end
    if opt_a > 1 && opt_mark == 0 
        break;
    end
    
    if opt_a > 1
        disp( 'Error! Convexity is not guaranteed!' )
       break;
    end
    
    % updating
    opt_x = real(opt_x) - (1 + opt_n * opt_a) / (opt_n + 1) * opt_P * opt_gn;
    opt_P = opt_n ^2 / (opt_n ^2 - 1) * (1 - opt_a ^2) * (opt_P - 2 * (1 + opt_n * opt_a) / (opt_n + 1) / (opt_a + 1) * opt_P * (opt_gn * opt_gn') * opt_P);
    
    opt_loopi = opt_loopi + 1;
%     if isnan(opt_x(1))
%         break;
%     end
end
obj_best = opt_objbest;

end

