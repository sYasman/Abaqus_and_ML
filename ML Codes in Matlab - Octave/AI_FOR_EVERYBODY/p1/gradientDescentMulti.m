function [theta, J_history] = gradientDescentMulti(X, y, theta, alpha, num_iters)
%GRADIENTDESCENTMULTI Performs gradient descent to learn theta
%   theta = GRADIENTDESCENTMULTI(x, y, theta, alpha, num_iters) updates theta by
%   taking num_iters gradient steps with learning rate alpha

% Initialize some useful values
m = length(y); % number of training examples
J_history = zeros(num_iters, 1);
n_theta = size(X,2);
n_der = n_theta;
dJ_dTheta = zeros(n_der, 1);
for iter = 1:num_iters

    % ====================== YOUR CODE HERE ======================
    % Instructions: Perform a single gradient step on the parameter vector
    %               theta. 
    %
    % Hint: While debugging, it can be useful to print out the values
    %       of the cost function (computeCostMulti) and gradient here.
    %

    dJ_dTheta = 1/m *X'*(X*theta - y);
    theta_old = theta;
    theta = theta - alpha* dJ_dTheta;
    J = computeCostMulti(X, y, theta);
    if J<1e-3,
      printf('The root is found %5f\n',theta)
      break
    endif
    % ============================================================

    % Save the cost J in every iteration    
     J_history(iter) = J;
end

end
