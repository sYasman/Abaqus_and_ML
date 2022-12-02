function [theta, J_history] = gradientDescent(X, y, theta, alpha, num_iters)
%GRADIENTDESCENT Performs gradient descent to learn theta
%   theta = GRADIENTDESCENT(X, y, theta, alpha, num_iters) updates theta by 
%   taking num_iters gradient steps with learning rate alpha

% Initialize some useful values
m = length(y); % number of training examples
J_history = zeros(num_iters, 1);

for iter = 1:num_iters

    % ====================== YOUR CODE HERE ======================
    % Instructions: Perform a single gradient step on the parameter vector
    %               theta. 
    %
    % Hint: While debugging, it can be useful to print out the values
    %       of the cost function (computeCost) and gradient here.
    %
    h = theta(1)*X(:,1) + theta(2)*X(:,2);
    dJ_dTheta0 = 1/m * sum( (h-y) .* X(:,1));
    dJ_dTheta1 = 1/m * sum( (h-y) .* X(:,2));
    dJ_dTheta = [dJ_dTheta0; dJ_dTheta1];
    
    theta = theta - alpha*dJ_dTheta;
    % ============================================================

    % Save the cost J in every iteration    
    J_history(iter) = computeCost(X, y, theta);

end
end
