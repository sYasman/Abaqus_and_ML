function [J grad] = nnCostFunction(nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, ...
                                   X, y, lambda)
%NNCOSTFUNCTION Implements the neural network cost function for a two layer
%neural network which performs classification
%   [J grad] = NNCOSTFUNCTON(nn_params, hidden_layer_size, num_labels, ...
%   X, y, lambda) computes the cost and gradient of the neural network. The
%   parameters for the neural network are "unrolled" into the vector
%   nn_params and need to be converted back into the weight matrices. 
% 
%   The returned parameter grad should be a "unrolled" vector of the
%   partial derivatives of the neural network.
%

% Reshape nn_params back into the parameters Theta1 and Theta2, the weight matrices
% for our 2 layer neural network
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
                 num_labels, (hidden_layer_size + 1));

% Setup some useful variables
m = size(X, 1);
         
% You need to return the following variables correctly 
J = 0;
Theta1_grad = zeros(size(Theta1));
Theta2_grad = zeros(size(Theta2));

% ====================== YOUR CODE HERE ======================
% Instructions: You should complete the code by working through the
%               following parts.
%
% Part 1: Feedforward the neural network and return the cost in the
%         variable J. After implementing Part 1, you can verify that your
%         cost function computation is correct by verifying the cost
%         computed in ex4.m
%
% Part 2: Implement the backpropagation algorithm to compute the gradients
%         Theta1_grad and Theta2_grad. You should return the partial derivatives of
%         the cost function with respect to Theta1 and Theta2 in Theta1_grad and
%         Theta2_grad, respectively. After implementing Part 2, you can check
%         that your implementation is correct by running checkNNGradients
%
%         Note: The vector y passed into the function is a vector of labels
%               containing values from 1..K. You need to map this vector into a 
%               binary vector of 1's and 0's to be used with the neural network
%               cost function.
%
%         Hint: We recommend implementing backpropagation using a for-loop
%               over the training examples if you are implementing it for the 
%               first time.
%
% Part 3: Implement regularization with the cost function and gradients.
%
%         Hint: You can implement this around the code for
%               backpropagation. That is, you can compute the gradients for
%               the regularization separately and then add them to Theta1_grad
%               and Theta2_grad from Part 2.
%

# one vs all
y_all = zeros(m, num_labels);
for i=1:m
  y_all(i, y(i)) = 1;
endfor  

# add bias unit to X matrix
X = [ones(m, 1) X];

# z2 = Theta1 *a0', a0=X
z2 = Theta1 * X'; 
# a2 = g(z2)
a2 = sigmoid( z2 );
# add bias term
#a2 = [1, ....; 1, ....; ...]
a2 = [ones(1, m); a2];

# z3 = Theta2 *a2'
z3 = Theta2 * a2; 
# a3 = g(z3)
# only in hidden layers add bias term
a3 = sigmoid( z3 );

# Typical cost function
J = sum((-1 / m) * sum(sum(y_all .* log(a3') + (1 - y_all) .* log(1 - a3'))));

# Remove bias terms
# Theta1 (i) = [BIAS TERM, v, v, v, v.. ]
Theta1_no_bias = Theta1(:, 2:end);
Theta2_no_bias = Theta2(:, 2:end);

# Adding regularization
J += lambda/2/m * ( sum(sum(Theta1_no_bias.^2)) + sum(sum(Theta2_no_bias.^2)) );

Delta_1 = 0;
Delta_2 = 0;
% TODO: CORRECT
for t=1:m
  # forward propagation
  a_1 = X(t,:)'; # 1. is already included in line 72
  z_2 = Theta1 * a_1; 
  a_2 = [1; sigmoid( z_2 )];
  z_3 = Theta2 * a_2; 
  a_3 = sigmoid( z_3 );
  
  # backward propagation
  d_3 = a_3 - y_all(t,:)';
  d_2 = (Theta2_no_bias'*d_3) .* sigmoidGradient(z_2);
  
  # D matrices
  Delta_2 += d_3 * a_2';
  Delta_1 += d_2 * a_1';
endfor

# Normal (NOT REGULARIZED) version
Theta1_grad = Delta_1 / m;
Theta2_grad = Delta_2 / m;

# Regularized version
Theta1_grad(:,2:end) += lambda / m * Theta1_no_bias;
Theta2_grad(:,2:end) += lambda / m * Theta2_no_bias;
% -------------------------------------------------------------

% =========================================================================

% Unroll gradients
grad = [Theta1_grad(:) ; Theta2_grad(:)];


end
