function [C, sigma] = dataset3Params(X, y, Xval, yval)
%DATASET3PARAMS returns your choice of C and sigma for Part 3 of the exercise
%where you select the optimal (C, sigma) learning parameters to use for SVM
%with RBF kernel
%   [C, sigma] = DATASET3PARAMS(X, y, Xval, yval) returns your choice of C and 
%   sigma. You should complete this function to return the optimal C and 
%   sigma based on a cross-validation set.
%

% You need to return the following variables correctly.
C = 1;
sigma = 0.3;

% ====================== YOUR CODE HERE ======================
% Instructions: Fill in this function to return the optimal C and sigma
%               learning parameters found using the cross validation set.
%               You can use svmPredict to predict the labels on the cross
%               validation set. For example, 
%                   predictions = svmPredict(model, Xval);
%               will return the predictions on the cross validation set.
%
%  Note: You can compute the prediction error using 
%        mean(double(predictions ~= yval))
%

C_test = [0.01, 0.03, 0.10, 0.30, 1.00, 3.00, 10.00, 30.00];
sigma_test = [0.01, 0.03, 0.10, 0.30, 1.00, 3.00, 10.00, 30.00];


##results = [];
##values = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30];
##for c = values
##   for s = values
##      model= svmTrain(X, y, c, @(x1, x2) gaussianKernel(x1, x2, s));
##	  predict = svmPredict(model, Xval);
##	  error = mean(double(predict ~= yval));
##	  results = [results; c, s, error];
##   end
##end
##
##[minError, minIndex] = min(results(:, 3));
##C = results(minIndex, 1);
##sigma = results(minIndex, 2);


x1 = X(:,1);
x2 = X(:,2);
err_best = 1e6;
i_best = -99;
j_best = -99;
res = [];
for i=1:length(C_test)
  for j=1:length(sigma_test)
    model = svmTrain(X, y, C_test(i), @(x1, x2) gaussianKernel(x1, x2, sigma_test(j))); 
    predictions = svmPredict(model, Xval);
    err = mean(double(predictions ~= yval))
    res = [res; C_test(i), sigma_test(j), err];
  endfor
endfor

[minError, minIndex] = min(res(:, 3));
C = res(minIndex, 1);
sigma = res(minIndex, 2);





% =========================================================================

end
