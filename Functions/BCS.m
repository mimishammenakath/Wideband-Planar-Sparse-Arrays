% Sparse Bayesian Learning (SBL) using Relevance Vector Machine (RVM)
% Assume we have training data X_train and corresponding labels y_train
function [weights] = ...
    BCS(PHI,t,beta,maxIts,delta)
% Initialize hyperparameters
% alpha = 1e-6;  % Initial precision
% beta = 1e-6;   % Initial noise precision

[N, M] = size(PHI);  % Number of samples and dimensionality

% Initialize parameters
w = pinv(PHI)*t;

% alpha=alpha*ones(M,1);

alpha=1./ (abs(w)).^2;
% B=beta *eye(N);

% max_iterations = 100;  % Maximum number of iterations

for iter = 1:maxIts
    fprintf(1,'SBL Algorithm  # iterations : %d \n',iter);
%     alpha_old = alpha;
    w_old = w;
    % Update posterior covariance
   A =diag(alpha)  +beta* (PHI'* PHI);
   

        A_inv=pinv(A);
    % Update posterior mean
    w =  beta*A_inv * (PHI' * t);
    
    % Update alpha
    alpha =1./( abs(w).^2);
    
    gamma = 1 - alpha .* diag(A_inv);
%     
    alpha =gamma./( abs(w).^2);
%     
%     % Update beta
    beta = (N - sum(gamma)) /norm(t - PHI *w)^2;
% %     
%     % Check convergence
% if  norm(w - w_old)/norm(w)< delta
%     break;
% end

figure(1)
plot(abs(w));
weights=w;
end


 
% % Predict on test data
% y_pred = X_test * w;
% 
% % Evaluation
% mse = mean((y_pred - y_test).^2);