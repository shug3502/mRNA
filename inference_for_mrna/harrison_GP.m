function [fBarStar, CVfStar] = harrison_GP(Xtrain,Xtest,y,k,sn)

%created 17/12/15 JH
%last edit 17/12/15 JH
%
%Inputs Xtrain training data
%       Xtest testing data? or possibly grid of pts to evaluate on
%       y train data outputs
%       k covariance function
%       sn variance of assumed additive noise 

%See Rasmussen and Wiliams 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = k(Xtrain,Xtrain);
Kstar = k(Xtest,Xtrain)';
L = chol(K+sn*eye(size(K)),'lower');
alpha = L'\(L\y);
fBarStar = Kstar'*alpha; %mean of f star (predictkive distribution)
v = L\Kstar;
CVfStar = k(Xtest,Xtest) - v'*v;   %covariance of f star (predictive distribution)

