% demonstrate usage of regression
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2013-10-16.
clear all, close all


%% SAY WHICH CODE WE WISH TO EXERCISE
id = [1,4];

seed = 197; randn('seed',seed), rand('seed',seed)


data = csvread('train.csv');
prop = 0.8;
ntr = size(data,1)*prop; nte = size(data,1)*(1-prop);                        % number of training and test points
for response_ind = 1:7

xtr = data(1:ntr,8:end);
sn = 0.2;
ytr = data(1:ntr,response_ind);
xte = data((1+ntr):end,8:end);
yte = data((1+ntr):end,response_ind);

cov = {@covSEiso}; sf = 1; ell = 0.4;                             % setup the GP
hyp0.cov  = log([ell;sf]);
mean = {@meanSum,{@meanLinear,@meanConst}}; a = 0; b = 1;       % m(x) = a*x+b
hyp0.mean = [a*zeros(size(data,2)-7,1);b];

lik_list = {'likGauss','likLaplace','likSech2','likT'};   % possible likelihoods
inf_list = {'infExact','infLaplace','infEP','infVB','infKL'};   % inference algs

Ncg = 50;                                   % number of conjugate gradient steps
sdscale = 0.5;                  % how many sd wide should the error bars become?
col = {'k',[.8,0,0],[0,.5,0],'b',[0,.75,.75],[.7,0,.5]};                % colors
ymu{1} = yte; ys2{1} = sn^2; nlZ(1) = -Inf;
for i=1:size(id,1)
  lik = lik_list{id(i,1)};                                % setup the likelihood
  if strcmp(lik,'likT')
    nu = 4;
    hyp0.lik  = log([nu-1;sqrt((nu-2)/nu)*sn]);
  else
    hyp0.lik  = log(sn);
  end
  inf = inf_list{id(i,2)};
  fprintf('OPT: %s/%s\n',lik_list{id(i,1)},inf_list{id(i,2)})
  if Ncg==0
    hyp = hyp0;
  else
    hyp = minimize(hyp0,'gp', -Ncg, inf, mean, cov, lik, xtr, ytr); % opt hypers
  end
  [ymu{i+1}, ys2{i+1}] = gp(hyp, inf, mean, cov, lik, xtr, ytr, xte);  % predict
  [nlZ(i+1)] = gp(hyp, inf, mean, cov, lik, xtr, ytr);
end

figure, hold on
%for i=1:size(id,1)+1
%  plot(xte,ymu{i},'Color',col{i},'LineWidth',2)
%  if i==1
%    leg = {'function'};
%  else
%    leg{end+1} = sprintf('%s/%s -lZ=%1.2f',...
%                                lik_list{id(i-1,1)},inf_list{id(i-1,2)},nlZ(i));
%  end
%end
%for i=1:size(id,1)+1
%  ysd = sdscale*sqrt(ys2{i});
%  fill([xte;flipud(xte)],[ymu{i}+ysd;flipud(ymu{i}-ysd)],...
%       col{i},'EdgeColor',col{i},'FaceAlpha',0.1,'EdgeAlpha',0.3);
%end
%for i=1:size(id,1)+1, plot(xte,ymu{i},'Color',col{i},'LineWidth',2), end
%plot(xtr,ytr,'k+'), plot(xtr,ytr,'ko'), legend(leg)

plot(ymu{1},ymu{2},'bo'),xlabel('y actual'),ylabel('y predicted');
print(sprintf('gpregress%d',response_ind),'-dpng');
end
