function [H,sampStat,CI,bootstrapStat]=bootstrapping3(x1,x2,x3,myStatistic,nReps)


CIrange = 90;
% nReps = 1000;
% n1 = 20;
% n2 = 100;

% myStatistic = @(x1,x2) var(x1)/var(x2);
sampStat = myStatistic(x1,x2,x3);
bootstrapStat = zeros(1,nReps);
for i=1:nReps
    resampX1 = x1(ceil(rand(size(x1))*length(x1)));
    resampX2 = x2(ceil(rand(size(x2))*length(x2)));
    resampX3 = x3(ceil(rand(size(x3))*length(x3)));
    bootstrapStat(i) = myStatistic(resampX1,resampX2,resampX3);
end
CI = prctile(bootstrapStat,[50-CIrange/2,50+CIrange/2]);
% fprintf('Ratio of variances: %5.2f\n',sampStat);
% fprintf('%d%% Confidence interval: [%5.2f,%5.2f]\n',CIrange,CI(1),CI(2));

if sampStat<=CI(2) && sampStat>=CI(1)
    H=0;
else
    H=1;
end

end