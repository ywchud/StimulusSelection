function [p,chi2stat]=chi2_test(ns,Ns)

%test if L>=2 sampled populations have same counts of group1 and group2 with chi square
%ns is a (n=L) vector containing counts of group1s for each population
%Ns is a (n=L) vector containing total counts (i.e. group1 + group2) for
%each population

groupN=2; %for now groupN has to be 2. Need to change code to include more

if length(ns)==length(Ns)
    L=length(ns);
else
    error('vector 1 and 2 should have equal lengths\n')
end

       % Pooled estimate of proportion
       p0 = sum(ns) / sum(Ns);
       % Expected counts under H0 (null hypothesis)
       Expected_Group1s= Ns .* p0;
       Expected_Group2s= Ns - Expected_Group1s;
       
       % Chi-square test, by hand
       observed = [ns(:);Ns(:)-ns(:)];
       expected = [Expected_Group1s(:);Expected_Group2s(:)];
       chi2stat = sum((observed-expected).^2 ./ expected);
       %degree of freedom
       nu=(groupN-1)*(L-1);
       %p-value
       p = 1 - chi2cdf(chi2stat,nu);

       
       
       %---  2x2 only
%               p0 = (n1+n2) / (N1+N2);
%        % Expected counts under H0 (null hypothesis)
%        n10 = N1 * p0;
%        n20 = N2 * p0;
%        % Chi-square test, by hand
%        observed = [n1 N1-n1 n2 N2-n2];
%        expected = [n10 N1-n10 n20 N2-n20];
%        chi2stat = sum((observed-expected).^2 ./ expected);
%        p = 1 - chi2cdf(chi2stat,1);
end