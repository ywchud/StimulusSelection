function D_stim=dprime(varargin)
%get dprime for 2 matrixes or 2 mu,std vectors along second dim
if nargin==2 && isfield(varargin{1},'PSTHmat')% SS,SS2
    SS=varargin{1};
    SS2=varargin{2};
    x1_mu(1,:) = mean(SS.PSTHmat{1});
    x2_mu(1,:) = mean(SS2.PSTHmat{1});
    x1_stdsq(1,:) = (std(SS.PSTHmat{1})).^2;
    x2_stdsq(1,:) = (std(SS2.PSTHmat{1})).^2;
elseif nargin==2 && ismatrix(varargin{1})% PSTHmat1 PSTHmat2
    x1_mu(1,:) = mean(varargin{1});
    x2_mu(1,:) = mean(varargin{2});
    x1_stdsq(1,:) = (std(varargin{1})).^2;
    x2_stdsq(1,:) = (std(varargin{2})).^2;
elseif nargin==4 %x1_mu x2_mu x1_std x2_std
    x1_mu(1,:) =varargin{1};
    x2_mu(1,:) =varargin{2};
    x1_stdsq(1,:) = (varargin{3}).^2;
    x2_stdsq(1,:) = (varargin{4}).^2;
else
    error('wrong input format')
end

D_stim(1,:) = (x1_mu(1,:) - x2_mu(1,:))./sqrt((x1_stdsq(1,:)+x2_stdsq(1,:))/2);
% if any(isnan(D_stim(1,:))) && any(x1_mu(1,isnan(D_stim(1,:))) - x2_mu(1,isnan(D_stim(1,:)))~=0)
%     disp('pause')   
%     %we assign nan dprime as 0 only if all delta(mu) observations are 0s. Otherwise large mu difference should be respected
% end
if ~all(isnan(D_stim(1,:)))
    D_stim(1,(isnan(D_stim(1,:))))=0;
    %we assign nan dprime as 0 only if not all instances are nan. All nans
    %could mean Nsig neurons
end

end