function dprime_plot(D_stim,PSTHrange,map)
%plot D_stim (Neuron*PSTHrange)

H=D_stim;
if ~exist('map','var')
    map = [[[1:-1/31:0]' zeros(32,1) zeros(32,1)]; [zeros(32,1) [0:1/31:1]' zeros(32,1)]];
end
h1 = heatmap(H,'GridVisible','off','Colormap',map,'ColorLimits',[-1.1 1.1]);
time_label = PSTHrange(1:end);
temp = 1:1:length(time_label);
temp(1:10:length(time_label)) = NaN;
time_label(rmmissing(temp)) = NaN;
h1.XDisplayLabels = time_label;
h1.YDisplayLabels = NaN(size(H,1),1);
h1.YLabel = ['neurons = ' num2str(size(H,1))];
h1.XLabel = 'time (s)';
% Damn code just to get a line on this heatmap
% Get underlying axis handle
origState = warning('query', 'MATLAB:structOnObject');
cleanup = onCleanup(@()warning(origState));
warning('off','MATLAB:structOnObject')
S = struct(h1); % Undocumented
ax = S.Axes;    % Undocumented
clear('cleanup')
% Remove grids
h1.GridVisible = 'off';
end