function [Exp,Neuron]=NeuralActivity(Exp,Neuron)
%% Determine neurons activity across trials
temp=nan(length(Neuron),1);
for i=1:length(Neuron)
    Neuron(i)=Neuron(i).findActiveTrials(Exp);
    temp(i)=mean(Neuron(i).ActiveTrials);
end
figure
histogram(temp,0:0.05:1)

temp1=find(temp<0.6);n=0;
while isempty(temp1)
    n=n+1;
    temp1=find(temp<0.6+n);
end
for i=1:1%length(temp1)
    figure
    Neuron(temp1(i)).findActiveTrials(Exp,1);
    temp(temp1(i))
end
figure
Mat=[];
for k=1:3
    subplot(1,3,k)
Ns=Neuron([Neuron.Shank]==k);disp(length(Ns))
[temp,I]=sort([Ns.Depth],'descend');
mat=false(length(Ns),Exp.TrN);
for i=1:length(Ns)
    mat(i,:)=Ns(I(i)).ActiveTrials;
end
imagesc(int8(~mat));
colormap(gray)
Mat=[Mat;mat];

hold on
II=find([Ns(I).Depth]<=-350,1,'first');
yline(II,'r');
II=find([Ns(I).Depth]<=-550,1,'first');
yline(II,'b');
end
filename=fullfile(Exp.Path.save,'NeuronActivityAcrossTrial.jpg');
saveas(gcf,filename)
% C = permute(Mat,[1 3 2]);
% C = reshape(C,[],size(mat,2),1);
figure,plot(smooth(mean(Mat))),hold on,ylabel('%ActiveNeurons'),xlabel('nth Trial')
temp=abs(diff(smooth(mean(Mat))));
yyaxis right
plot(temp)
filename=fullfile(Exp.Path.save,'NeuronActivityAcrossTrial_avg.jpg');
saveas(gcf,filename)
[~,Exp.nthTrialCutoff]=max(temp);
Exp.ActivityAcrossTrials=smooth(mean(Mat));
end