function [y,y_um]=CSDplot(Exp,Neuron,CSD_data)
getNormDepth=1;
% landMark_um=[90 200 360 600 700];

%% hehe ispline CSD code yoinked (from CSDplotter)
G=findgroups(Exp.KSort.x);
diam=5.0000e-04;
cond=0.3;cond_top=0.3;
gauss_sigma=1.0000e-04;
filter_range=5.0000e-04;
dt=1/Exp.KSort.LFP_fs*1000;
h=figure('WindowState','maximized');  
h_CSDonly=figure('WindowState','maximized');  
shk_i_ypos=cell(1,Exp.KSort.ShankN);
for i=1: Exp.KSort.ShankN
    topDepth=Exp.KSort.top(i);
    offset=-(-topDepth-Exp.KSort.ChGap*[0:Exp.KSort.ChN(i)-1]);
%     el_pos=offset;

    % not using Exp.KSort.y cos it could be inverted, we know CSD_data
    % starts from top to bot
%     el_pos=(topDepth+(Exp.KSort.ChGap:Exp.KSort.ChGap:Exp.KSort.ChN(i)*Exp.KSort.ChGap))*1e-6;
    el_pos=(Exp.KSort.ChGap:Exp.KSort.ChGap:Exp.KSort.ChN(i)*Exp.KSort.ChGap)*1e-6;

    data=CSD_data{i};
    % compute spline iCSD:
    Fcs = F_cubic_spline(el_pos,diam,cond,cond_top);
    [zs,CSD_cs] = make_cubic_splines(el_pos,data,Fcs);
    if gauss_sigma~=0 %filter iCSD
      [zs,CSD_cs]=gaussian_filtering(zs,CSD_cs,gauss_sigma,filter_range);  
    end
    % plot CSD
    figure(h)
    subplot(1,Exp.KSort.ShankN*2,i*2-1);
    plot_CSD(CSD_cs,zs,dt,1,0)
    
    observedPosition=-round(zs*1e6+topDepth-Exp.KSort.ChGap);
    
    y=yticks;   
    yy=observedPosition(y);
    
    set(gca,'YTickLabel',yy)
    shk_i_ypos{i}=observedPosition;
    % plot Neuron Onset
    subplot(1,Exp.KSort.ShankN*2,i*2)
    DepthXPSTHOnset(Exp,Neuron,i)
%     ylim([zs(end) zs(1)]*-1e6-topDepth)
    ylim(-round([zs(end) zs(1)]*1e6+topDepth-Exp.KSort.ChGap))
    figure(h_CSDonly)
    subplot(1,Exp.KSort.ShankN,i);
    plot_CSD(CSD_cs,zs,dt,1,0)
    y=yticks;
    yy=-round(zs(y)*1e6+topDepth-Exp.KSort.ChGap);
    set(gca,'YTickLabel',yy)

end 
filename=fullfile(Exp.Path.save,'touchCSD.jpg');
saveas(gcf,filename)

%%
if getNormDepth
    figure(h_CSDonly)
    
    Nlandmarks=5;
    y=zeros(Nlandmarks,Exp.KSort.ShankN);
    y_um=zeros(Nlandmarks,Exp.KSort.ShankN);
    for k=1:Exp.KSort.ShankN
        for i=1:Nlandmarks
            str=sprintf('Shank %d landmark %d/%d',k,i,Nlandmarks);
            title(str)
            [~,y(i,k)] = ginputWhite(1);
            hold on
            yline(y(i,k),'w');
        end
        y_um(:,k)=shk_i_ypos{k}(ceil(y(:,k)));
    end
    
    
    
else 
    y=nan;
    y_um=nan;
end

end

