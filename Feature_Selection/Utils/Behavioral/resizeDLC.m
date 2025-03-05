function [Xmat,Ymat]=resizeDLC(Exp,Xmat,Ymat)
fprintf('Resizing Xmat Ymat to accomodate inactive pistons and/or disordered PWOI during DLC...\n')
for k=1:length(Xmat)
    temp=Xmat{k};temp1=Ymat{k};
    Xmat_padded=nan(size(temp,1),size(temp,2)+sum(~Exp.Stim.Piston.active)*Exp.DLC.LblPerW);
    Ymat_padded=nan(size(temp,1),size(temp,2)+sum(~Exp.Stim.Piston.active)*Exp.DLC.LblPerW);
    for i=1:Exp.Stim.Piston.tot
        if Exp.Stim.Piston.active(i)
            j=Exp.Stim.Piston.PWOI(i);
            Xmat_padded(:,Exp.DLC.LblPerW*(i-1)+1:Exp.DLC.LblPerW*i)=temp(:,Exp.DLC.LblPerW*(j-1)+1:Exp.DLC.LblPerW*j);
            Ymat_padded(:,Exp.DLC.LblPerW*(i-1)+1:Exp.DLC.LblPerW*i)=temp1(:,Exp.DLC.LblPerW*(j-1)+1:Exp.DLC.LblPerW*j);
        end
    end
    %snout
    if rem(size(Xmat_padded,2),Exp.DLC.LblPerW)==1
        Xmat_padded(:,end)=temp(:,end);
        Ymat_padded(:,end)=temp1(:,end);
    end
    Xmat{k}=Xmat_padded;
    Ymat{k}=Ymat_padded;
    
    clear Xmat_padded Ymat_padded
end
end