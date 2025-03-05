function Whisker=findWhiskerbase(Exp,Whisker)

    disp('Obtaining whiskerbase...')
    for i=1:length(Whisker)
        SampleSize=1;
        if isempty(Whisker(i).whiskerbaseX) || isempty(Whisker(i).whiskerbaseY)
            Whisker(i)=Whisker(i).Whiskerbase(Exp,SampleSize);  
        else
            x = input(sprintf('Whiskerbase already exists for W%d. Do you wish to redo it?(y/n)',Whisker(i).W_id));
            if x(1)=='y' || x(1)=='Y'
                fprintf('Finding Whiskerbase for whisker %d\n',i)
                Whisker(i)=Whisker(i).Whiskerbase(Exp,SampleSize);
            else
                continue;
            end
        end
    end
    disp('Whiskerbase obtained') 

end