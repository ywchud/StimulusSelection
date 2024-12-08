function h=addline(z,value,color,lineStyle)
if isempty(color)
    color='k';
end
if isempty(lineStyle)
    lineStyle='-';
end
hold on 
    
for i=1:length(value)
    
    if z=='x'       
            h=line([value(i), value(i)],ylim(), 'LineWidth', 0.5, 'Color', color,'LineStyle',lineStyle);   
    elseif z=='y'    
            h=line(xlim(), [value(i),value(i)], 'LineWidth', 0.5, 'Color', color,'LineStyle',lineStyle);
    end
end

end