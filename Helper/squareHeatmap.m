function squareHeatmap(h)
% Temporarily change axis units 
originalUnits = h.Units;  % save original units (probaly normalized)
h.Units = 'centimeters';  % any unit that will result in squares
% Get number of rows & columns
sz = size(h.ColorData); 
% Change axis size & position;
originalPos = h.Position; 
% make axes square (not the table cells, just the axes)
h.Position(3:4) = min(h.Position(3:4))*[1,1]; 
if sz(1)>sz(2)
    % make the axis size more narrow and re-center
    h.Position(3) = h.Position(3)*(sz(2)/sz(1)); 
    h.Position(1) = (originalPos(1)+originalPos(3)/2)-(h.Position(3)/2);
else
    % make the axis size shorter and re-center
    h.Position(4) = h.Position(4)*(sz(1)/sz(2));
    h.Position(2) = (originalPos(2)+originalPos(4)/2)-(h.Position(4)/2);
end
% Return axis to original units
h.Units = originalUnits; 
end