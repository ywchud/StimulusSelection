function bool=PistonComb(PistonArray,Piston)
gap=2;  %check PistonArray(i) with  Piston(i*gap)
if ~islogical(PistonArray)
    PistonArray=logical(PistonArray);
end
PistonArray=toColumn(PistonArray)';
bool=[];
for i=1:length(Piston)
    bool=[bool; (sum(~(~isnan(Piston(i,1*gap:1*gap:end))==PistonArray))==0)];           %returns 1 if piston(i)==user defined PistonArray
end
bool=logical(bool);
end