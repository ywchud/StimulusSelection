function Exp=salvage_analog_from_script(Exp,pathfile)
% pathfile='E:\Darren\DarrenM3LMC_210319\DarrenM3LMCrecording_210319_164440\light_log.xlsx';
trash = readmatrix(pathfile);
trash=trash(~isnan(trash));
trash(trash==0)=nan;
if length(trash)~=Exp.TrN
    error('Incomplete trial info')
else
    Exp.Stim.Optic.f=trash;
    Exp.Stim.Optic.On=~isnan(trash);
    Exp.Stim.Optic.Delay=0.5;
end
end