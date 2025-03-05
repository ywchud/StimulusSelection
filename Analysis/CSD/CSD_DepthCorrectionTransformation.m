function  [correctionMatrix,Vector2correct]=CSD_DepthCorrectionTransformation(y_um,landMark_um,Exp)
%size(y_um,1)==length(landMark_um)  i.e.n=5
%y_um is the observed depth of n landmarks * Nshank, this function attempts
%to interpolate them so the n landmarks align with given landMark_um

if size(y_um,1)~=length(landMark_um)
    error('Mismatch landmarks number')
end

ys_um=-1*y_um; %change to positive for interpolation

Vector2correct=[0:1:Exp.KSort.ChDepth+Exp.KSort.ChGap];
correctionMatrix=zeros(length(Vector2correct),Exp.KSort.ShankN);

for k=1:Exp.KSort.ShankN
correctedV = smooth(interp1(ys_um(:,k),landMark_um,Vector2correct,'linear','extrap'),100);  %interpolate and smooth by 100units
correctionMatrix(:,k)=correctedV;
end
Vector2correct=-1*Vector2correct(:);
correctionMatrix=-1*correctionMatrix;

end