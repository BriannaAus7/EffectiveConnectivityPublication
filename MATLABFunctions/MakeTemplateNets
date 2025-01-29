
% Make all the template networks needed - 
%{
RedNetNames={"Vis","SoMat","DorsAttn","Sal/VentAttn","Limb","Cntrl","DMN","TempPar","SubCor"}';


%}
addpath('C:\Users\dk818\OneDrive - Imperial College London\NCORE\miFCdFC\AtlasAndOtherInputs')
load templateEigAndMat.mat
clearvars -except VMN DMN TmpPar SubCor SoMat SalVent Limb DorsAttn Cont Vis
names=readtable('ComboNames.xlsx');

% See what networks have VMN ROIs
for ii=1:size(VMN)
names(VMN(ii,1),1)
end

% Now need to get rid of numbers in other networks - 
toDel=ismember(SubCor,VMN);
SubCor(toDel,:)=[];

toDel=ismember(Limb,VMN);
Limb(toDel,:)=[];

toDel=ismember(DMN,VMN);
DMN(toDel,:)=[];

Combo=[VMN;DMN;TmpPar;SubCor;SoMat;SalVent;Limb;DorsAttn;Cont;Vis];

VisCount=0;
VMNCount=0;
DMNCount=0;
TmpParCount=0;
SubCorCount=0;
SoMatCount=0;
SalVentCount=0;
LimbCount=0;
DorsAttnCount=0;
ContCount=0;
VisCount=0;

clearvars ii toDel

save('TemplateNets.mat','-v7.3')



