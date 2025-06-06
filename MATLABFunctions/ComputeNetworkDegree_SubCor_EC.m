% Please download this function and others to ensure that you can run the ECPermutation.m and EdgeCountsAndDigraph.m scripts.

function [WorkingTable,PropRedNetList,RedNetNames] = ComputeNetworkDegree(ResDF)
% No figure made here since it's for EC
% SubCor ROIs added


%%
load TemplateNets.mat

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making a table to store results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All I actually need is ROI1, ROI2
% and then what's computed first is 
% NetworkID 17 Networks ROI, NetworkID 17 Networks ROI, 
% and then the same but for Network ID 7 Networks next
if exist('WorkingTable','var')==1
    clear WorkingTable
end
WorkingTable(:,1:2)=ResDF(:,1:2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Counting ROI1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:size(WorkingTable,1)
   
    if ismember(WorkingTable(ii,1),Vis) == 1
       VisCount=VisCount+1 ;
       WorkingTable(ii,3)=1;
       
   elseif ismember(WorkingTable(ii,1),SoMat) == 1
       SoMatCount=SoMatCount+1 ;
       WorkingTable(ii,3)=2;
        
   elseif ismember(WorkingTable(ii,1),DorsAttn) == 1
       DorsAttnCount=DorsAttnCount+1  ;
       WorkingTable(ii,3)=3;
        
   elseif ismember(WorkingTable(ii,1),SalVent) == 1
       SalVentCount=SalVentCount+1 ;
       WorkingTable(ii,3)=4;
        
   elseif ismember(WorkingTable(ii,1),Limb) == 1
       LimbCount=LimbCount+1 ;
       WorkingTable(ii,3)=5;
             
   elseif ismember(WorkingTable(ii,1),Cont) == 1
       ContCount=ContCount+1 ;    
       WorkingTable(ii,3)=6;
       
   elseif ismember(WorkingTable(ii,1),DMN) == 1
       DMNCount=DMNCount+1  ;
       WorkingTable(ii,3)=7;
        
   elseif ismember(WorkingTable(ii,1),TmpPar) == 1
       TmpParCount=TmpParCount+1 ;    
       WorkingTable(ii,3)=8;  
          
   elseif ismember(WorkingTable(ii,1),SubCor) == 1
       SubCorCount=SubCorCount+1 ;    
       WorkingTable(ii,3)=9;  

  elseif ismember(WorkingTable(ii,1),VMN) == 1
       VMNCount=VMNCount+1 ;    
       WorkingTable(ii,3)=10;  

    else 
        disp(strcat('no match for ',num2str(ii)))
              
   end
       
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Counting ROI2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:size(WorkingTable,1)
   
    if ismember(WorkingTable(ii,2),Vis) == 1
       VisCount=VisCount+1 ;
       WorkingTable(ii,4)=1;
       
   elseif ismember(WorkingTable(ii,2),SoMat) == 1
       SoMatCount=SoMatCount+1 ;
       WorkingTable(ii,4)=2;
        
   elseif ismember(WorkingTable(ii,2),DorsAttn) == 1
       DorsAttnCount=DorsAttnCount+1  ;
       WorkingTable(ii,4)=3;
        
   elseif ismember(WorkingTable(ii,2),SalVent) == 1
       SalVentCount=SalVentCount+1 ;
       WorkingTable(ii,4)=4;
        
   elseif ismember(WorkingTable(ii,2),Limb) == 1
       LimbCount=LimbCount+1 ;
       WorkingTable(ii,4)=5;
             
   elseif ismember(WorkingTable(ii,2),Cont) == 1
       ContCount=ContCount+1 ;    
       WorkingTable(ii,4)=6;
       
   elseif ismember(WorkingTable(ii,2),DMN) == 1
       DMNCount=DMNCount+1  ;
       WorkingTable(ii,4)=7;
        
   elseif ismember(WorkingTable(ii,2),TmpPar) == 1
       TmpParCount=TmpParCount+1 ;    
       WorkingTable(ii,4)=8;  
          
   elseif ismember(WorkingTable(ii,2),SubCor) == 1
       SubCorCount=SubCorCount+1 ;    
       WorkingTable(ii,4)=9;  

   elseif ismember(WorkingTable(ii,2),VMN) == 1
       VMNCount=VMNCount+1 ;    
       WorkingTable(ii,4)=10;  

   else 
        disp(strcat('no match for ',num2str(ii)))
              
   end
       
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute degree per 10 networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RedNetList=[VisCount,SoMatCount,DorsAttnCount,SalVentCount,LimbCount,ContCount,DMNCount,TmpParCount,SubCorCount,VMNCount]';
RedNetNames={"Vis","SoMat","DorsAttn","SalVent","Limb","Cntrl","DMN","TempPar","SubCor","VMN"}';

Redmm=mean(RedNetList)
Redmd=median(RedNetList)
Redsd=std(RedNetList)

PropRedNetList=RedNetList/sum(RedNetList);

end
