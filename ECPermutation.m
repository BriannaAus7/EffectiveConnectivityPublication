% This is the permutation script that seeks to identify which edges had significantly different EC between MD vs HC participants. 
% A minimum of 1000 permutations is run, after which a Kolmogorov-Smirnov test evaluates whether the null distribution was significantly differently from a normal distribution. 
% If it was different, an additional 1000 permutations are run until the null distribution was not significantly different from a normal distribution. 
% Permutation testing was chosen due to its applicability to small samples and previous use in Effective Connectivity studies as described in Goulden et al., 2010; https://www.sciencedirect.com/science/article/pii/S1053811909012452?via%3Dihub

clear all
addpath('/Users/briannaaustin/Desktop/lsngc(2)/EC_Brianna(2)/MIDData/HC_MID')
outputpath='/Users/briannaaustin/Desktop/lsngc(2)/EC_Brianna(2)/MIDData';
addpath(outputpath)


%%

% HC
baseloc='/Users/briannaaustin/Desktop/lsngc(2)/EC_Brianna(2)/MIDData/HC_MID';
% baseloc='C:\Users\dk818\OneDrive - Imperial College London\NCORE\EC_Brianna\Data4LsNGC\lsngc\EC_Brianna(3)\CueData-2\HC_Cue';
addpath(baseloc)
Files=dir(strcat(baseloc,'/*Aff*'));
n_subj=size(Files,1);
num_rois=214;
% What I need to do - get all GC values per ROI pair. 
for subj=1:n_subj
Aff=readtable(Files(subj).name);
% Remove the ROI names
Aff(:,1)=[];
Aff=table2array(Aff);
for rr1=1:num_rois
    for rr2=1:num_rois
        BigAffHC(rr1,rr2,subj)=Aff(rr1,rr2);
    end
end
end


%%
% MD
% baseloc='C:\Users\dk818\OneDrive - Imperial College London\NCORE\EC_Brianna\Data4LsNGC\lsngc\EC_Brianna(3)\CueData-2\Patients_Cue';
baseloc='/Users/briannaaustin/Desktop/lsngc(2)/EC_Brianna(2)/MIDData/Patients_MID';
addpath(baseloc);
Files=dir(strcat(baseloc,'/*Aff*'));
n_subj=size(Files,1);

% What I need to do - get all GC values per ROI pair. 
for subj=1:n_subj
Aff=readtable(Files(subj).name);
Aff(:,1)=[];
Aff=table2array(Aff);
for rr1=1:num_rois
    for rr2=1:num_rois
        BigAffPT(rr1,rr2,subj)=Aff(rr1,rr2);
    end
end
end

%%
% Let's do some tesstttttsss!!! First, let's get Brianna the mean and std
% per group. Then, let's run a sign-rank test for each pair, storing
% outputs in diff tables. Then, let's look at all the networks with
% effects. 
meanHC=nanmean(BigAffHC,3);
meanPT=nanmean(BigAffPT,3);
writematrix(meanHC,strcat(outputpath,'/MeanHCMIDAff.csv'));
writematrix(meanPT,strcat(outputpath,'/MeanMDMIDAff.csv'));


%%
% Permutation testing
% https://www.sciencedirect.com/science/article/pii/S1053811909012452?via%3Dihub#aep-section-id13
% For each region - 
% Get the diff between the mean of the HC and PT groups
% For minimum 1,000 permutations - 
% Compile all groups, and reshuffle the labels
% Take diffs of mean of fake HC and fake PT groups
% Store this value in a vector
% Test distribution of this vector
% If normal, then see what % of values are more extreme than the empirical
% If not normal, add another 1000 permutations, and repeat
h=0;
for ii=1:num_rois
    for jj=1:num_rois

        if ii ~= jj

        HCROI=squeeze(BigAffHC(ii,jj,:));
        PTROI=squeeze(BigAffPT(ii,jj,:));

        % Get the empirical difference - our test stat
        Tstat=nanmean(HCROI)-nanmean(PTROI);
        TstatMatForFig(ii,jj)=Tstat;

        % Combine everyone - we know the real index is that 1:22 are HC,
        % and 23:46 are PT
        Combo=[HCROI;PTROI];
        
        % Start permutations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for perm=1:2000
        % Get the index for the fake labels
        indx=randperm(size(Combo,1));

        % Get the fake HC and fake PT groups
        FakeHC=Combo(indx(1:22),1);
        FakePT=Combo(indx(23:end),1);

        % Get the diff of the mean of both 
        NullDistrib(perm,1)=nanmean(FakeHC)-nanmean(FakePT);
        end

        % Test normality of null distribution
        % If it's not normal, record which roi pair it is, and deal with it later
        if kstest(NullDistrib)==0
        NotNormTest(ii,jj)=1;

        % If it is, then see how many values are more extreme than our test
        % set
        else

        % Significance testing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % If the Tstat is on the neg side, then it's sig lower; if pos,
        % then extreme vals are those that are greater
        if Tstat<0
        p=sum(Tstat>NullDistrib)/size(NullDistrib,1);
        else
        p=sum(Tstat<NullDistrib)/size(NullDistrib,1);
        end

        % Format results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ResTbl_p(ii,jj)=p;

        % Effect size for two independent samples uses pooled SD
        % Link for how to do this:
        % https://www.fieldtriptoolbox.org/example/effectsize/#computing-the-effect-size-by-hand
        % n = number of samples in empirical vs permutation data 
        n1=size(HCROI,1);
        n2=size(PTROI,1);
        pooled_sd = sqrt( ((n1-1)*std(HCROI)^2 + (n2-1)*std(PTROI)^2) / (n1+n2-1) );
        cohensd = abs(nanmean(HCROI)-nanmean(PTROI)) / pooled_sd;
        ResTbl_p(ii,jj)=cohensd;

        % BigResTbl has cols for:
        % ROI1, ROI2, p, FDR p, tstat, effect size
        h=h+1;
        BigRes(h,1)=ii;
        BigRes(h,2)=jj;
        BigRes(h,3)=p;
        BigRes(h,5)=Tstat;
        BigRes(h,6)=cohensd;
        end
        end
    end
end

% FDR correction
[~,~,BigRes(:,4)]=mafdr(BigRes(:,3));

%%
SigEffectsFDR=sum(BigRes(:,4)<0.05);
SigEffects=sum(BigRes(:,3)<0.05);
StatCount=0;
% What prop of all analyis are sig?
for ii=1:num_rois
for jj=1:num_rois
    if ii ~= jj
    StatCount=StatCount+1;
    end
end
end

SigEffects/StatCount
 

save('MID_EC.mat','-v7.3')

% %%
% % Note - BigResTbl is in the correct input for the function
% [WorkingTable,PropRedNetList,RedNetNames]=ComputeNetworkDegree_SubCor_EC(BigRes);

% %%
% BigResTbl=array2table(BigRes);
% BigResTbl.Properties.VariableNames={'r1','r2','p','FDRp','zstat','EffectSize'};

%%
%{
DEPRECATED
Because all the files like cue_sub-NK1004T.csv_Aff.csv, I had to rename them
all since MATLAB wouldnt read that in, so I did the below - 
for subj=1:n_subj
name=extractBefore(Files(subj).name,'.csv');
movefile(Files(subj).name,strcat(name,'_Aff.csv'))
end
%}
