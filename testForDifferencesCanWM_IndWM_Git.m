clear; clc;
%Name paths

%Origin, DARTEL normalized CoLaus MPM DATA
myDataPath =('//filerch/data/LREN/USERS_DATA/LKhenissi/MPMDartel');

%Destination, image outputs to test for differences between canonical white
%matter mask and individual WM MT
diffCanIndWMTestPath = ('\\filerch\data\LREN\USERS_DATA\LKhenissi\CoLaus_CanMaskTest\WM_Test');

%Now harmonize subject list
cd(myDataPath)
dataFolders = dir;
for i = 1:length(dataFolders)
if strfind(dataFolders(i).name,'PR0')
dataFolderNames{i} = dataFolders(i).name;
end
end
nulls = cellfun(@isempty, dataFolderNames) == 0;

%Remove additional outliers flagged by AL ca 01/19
%morePoorSubjects = {''};
dataFolderNames= dataFolderNames(nulls);
for i = 1:length(morePoorSubjects)
indexPoorSubj(i) =find(strcmpi(morePoorSubjects{i}, dataFolderNames));
end
dataFolderNames = dataFolderNames';
Goods = ones(length(dataFolderNames), 1);
Goods(indexPoorSubj) = 0;
dataFolderNames = dataFolderNames(logical(Goods)); 

% %Ensure subject list matches those in my SES data file 
% SES_Somatic_PR = readtable('//filerch/data/LREN/USERS_DATA/LKhenissi/CuratedCoLausData/SES_Somatic_PRNumber_20190128.csv');
% 
% samePRs =ismember(SES_Somatic_PR.PRNumber, dataFolderNames);
% SES_Somatic_PR =SES_Somatic_PR(samePRs, :);
% clear samePRs
% 
% 
% samePRs =ismember(dataFolderNames, SES_Somatic_PR.PRNumber);
% dataFolderNames =dataFolderNames(samePRs, :);
% clear samePRs
% 
% dataFolderNames = SES_Somatic_PR.PRNumber;

%Make correct list of subject directories in destination folder (do this
%only once, comment/uncomment as needed)

cd(diffCanIndWMTestPath)
for i = 1:length(dataFolderNames)
    mkdir (char(dataFolderNames(i)))
end

%Initialize canonical WM mask 
canonicalWM = 'C:\DATA\software\SPM\spm12\tpm\mask_WM.nii';
%Initialize individual MT WM files
indWM = {};
binIndWM = {};

%Now the main part of the script.
%The aim is to find individual brains that are missing parts relative to
%the canonical mask. 

%1.First, collate individual WM MT
%maps 

%2.Create an individual binary mask of individual WM MT maps, and, in
%addition mask this image by the canonical image. Why? Because we are
%looking for missing brain in individuals relative to the canonical mask
%and not missing brain in the canonical relative to the individual mask
%(yet).

%3.Once an individual binary image of MT WM is obtained, subtract it from
%the canonical mask. Ideal result is an empty image. Expected results are a
%few stray voxels > 0. Unexpected results are a large number of voxels > 0.
%How large is large? To be discussed.

%4.Next, check for these stray voxels in the resulting image.

%5.Finally, sum these voxels along the 3 dimensions to determine if they
%constitute a large number. 

for i = 1:length(dataFolderNames)
    hManToGo = length(dataFolderNames) - i;
    fprintf('%d subjects to test...\n', hManToGo);
    thisSubjectFolder = fullfile(myDataPath, char(dataFolderNames{i}));
    thisSubjectFolderTest = fullfile(diffCanIndWMTestPath, char(dataFolderNames{i}));
    cd(thisSubjectFolder)
    thisSubjectFiles = dir(thisSubjectFolder);
    thisSubjectFileNames={};
    for j = 1:length(thisSubjectFiles)
        if startsWith(thisSubjectFiles(j).name,'swc2s') && endsWith(thisSubjectFiles(j).name,'01_MT.nii')
        indWM{i} = fullfile(thisSubjectFolder,thisSubjectFiles(j).name);
        end
    end
    spm_jobman('initcfg');
    matlabbatch{1}.spm.util.imcalc.input = {
                                           indWM{i}
                                           canonicalWM
                                           };
    matlabbatch{1}.spm.util.imcalc.output = 'binIndWM';
    matlabbatch{1}.spm.util.imcalc.outdir = {thisSubjectFolderTest};
    matlabbatch{1}.spm.util.imcalc.expression = 'i2.*(i1 >0.2)';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    matlabbatch{2}.spm.util.imcalc.input = {
                                        canonicalWM
                                        fullfile(thisSubjectFolderTest,'binIndWM.nii')
                                        };
    matlabbatch{2}.spm.util.imcalc.output = 'diffCanonicalWMbinIndWM';
    matlabbatch{2}.spm.util.imcalc.outdir = {thisSubjectFolderTest};
    matlabbatch{2}.spm.util.imcalc.expression = 'i1 - i2';
    matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{2}.spm.util.imcalc.options.mask = 0;
    matlabbatch{2}.spm.util.imcalc.options.interp = 1;
    matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run', matlabbatch);
    
    thisImg =load_nii(fullfile(diffCanIndWMTestPath, char(dataFolderNames{i}),'diffCanonicalWMbinIndWM.nii'));
    chkWMDiffANY.subject{i} = any(thisImg.img);
    chkWMDiffHOWMANY.subject{i} = sum(chkWMDiffANY.subject{i});
end    

%sum, for each subject, their total missing voxels across dim
for i = 1:length(chkWMDiffHOWMANY.subject)
totalMissingVoxels(i) = sum(chkWMDiffHOWMANY.subject{1,i});
end

totalMissingVoxels= totalMissingVoxels';

%Find the worst case
%Find the worst case
biggestLoss= max(totalMissingVoxels);
rangeLoss= biggestLoss - min(totalMissingVoxels);
modeLoss = mode(totalMissingVoxels);
meanLoss = mean(totalMissingVoxels);

display(biggestLoss)
display(rangeLoss)
display(modeLoss)
display(meanLoss)

pause 

biggestLossi =find(totalMissingVoxels==biggestLoss);

PrimeSuspect = dataFolderNames(biggestLossi);

otherBigLosses = find(totalMissingVoxels>(meanLoss + 2*(std(totalMissingVoxels))));

otherSuspects = dataFolderNames(otherBigLosses);

resFile =fullfile(diffCanIndWMTestPath,'resComparisonCanonicalWMMask_IndSubjWM_20190215.mat');
save(resFile, 'chkWMDiffANY', 'chkWMDiffHOWMANY', 'totalMissingVoxels', 'biggestLoss', 'PrimeSuspect','otherSuspects');
