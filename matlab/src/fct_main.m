% FCT pipeline

%%%
function fct_main(inp)

t1_file = inp.t1_niigz;
fmri_file = inp.fmri_niigz;
xnat_project   = inp.xnat_project;
xnat_subject   = inp.xnat_subject;
xnat_session   = inp.xnat_session;
patho = inp.out_dir;

% Convert chars to strings if deployed
%if ischar(t1_file)
%    t1_file=convertCharsToStrings(t1_file);
%    xnat_project=convertCharsToStrings(xnat_project);
%    xnat_subject=convertCharsToStrings(xnat_subject);
%    xnat_session=convertCharsToStrings(xnat_session);
%    patho=convertCharsToStrings(patho);
%end

% root_dir = inptpara.rootDir;
% xnat_project = inptpara.xnat_project;
% pathi = inptpara.path01;
% patho = inptpara.path02;
% 
% if ~isdeployed
%     addpath(genpath([root_dir,'/add2path/spm12']));
%     addpath(genpath([root_dir,'/z_batch_pipeline/add2path/DPABI_V2.3_170105']));
%     addpath(genpath([root_dir '/myfunc/']))
% end
% 
% path = [root_dir '/z_batch_pipeline/'];


% %% check existence of fMRI and T1
% if isempty(dir([pathi '/T1Img/']))
%     disp('no T1 scan'); return;
% end
% if isempty(dir([pathi '/FunImg/']))
%     disp('no fMRI scan'); return;
% end


% %% delete residues
% if exist([patho '/preprocess'], 'dir')
%     rmdir([patho '/preprocess'], 's');
% end

%% Prepare input files and move to output directory
disp('Preparing inputs')

if ~exist(patho,'dir')
    mkdir(patho);
end

%% preprocessing

fprintf('Preprocessing:')
tic

i=1;
mkdir([patho '/preprocess/FunImg/' num2str(i)]);
mkdir([patho '/preprocess/T1Img/' num2str(i)]);

tmp = split(fmri_file,'/');
fmri_name = tmp{end};
tmp = split(t1_file,'/');
t1_name = tmp{end};
%check for json
tmp=dir(fmri_file);
if isempty(dir([tmp.folder,'/*.json']))
    fprintf('No .json files found. Extracting header info from nifti.');
else
    % copy json files
    fmri_js = replace(fmri_file,'.nii.gz','.json'); 
    tmp = split(fmri_js,'/');
    json_name = tmp{end};
    copyfile(fmri_js,sprintf('%s/preprocess/FunImg/%s/%s',patho,num2str(i),json_name));
%            t1_js = [t1_file(1:end-7) '.json'];
%            copyfile(t1_js,sprintf('%s/preprocess/T1Img/%s/%s',patho,num2str(i),[t1_name(1:end-4) '.json']));
end
copyfile(fmri_file,sprintf('%s/preprocess/FunImg/%s/%s',patho,num2str(i),fmri_name)); % original T1 and rsfMRI
copyfile(t1_file,sprintf('%s/preprocess/T1Img/%s/%s',patho,num2str(i) ,t1_name));% original T1 and rsfMRI
gunzip(fullfile([patho '/preprocess/T1Img/' num2str(i)],t1_name));
gunzip(fullfile([patho '/preprocess/FunImg/' num2str(i)],fmri_name));

toc

%% ==== load config file =============================
% Change this path once code is containerized
LOADFILENAME1=which(fullfile('prepro_saved_v5.mat'));

% Load the data file from current working directory
load(which(LOADFILENAME1),'Cfg');


%% ######################## set parameters #######################%
% set dir etc.
fprintf('Set Parameters:')
tic
Cfg.WorkingDir=[patho,'/preprocess'];
Cfg.DataProcessDir=[patho,'/preprocess'];
list=dir([patho,'/preprocess/FunImg/']);
Cfg.SubjectID={};
for i=3:length(list) % multiple scans
    Cfg.SubjectID{1,1}=list(i).name;

    % get num of TRs and num of slices from fMRI header
    tmp = dir([patho,'/preprocess/FunImg/' list(i).name '/*.nii.gz']);
    info = niftiinfo([patho,'/preprocess/FunImg/' list(i).name '/' tmp.name]);

    pdim = info.PixelDimensions;
    if length(pdim)<4
        js = jsondecode(fileread([patho,'/preprocess/FunImg/' list(i).name '/' tmp.name(1:end-7) '.json']));
        TR = js.RepetitionTime;
    else
        TR = pdim(4);
    end
    disp(['TR =' num2str(TR) 's']);

    slnum = info.ImageSize(3); disp(['num of slices =' num2str(slnum)]);
    tpnum = info.ImageSize(4); disp(['num of time points =' num2str(tpnum)]);

    % set default slice order according to the XNAT project
    if strcmp(xnat_project, 'ADNI_23')
        slorder = [1:2:slnum 2:2:slnum];
        csfthresh = 0.9999999;
    elseif strcmp(xnat_project, 'BLSA')
        slorder = 1:slnum;
        csfthresh = 0.995;
    elseif strcmp(xnat_project, 'OASIS-3')
        slorder = 1:slnum;
        csfthresh = 0.99;
    end
    % check whether there is JSON file (BLSA-no; OASIS3-yes; ADNI23-most yes)
    % "SliceTiming" in JSON file indicates slice order
    % if there is no such info, then keep using default slice order
    json_flag = 0;
    if isfile([patho,'/preprocess/FunImg/' list(i).name '/' tmp.name(1:end-7) '.json'])
        js = jsondecode(fileread([patho,'/preprocess/FunImg/' list(i).name '/' tmp.name(1:end-7) '.json']));
        % check whether there is SliceTiming filed in JSON file (OASIS-part yes)
        if isfield(js, 'SliceTiming')
            [~, slorder] = sort(js.SliceTiming); 
            slorder = slorder'; json_flag = 1;
        end
    end

    
    % Set path variable
    tmp = which('Template_0_IXI555_MNI152_GS.nii');
    tmp=split(tmp,'/');
    path = ['/', strjoin(tmp(2:end-1),'/')];
    
    
    % if there is no SliceTiming in JSON, try to find it in summary CSV file
    if strcmp(xnat_project, 'ADNI_23') && json_flag == 0 
        % find slice order in spreadsheet 
        T = readtable([path '/MAYOADIRL_MRI_FMRI_NFQ_05_09_22_modifiedByGao_ADNI23.csv']);
        T_sid = T.RID;   T_vid = T.VISCODE0;  T_vid1 = T.VISCODE00; 
        T_st = T.SLICETIMING;
        sID = str2double(xnat_subject(7:10));    vID = str2double(xnat_subject(17:end));
        ind = find(T_sid==sID & (T_vid==vID | T_vid1==vID), 1);
        if ~isempty(ind)
            st_n = str2double(split(T_st{ind}, '_')); % array of slice timing
            [~, so_i] = sort(st_n); 
            slorder = so_i'; 
        end
    end
    
    clear nii;

    Cfg.TimePoints=tpnum; % set number of time points
    Cfg.TR=TR; % set TR
    % set whether to do slice timing (skip slice timing if TR<1)
    if strcmp(xnat_project, 'ADNI_23') && TR < 1
        Cfg.IsSliceTiming = 0; %  
        disp('No Slice Timing Correction for Multi-Band');
    else
        Cfg.IsSliceTiming = 1; % 
        disp(['slice order= ' num2str(slorder)]);
    end
    Cfg.SliceTiming.SliceNumber=slnum; % set number of slices
    Cfg.SliceTiming.SliceOrder=slorder; % set slice order
    Cfg.SliceTiming.ReferenceSlice=round(slnum/2); % set reference slice
    Cfg.Covremove.CSF.MaskThreshold=csfthresh;

    save([patho '/prepro_saved_v5_' list(i).name '.mat'], 'Cfg', '-v7.3'); % modified in v5
    DPARSFA_run_XNAT_v5([patho '/prepro_saved_v5_' list(i).name '.mat']); % modified in v5

    toc
    %% fMRI->T1, T1 segmentation and normalization
    fprintf('T1 segmentation:')
    tic
    
    mkdir([patho,'/preprocess/CatNormalization/' list(i).name]);

% Copy files to CatNormalization folder
    tmp = dir([patho,'/preprocess/T1Img/' list(i).name '/*.nii']);
    copyfile( [patho,'/preprocess/T1Img/' list(i).name '/' tmp.name], ...
        [patho,'/preprocess/CatNormalization/' list(i).name '/T1.nii']);

    if exist([patho,'/preprocess/FunImgRCD/'], 'dir')
        copyfile( [patho,'/preprocess/FunImgRCD/' list(i).name '/Detrend_4DVolume.nii'],...
            [patho,'/preprocess/CatNormalization/' list(i).name '/Detrend_4DVolume.nii']);
        copyfile( [patho,'/preprocess/FunImgRCDF/' list(i).name '/Filtered_4DVolume.nii'],...
            [patho,'/preprocess/CatNormalization/' list(i).name '/Filtered_4DVolume.nii']);
    else
        copyfile( [patho,'/preprocess/FunImgARCD/' list(i).name '/Detrend_4DVolume.nii'],...
            [patho,'/preprocess/CatNormalization/' list(i).name '/Detrend_4DVolume.nii']);
        copyfile( [patho,'/preprocess/FunImgARCDF/' list(i).name '/Filtered_4DVolume.nii'],...
            [patho,'/preprocess/CatNormalization/' list(i).name '/Filtered_4DVolume.nii']);
    end
    % Segmentation based on T1, output = tissue-masks
    LOADFILENAME1=which(fullfile('saved_cat12.mat'));

    % Load the data file from current working directory
    load(which(LOADFILENAME1),'matlabbatch');
    
    matlabbatch{1,1}.spm.tools.cat.estwrite.data{1,1}=[patho,'/preprocess/CatNormalization/' list(i).name '/T1.nii,1'];
    matlabbatch{1,1}.spm.tools.cat.estwrite.opts.tpm{1,1}=[path,'/TPM.nii'];
    matlabbatch{1,1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm{1,1}=[path,'/Template_0_IXI555_MNI152_GS.nii'];
    spm_jobman('run',matlabbatch);

    % Coregister fMRI to T1
    
    %load([path,'/saved_coreg.mat'],'matlabbatch');
    
    LOADFILENAME1=which(fullfile('saved_coreg.mat'));

    % Load the data file from current working directory
    load(which(LOADFILENAME1),'matlabbatch');
    
    matlabbatch{1,1}.spm.spatial.coreg.estimate.ref{1,1}=[patho,'/preprocess/CatNormalization/' list(i).name '/T1.nii,1'];
    lstmp = dir([patho,'/preprocess/RealignParameter/' list(i).name '/mean*']);
    matlabbatch{1,1}.spm.spatial.coreg.estimate.source{1,1}=[patho,'/preprocess/RealignParameter/' list(i).name '/' lstmp.name];
    for nn=1:tpnum
        matlabbatch{1,1}.spm.spatial.coreg.estimate.other{nn,      1}=[patho,'/preprocess/CatNormalization/' list(i).name '/Detrend_4DVolume.nii,',num2str(nn)];
    end
    for nn=1:tpnum
        matlabbatch{1,1}.spm.spatial.coreg.estimate.other{nn+tpnum,1}=[patho,'/preprocess/CatNormalization/' list(i).name '/Filtered_4DVolume.nii,',num2str(nn)];
    end
    spm_jobman('run',matlabbatch);

    % Normalize T1/tissue-masks/ALFF/fALFF/Reho/fMRI (->MNI space)
    
    LOADFILENAME1=which(fullfile('saved_warp.mat'));

    % Load the data file from current working directory
    load(which(LOADFILENAME1),'matlabbatch');
    %load([path,'/saved_warp.mat'],'matlabbatch');
     matlabbatch{1,1}.spm.spatial.normalise.write.subj.def{1,1}=[patho,'/preprocess/CatNormalization/' list(i).name '/mri/y_T1.nii'];
    matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{1,1} =[patho,'/preprocess/CatNormalization/' list(i).name '/mri/p0T1.nii,1'];
    matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{2,1} =[patho,'/preprocess/CatNormalization/' list(i).name '/mri/p1T1.nii,1'];
    matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{3,1} =[patho,'/preprocess/CatNormalization/' list(i).name '/mri/p2T1.nii,1'];
    matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{4,1} =[patho,'/preprocess/CatNormalization/' list(i).name '/mri/p3T1.nii,1'];
    matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{5,1} =[patho,'/preprocess/CatNormalization/' list(i).name '/T1.nii,1'];

    if exist([patho,'/preprocess/FunImgRCD/'], 'dir')
        matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{6,1} =[patho,'/preprocess/Results/ALFF_FunImgRCD/ALFFMap_' list(i).name '.nii,1'];
        matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{7,1} =[patho,'/preprocess/Results/DegreeCentrality_FunImgRCDF/DegreeCentrality_PositiveBinarizedSumBrainMap_' list(i).name '.nii,1'];
        matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{8,1} =[patho,'/preprocess/Results/DegreeCentrality_FunImgRCDF/DegreeCentrality_PositiveWeightedSumBrainMap_' list(i).name '.nii,1'];
        matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{9,1} =[patho,'/preprocess/Results/fALFF_FunImgRCD/fALFFMap_' list(i).name '.nii,1'];
        matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{10,1}=[patho,'/preprocess/Results/ReHo_FunImgRCDF/ReHoMap_' list(i).name '.nii,1'];
    else
        matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{6,1} =[patho,'/preprocess/Results/ALFF_FunImgARCD/ALFFMap_' list(i).name '.nii,1'];
        matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{7,1} =[patho,'/preprocess/Results/DegreeCentrality_FunImgARCDF/DegreeCentrality_PositiveBinarizedSumBrainMap_' list(i).name '.nii,1'];
        matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{8,1} =[patho,'/preprocess/Results/DegreeCentrality_FunImgARCDF/DegreeCentrality_PositiveWeightedSumBrainMap_' list(i).name '.nii,1'];
        matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{9,1} =[patho,'/preprocess/Results/fALFF_FunImgARCD/fALFFMap_' list(i).name '.nii,1'];
        matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{10,1}=[patho,'/preprocess/Results/ReHo_FunImgARCDF/ReHoMap_' list(i).name '.nii,1'];
    end
    for nn=1:tpnum
        disp(nn);
        matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{nn+10,      1}=[patho,'/preprocess/CatNormalization/' list(i).name '/Detrend_4DVolume.nii,',num2str(nn)];
    end
    for nn=1:tpnum
        disp(nn)
        matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample{nn+10+tpnum,1}=[patho,'/preprocess/CatNormalization/' list(i).name '/Filtered_4DVolume.nii,',num2str(nn)];
    end

    spm_jobman('run',matlabbatch);

    delete([patho,'/preprocess/CatNormalization/' list(i).name '/T1.nii']);
    delete([patho,'/preprocess/CatNormalization/' list(i).name '/Detrend_4DVolume.nii']);
    delete([patho,'/preprocess/CatNormalization/' list(i).name '/Filtered_4DVolume.nii']);
    
    toc
end

fprintf('fmriprocessing:')
tic

% saveDir = [path01 '/FCT/'];
%%%brain mask
if ~exist([patho '/results'],'dir')
    mkdir([patho '/results']);
end
list1 = dir([patho '/preprocess/CatNormalization']);
list1(1:2,:) = [];
for ii = 1 : length(list1)
    if ~exist([patho '/results/' list1(ii).name],'dir')
        mkdir([patho '/results/' list1(ii).name]);
    end
    copyfile([patho '/preprocess/CatNormalization/' list1(ii).name '/mni_Detrend_4DVolume.nii']...
        ,[patho '/results/' list1(ii).name,'/']);
end





%%%%%%%%%%%%%

flags = struct('mask', false, 'mean', false, 'interp', 1, 'which', 1, 'wrap', [0 0 0], 'prefix', 'r');
path_brod = which('rBrodmann_YCG.nii');
% path_mask = [patho '/preprocess/Masks/AllResampled_BrainMask_05_91x109x91.nii'];


cd(patho);
for i = 1 : length(list1)

    path_mask = [patho '/preprocess/CatNormalization/' list1(i).name '/mri/mni_p0T1.nii'];
    spm_reslice_quiet({path_brod path_mask},flags);
    path1 = [[patho '/results/'] list1(i).name '/mni_Detrend_4DVolume.nii'];
    spm_reslice_quiet({path_brod path1},flags);
    copyfile([patho '/preprocess/CatNormalization/' list1(i).name '/mri/rmni_p0T1.nii'],...
    [patho '/brainmask.nii']);
    %%%% save dir
%     mkdir([patho '/results/',list1(i).name,'/']);
    saveDir = [patho '/results/',list1(i).name,'/'];
    %%% set up parameters for FTI functions
    %     info=niftiinfo([outDir list1(i).name '/rDetrend_4DVolume.nii']);
    nhood = 7;
    rpower = 2;


%%% Step 0: preprocess the fmri
    %     GMnii  = load_nii(strcat(path01,'/brainmask.nii'));
    GM  = spm_vol(strcat(patho,'/brainmask.nii')); GM = spm_read_vols(GM);
    GM = ones(size(GM)).*(GM>0);
    GM = flip(GM,1);
    %     fMRnii = load_nii(strcat(rootDir,list1(i).name,'/rest1.nii'));
    fMR = spm_vol(strcat([patho '/results/'],list1(i).name,'/rmni_Detrend_4DVolume.nii')); fMR = spm_read_vols(fMR);
    fMR = flip(fMR,1);
    xcell_mask = GM;
    xcell = fMR;
    xcell = double(xcell);

    % space smooth
    FWHM = 14;
    voxelsize = 2;

    for time = 1 : size(xcell,4)
        xcell(:,:,:,time) = imgaussfilt3(xcell(:,:,:,time),FWHM/voxelsize/2.355);
        %     xcell(:,:,:,time) = smooth3(xcell(:,:,:,time),'gaussian',[3 15 3],0.65);
    end

    % generate signal matrix
    [xcellmatrix, I, J, K,~,size01] = mypreprocessxcell(xcell, xcell_mask);
    xcellmatrix1 = xcellmatrix;
    % regress motion and physiological parameters

    % remove time points
    N_delete = 0;
    xcellmatrix1(1:N_delete,:) = [];
    size01_new = [size01(1),size01(2),size01(3),size01(4)-N_delete];

    Fs = 1/TR;
    xcellmatrix1 = xcellmatrix1-mean(xcellmatrix1);
    xcellmatrix1 = detrend(xcellmatrix1);
    xcellmatrix2 = mylpfilt(xcellmatrix1,TR);

    % normalize and back to 4D
    fMR1 = back2fourd(xcellmatrix2,I,J,K,size01_new);
    y0 = 1;
    %%% calculate FCT 5D
    fti5d= reconFTI5d_m4_sm_saveFTI5d(GM, double(fMR1), nhood,rpower,y0);
    save(strcat(saveDir,'fti5d.mat'),'fti5d');
    %%% save FTI tensor figs for following QA steps
    %fMRmean = mean(fMR1,4,'omitnan');
    %FA1 = ones(size(fMR1,1),size(fMR1,2),size(fMR1,3));
    %f1 = plottensortotal1(fMRmean,fti5d,GM,1,44,FA1); % axial slice
    %f2 = plottensortotal1(fMRmean,fti5d,GM,2,56,FA1); % coronal slice
    %f3 = plottensortotal1(fMRmean,fti5d,GM,3,46,FA1); % sagital slice
    %savefig(f1,strcat(saveDir,'axialslice.fig'));
    %savefig(f2,strcat(saveDir,'coronalslice.fig'));
    %savefig(f3,strcat(saveDir,'sagitalslice.fig'));
    %close all;
end
toc
%% remove preprocess to save some space in ACCRE
if exist([patho '/preprocess/FunImg/'], 'dir'), rmdir([patho '/preprocess/FunImg/'], 's'); end
if exist([patho '/preprocess/FunImgA/'], 'dir'), rmdir([patho '/preprocess/FunImgA/'], 's'); end
if exist([patho '/preprocess/FunImgAR/'], 'dir'), rmdir([patho '/preprocess/FunImgAR/'], 's'); end
if exist([patho '/preprocess/FunImgARC/'], 'dir'), rmdir([patho '/preprocess/FunImgARC/'], 's'); end
% if exist([patho '/preprocess/Masks/'], 'dir'), rmdir([patho '/preprocess/Masks/'], 's'); end
if exist([patho '/preprocess/T1Img/'], 'dir'), rmdir([patho '/preprocess/T1Img/'], 's'); end
if exist([patho '/preprocess/T1ImgCoreg/'], 'dir'), rmdir([patho '/preprocess/T1ImgCoreg/'], 's'); end
if exist([patho '/preprocess/T1ImgNewSegment/'], 'dir'), rmdir([patho '/preprocess/T1ImgNewSegment/'], 's'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Zip all nifti files

cd(patho)

A = dir(patho);
for i = 3:length(A)
    if A(i).isdir == 0
        continue
    else 
        cd(A(i).folder);
        B = dir(A(i).name);
        for j = 3:length(B)
            if endsWith(B(j).name,'.nii')
                cd(B(j).folder);
                gzip(B(j).name);
                delete(B(j).name);
            elseif B(j).isdir == 0
                continue
            else
                cd(B(j).folder);
                C = dir(B(j).name);
                for k = 3:length(C)
                    if C(k).isdir == 1
                        cd(C(k).folder);
                        D = dir(C(k).name);
                        for l=3:length(D)
                            if D(l).isdir == 1
                                cd(D(l).folder);
                                E = dir(D(l).name);
                                for m=3:length(E)
                                    if endsWith(E(m).name,'.nii')
                                        cd(E(m).folder);
                                        gzip(E(m).name);
                                        delete(E(m).name);
                                    else
                                        continue
                                    end
                                end
                            elseif endsWith(D(l).name,'.nii')
                                cd(D(l).folder);
                                gzip(D(l).name);
                                delete(D(l).name);
                            else
                                continue
                            end
                        end
                    elseif endsWith(C(k).name,'.nii')
                        cd(C(k).folder);
                        gzip(C(k).name);
                        delete(C(k).name);
                    else
                        continue
                    end
                end
            end
        end
    end
end

disp('all done');
exit

end

%%%%%%END!!

