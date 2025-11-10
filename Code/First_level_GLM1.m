clear
clc
spm fmri

FolderPath = 'E:\CBCS_fMRI\GIIN\data_fMRI\Nifti_New'; %Change it to your nifti path
OutputPath = 'E:\CBCS_fMRI\GIIN\data_fMRI\First_level\First_level_0'; %Change it to your first level path

cd('E:\CBCS_fMRI\GIIN\data_fMRI\First_level'); %Change it to your first level path

Onset_0  %onset file name

NSlice = 62; %% Number of slices
refSlice = 31; %% Reference slice
TR = 2.0;

SubList = dir(fullfile(FolderPath,'GIIN*')); %text that can recognize folders for all participants

% Regressor list
CondList = {'A1S3' 'A1S1' 'A3S1' 'Cost' 'Allocation' 'Choice_Miss' 'Allocation_Miss'};  %% conresponding to onset condition

% parametric regressor name
Parametric_under_condition = { };

% contrast name
cname = { 'A1S3', 'A1S1', 'A3S1'};
% Usually I prefer to set 4 simple effects to extract beta values and for
% future second-level flexible analysis

ctype = { 'T', 'T', 'T'}; % T test for usual contrasts, and F contrast for complex contrasts and PPI/DCM analysis

% The number of columes is the same as the number of regressors
% The number of rows is the same as the number of contrasts
simple_cons = [ 1  0  0  0  0  0  0;   
                0  1  0  0  0  0  0;   
                0  0  1  0  0  0  0;      
                 ];

para_cons = [ ];

for nsub = 1:length(SubList)
    
    s = regexp(SubList(nsub).name,'_','split');
    subnum = s(2);
    subnum = subnum{1};
    
    %     substr =  regexp(s(1),'I','split');
    %     subnum = substr{1,1}{2};
    
    SubList(nsub).name
    SubOutput = fullfile(OutputPath,SubList(nsub).name);
    if ~exist(SubOutput,'dir')
        mkdir(SubOutput);
    end
    
    
    for c= 1:length(cname)
        contrast_dir{c} = fullfile(SubOutput,['contrast',num2str(c),'_',char(cname{c})]);
        if ~exist(contrast_dir{c})
            mkdir(contrast_dir{c});
        end
    end
    
    Checkfolder = fullfile(SubOutput,'ChkMsg2EstimateDone');
    
    RunList = dir(fullfile(FolderPath,SubList(nsub).name,'*run*')); %text that can recognize folders for all runs
    
    
    clear matlabbatch;
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = {SubOutput};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';  %Or you can set it as 'secs'. Be careful your onsets and durations should be constant. 
    
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = NSlice;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = refSlice;
    
    for nrun = 1:length(RunList)
        
        clear Files FileList FilePath;
        
        FilePath = fullfile(FolderPath,SubList(nsub).name,RunList(nrun).name);
        FileList = dir(fullfile(FilePath,'swra*.nii'));
        
        %%
        for nfile = 1:length(FileList)
            Files(nfile) = {fullfile(FilePath,FileList(nfile).name)};
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).scans = cellstr(Files)';
        
        
        runnum = nrun;
        cond_num = length(CondList);
        miss_plus = [999];

% Regressor list
% CondList = {'A1S3' 'A1S1' 'A3S1' 'Cost' 'Allocation' 'Choice_Miss' 'Allocation_Miss'};  %% conresponding to onset condition

        for ncond = 1:cond_num
            matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).cond(ncond).name = CondList{ncond};
            
            %%% Onset %%%
            if ncond == 6|| ncond == 7
                matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).cond(ncond).onset = [eval([sprintf('sub%s',subnum) '_' sprintf('run%d',runnum) '_' CondList{ncond} '_onset'])/2  miss_plus];
            else
                matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).cond(ncond).onset = eval([sprintf('sub%s',subnum) '_' sprintf('run%d',runnum) '_' CondList{ncond} '_onset'])/2;
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).cond(ncond).tmod = 0;
            
            %%% Duration %%%
            if ncond == 6 || ncond == 7
                matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).cond(ncond).duration = 10.5/2;
            else
                matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).cond(ncond).duration = eval([sprintf('sub%s',subnum) '_' sprintf('run%d',runnum) '_' CondList{ncond} '_duration'])/2;
            end

            %%% Parametric %%%
            if isempty(Parametric_under_condition)
                matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).cond(ncond).pmod = struct('name', {}, 'param', {}, 'poly', {});
            else
                
                if ncond == 99  %if you have parametirc regressor, indicate the condition to add it; if not, change it to a very large number
                    for p_c = 1:length(Parametric_under_condition)
                        p_c
                        
                        matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).cond(ncond).pmod(p_c).name = char(Parametric_under_condition{p_c});
                        matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).cond(ncond).pmod(p_c).param = beta_sac*eval([sprintf('sub%s',subnum) '_' sprintf('run%d',runnum) '_sac_nosac0'])+beta_bene*eval([sprintf('sub%s',subnum) '_' sprintf('run%d',runnum) '_bene_nosac0']);
                        matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).cond(ncond).pmod(p_c).poly = 1;
                    end
                    matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).cond(ncond).orth = orthogonalize;
                else
                    matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).cond(ncond).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).cond(ncond).orth = orthogonalize;
                end
            end
        end
        
        matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).regress = struct('name', {}, 'val', {});
        MultiReg = dir(fullfile(FilePath,'rp_a*.txt'));
        matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).multi_reg = {fullfile(FilePath,MultiReg.name)};
        matlabbatch{1}.spm.stats.fmri_spec.sess(nrun).hpf = 128;
        
        
    end
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; %  [1 0]  time derives
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.2;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    if ~exist(Checkfolder,'dir')
        spm_jobman('run',matlabbatch);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save design matrix %%%%%%%%%%%%%%%%%%%%%
    clear matlabbatch;
    cd(SubOutput);
    matlabbatch{1}.spm.util.print.fname = 'DsgnMatrix';
    matlabbatch{1}.spm.util.print.opts.opt = {'-dpsc2'; '-append'}';
    matlabbatch{1}.spm.util.print.opts.append = true;
    matlabbatch{1}.spm.util.print.opts.ext = '.ps';
    spm_jobman('run',matlabbatch);
    save(fullfile(SubOutput,sprintf('CntrP.mat')));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% estimate model %%%%%%%%%%%%%%%%%%%%%
    
    clear matlabbatch;
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(SubOutput,'SPM.mat')};
    spm_jobman('run',matlabbatch);
    mkdir(Checkfolder);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% define contrasts%%%%%%%%%%%%%%%%%%%%
    Checkfolder = fullfile(SubOutput,'ChkMsg2ContrastDone');
    if exist(Checkfolder,'dir')
    else
        SPMest=load(fullfile(SubOutput,'SPM.mat'));
        SPMest=SPMest.SPM;
        SPMest.xCon = [];
        headmotion_constant = [0 0 0 0 0 0];
        for c= 1:length(cname)
            cons = [];
            
            for r = 1:length(RunList)
                
                combined_cons = [];
                run_cons = [];
                
                r_in = r;
                
                for s_c = 1:length(simple_cons(c,:))
                    if s_c == 99  %if you have parametirc regressor, indicate the condition to add it; if not, change it to a very large number
                        combined_cons = [combined_cons simple_cons(c,s_c)  para_cons(c,:)  ];
                    else
                        combined_cons = [combined_cons simple_cons(c,s_c)  ];
                    end
                end
                
                run_cons = [combined_cons  headmotion_constant];
                
                cons = [cons run_cons];
            end
            
            cons = [cons zeros(1,length(RunList))];
            
            contrast(c).cname = char(cname(c));
            contrast(c).ctype = char(ctype(c));
            
            contrast(c).cons = cons';
            
            if isempty(SPMest.xCon)
                SPMest.xCon = spm_FcUtil('Set',contrast(c).cname, contrast(c).ctype,'c',contrast(c).cons,SPMest.xX.xKXs);
            else
                SPMest.xCon (end+1) = spm_FcUtil('Set',contrast(c).cname, contrast(c).ctype,'c',contrast(c).cons,SPMest.xX.xKXs);
            end
            
            dlmwrite(fullfile(SubOutput,[SubList(nsub).name,'_contrast',num2str(c),'.txt']), cons','delimiter','\t');
        end
        spm_contrasts(SPMest);
        save(fullfile(SubOutput,[SubList(nsub).name,'_','1stLevel_contrast.mat']),'contrast');
        
        
        %% copy contrast files
        for c= 1:length(cname)
            sourcefile = ['con_',strrep(num2str(c+100000000),'10000','')];
            copyfile(fullfile(SubOutput,[sourcefile,'.nii']),fullfile(contrast_dir{c},[SubList(nsub).name,'_',sourcefile,'.nii']));
        end
        
        mkdir(Checkfolder);
    end
end