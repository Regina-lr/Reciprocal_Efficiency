clear;
clc;
spm fmri

func_path = 'E:\CBCS_fMRI\GIIN\data_fMRI\Nifti_New';

folderstruct = dir(fullfile(func_path,'GIIN*'));
FunImg_sublist = {};

NSlice = 62;
TR = 2;
dcm_file = 'E:\CBCS_fMRI\GIIN\data_fMRI\Raw\GIIN_301_LuoAoxuan\1005-sms_bold_2mm_run1\GIIN_301_LuoAoxuan-1005-sms_bold_2mm_run1-00001.dcm';  %%%One priginal image of one participant
Slice0 = dicominfo(dcm_file);
SliceOrder = Slice0.Private_0019_1029;
s0 = sort(SliceOrder);
refSlice = s0(31);
VoxSize = [3 3 3];
smooth_fwh = [8 8 8];
start_image = 1;  %%whether del first 5 volumes; if delete, 6
spm_path = 'F:\MATLAB_Toolbox\spm12';

    
for i = 1:length(folderstruct)
    subj{i} = folderstruct(i).name;
end;
%%
for Si = 1:length(subj)
%for Ti = 1:length(sess)
%   
% if Si ~= 8 
%     continue
% end


    clear matlabbatch;
% %     
    data_f1 = [];
    data_a1 = [];
    data_r1 = [];
    data_w1 = [];
    
    data_f2 = [];
    data_a2 = [];
    data_r2 = [];
    data_w2 = [];
    
    data_f3 = [];
    data_a3 = [];
    data_r3 = [];
    data_w3 = [];

   
     
    % %     
    run_path = fullfile(func_path,subj{Si});
    
    
    run_name = dir(fullfile(func_path,subj{Si},'*run1*'));
    run1_files = dir(fullfile(run_path,run_name.name,'GIIN*.nii'));
    for i = 1:length(run1_files)
        data_f1 = [data_f1;fullfile(run_path,run_name.name),'/',[run1_files(i).name],',1'];
        data_a1 = [data_a1;fullfile(run_path,run_name.name),'/a',[run1_files(i).name],',1'];
        data_r1 = [data_r1;fullfile(run_path,run_name.name),'/ra',[run1_files(i).name],',1'];
        data_w1 = [data_w1;fullfile(run_path,run_name.name),'/wra',[run1_files(i).name],',1'];
    end;
    
%         data_f1 (1:5,:) = [];
%         data_a1 (1:5,:) = [];
%         data_r1 (1:5,:) = [];
%         data_w1 (1:5,:) = [];

    
    run_name = dir(fullfile(func_path,subj{Si},'*run2*'));
    run2_files = dir(fullfile(run_path,run_name.name,'GIIN*.nii'));
    for i = 1:length(run2_files)
        data_f2 = [data_f2;fullfile(run_path,run_name.name),'/',[run2_files(i).name],',1'];
        data_a2 = [data_a2;fullfile(run_path,run_name.name),'/a',[run2_files(i).name],',1'];
        data_r2 = [data_r2;fullfile(run_path,run_name.name),'/ra',[run2_files(i).name],',1'];
        data_w2 = [data_w2;fullfile(run_path,run_name.name),'/wra',[run2_files(i).name],',1'];
    end;
    
%         data_f2 (1:5,:) = [];
%         data_a2 (1:5,:) = [];
%         data_r2 (1:5,:) = [];
%         data_w2 (1:5,:) = [];
%     
    run_name = dir(fullfile(func_path,subj{Si},'*run3*'));
    run3_files = dir(fullfile(run_path,run_name.name,'GIIN*.nii'));
    for i = 1:length(run3_files)
        data_f3 = [data_f3;fullfile(run_path,run_name.name),'/',[run3_files(i).name],',1'];
        data_a3 = [data_a3;fullfile(run_path,run_name.name),'/a',[run3_files(i).name],',1'];
        data_r3 = [data_r3;fullfile(run_path,run_name.name),'/ra',[run3_files(i).name],',1'];
        data_w3 = [data_w3;fullfile(run_path,run_name.name),'/wra',[run3_files(i).name],',1'];
    end;
    
%         data_f3 (1:5,:) = [];
%         data_a3 (1:5,:) = [];
%         data_r3 (1:5,:) = [];
%         data_w3 (1:5,:) = [];

     
        
        data_a = [];
        data_r = [];
        data_w = [];
%     for runi = 1:10
      for runi = 1:3
        data_a = [data_a;cellstr(eval(strcat('data_a',num2str(runi))))];
        data_r = [data_r;cellstr(eval(strcat('data_r',num2str(runi))))];
        data_w = [data_w;cellstr(eval(strcat('data_w',num2str(runi))))];
      end
%     
% %     
    %_______SLICE TIMING_______
    
%%%%%%    spm8   %%%%%%%
   
%     zatlabbatch{1,1}.spm.temporal.st.scans{1,1} = cellstr(data_f1);
%     zatlabbatch{1,1}.spm.temporal.st.scans{1,2} = cellstr(data_f2);
%     zatlabbatch{1,1}.spm.temporal.st.scans{1,3} = cellstr(data_f3);


%     zatlabbatch{1,1}.spm.temporal.st.nslices = 35;
%     zatlabbatch{1,1}.spm.temporal.st.tr = 2;
%     zatlabbatch{1,1}.spm.temporal.st.ta = 2-2/35;
%     zatlabbatch{1,1}.spm.temporal.st.so = [1:2:35 2:2:35];
%     zatlabbatch{1,1}.spm.temporal.st.refslice = 1;
%     zatlabbatch{1,1}.spm.temporal.st.prefix = 'a';


%%%%%%    spm12   %%%%%%%
    

    matlabbatch{1}.spm.temporal.st.scans = {cellstr(data_f1) cellstr(data_f2) cellstr(data_f3)}';
    
    matlabbatch{1}.spm.temporal.st.nslices = NSlice;
    matlabbatch{1}.spm.temporal.st.tr = TR;
    matlabbatch{1}.spm.temporal.st.ta = TR*(NSlice-1)/NSlice;
    matlabbatch{1}.spm.temporal.st.so = SliceOrder;
    matlabbatch{1}.spm.temporal.st.refslice = refSlice;
    matlabbatch{1}.spm.temporal.st.prefix = 'a';

    spm_jobman('run',matlabbatch)
    %_______REALIGN: Estimation and Reslice_______
    
%%%%%   spm8 %%%%%%%%    
%     zatlabbatch{1,2}.spm.spatial.realign.estimate.data{1,1} = cellstr(data_a1);
%     zatlabbatch{1,2}.spm.spatial.realign.estimate.data{1,2} = cellstr(data_a2);
%     zatlabbatch{1,2}.spm.spatial.realign.estimate.data{1,3} = cellstr(data_a3);
%     zatlabbatch{1,2}.spm.spatial.realign.estimate.data{1,4} = cellstr(data_a4);


%     zatlabbatch{1,2}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
% 	zatlabbatch{1,2}.spm.spatial.realign.estimate.eoptions.sep = 4;
% 	zatlabbatch{1,2}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
% 	zatlabbatch{1,2}.spm.spatial.realign.estimate.eoptions.rtm = 0; 
% 	zatlabbatch{1,2}.spm.spatial.realign.estimate.eoptions.interp = 2; 
% 	zatlabbatch{1,2}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0]; 
% 	zatlabbatch{1,2}.spm.spatial.realign.estimate.eoptions.weight = {};
% 	
% %     zatlabbatch{1,3}.spm.spatial.realign.write.data = cellstr(data_a);
%     zatlabbatch{1,3}.spm.spatial.realign.write.data = data_a;
% 	zatlabbatch{1,3}.spm.spatial.realign.write.roptions.which = [2 1]; 
% 	zatlabbatch{1,3}.spm.spatial.realign.write.roptions.interp = 4; 
% 	zatlabbatch{1,3}.spm.spatial.realign.write.roptions.wrap = [0 0 0]; 
% 	zatlabbatch{1,3}.spm.spatial.realign.write.roptions.mask = 1; 
%     zatlabbatch{1,3}.spm.spatial.realign.write.roptions.prefix = 'r';
%     
%%%%%   spm12 %%%%%%%% 
    clear matlabbatch;   

    matlabbatch{1}.spm.spatial.realign.estwrite.data = {cellstr(data_a1) cellstr(data_a2) cellstr(data_a3)}';
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1; % realign to mean
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';    
    
    spm_jobman('run',matlabbatch)
    
    
%     	%________________________SEGMENTATION (deriving normalisation parameters) _______________________
% run_name = dir(fullfile(func_path,subj{Si},'*para1*'));
% run1_path = fullfile(run_path,run_name.name)
% 
% 
% % 	mean_realign=[run1_path,'/meana',[run1_files(6).name],',1'];
%     mean_realign=[run1_path,'/meana',[run1_files(1).name],',1'];%Siemens
% 
%     zatlabbatch{1,4}.spm.spatial.preproc.data=cellstr(mean_realign);
% 	zatlabbatch{1,4}.spm.spatial.preproc.output.GM=[0 0 1];
% 	zatlabbatch{1,4}.spm.spatial.preproc.output.WM=[0 0 1];
% 	zatlabbatch{1,4}.spm.spatial.preproc.output.CSF=[0 0 0];
% 	zatlabbatch{1,4}.spm.spatial.preproc.output.biascor=1;
% 	zatlabbatch{1,4}.spm.spatial.preproc.output.cleanup=0;
% 	zatlabbatch{1,4}.spm.spatial.preproc.opts.tpm={['/home/Toolbox/spm8/tpm/grey.nii'];['/home/Toolbox/spm8/tpm/white.nii'];['/home/Toolbox/spm8/tpm/csf.nii']};
% 	
% 		
% 	zatlabbatch{1,4}.spm.spatial.preproc.opts.ngaus= [2 2 2 4];
% 	zatlabbatch{1,4}.spm.spatial.preproc.opts.regtype='mni';
% 	zatlabbatch{1,4}.spm.spatial.preproc.opts.warpreg=1;
% 	zatlabbatch{1,4}.spm.spatial.preproc.opts.warpco=25;
% 	zatlabbatch{1,4}.spm.spatial.preproc.opts.biasreg=0.0001;
% 	zatlabbatch{1,4}.spm.spatial.preproc.opts.biasfwhm=60;
% 	zatlabbatch{1,4}.spm.spatial.preproc.opts.samp=3;
% 	zatlabbatch{1,4}.spm.spatial.preproc.opts.msk={''};
    
    
    %_____________________NORMALIZE: Estimate and Write_____________________
    
    %%%%%   spm8 %%%%%%%% 
    %sour_file = [run1_path,'/meana',[run1_files(6).name],',1'];%GE
    
%     sour_file = [run1_path,'/meana',[run1_files(1).name],',1'];%Siemens
% 
%     zatlabbatch{1,5}.spm.spatial.normalise.est.subj.source = cellstr(sour_file);
%     zatlabbatch{1,5}.spm.spatial.normalise.est.subj.wtsrc = {};
%     
% 	zatlabbatch{1,5}.spm.spatial.normalise.est.eoptions.template{1,1} = '/home/Toolbox/spm8/templates/EPI.nii';
%     zatlabbatch{1,5}.spm.spatial.normalise.est.eoptions.weight = {};
%     zatlabbatch{1,5}.spm.spatial.normalise.est.eoptions.smosrc = 8;
%     zatlabbatch{1,5}.spm.spatial.normalise.est.eoptions.smoref = 0;
%     zatlabbatch{1,5}.spm.spatial.normalise.est.eoptions.regtype = 'mni';
%     zatlabbatch{1,5}.spm.spatial.normalise.est.eoptions.cutoff = 25;
%     zatlabbatch{1,5}.spm.spatial.normalise.est.eoptions.nits = 16;
%     zatlabbatch{1,5}.spm.spatial.normalise.est.eoptions.reg = 1;
%     
%     zatlabbatch{1,6}.spm.spatial.normalise.write.subj.matname = cellstr([sour_file(1,1:end-6),'_seg_sn.mat']);
% %     zatlabbatch{1,6}.spm.spatial.normalise.write.subj.resample = cellstr(data_r);
%     zatlabbatch{1,6}.spm.spatial.normalise.write.subj.resample = data_r;
%     
%     zatlabbatch{1,6}.spm.spatial.normalise.write.roptions.preserve = 0;
%     zatlabbatch{1,6}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50;78 76 85];
%     zatlabbatch{1,6}.spm.spatial.normalise.write.roptions.vox = [3 3 3];
%     zatlabbatch{1,6}.spm.spatial.normalise.write.roptions.interp = 1;
%     zatlabbatch{1,6}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
%     zatlabbatch{1,6}.spm.spatial.normalise.write.roptions.prefix = 'w';
    
    %%%%%   spm12 %%%%%%%% 
    clear matlabbatch;
    
    run_name = dir(fullfile(func_path,subj{Si},'*run1*'));
    run1_path = fullfile(run_path,run_name.name)
    
    sour_file = [run1_path,'/meana',[run1_files(start_image).name],',1'];%Siemens
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = cellstr(sour_file);
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = cellstr(data_r);
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;

    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {fullfile(spm_path,'tpm','TPM.nii')};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                                 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = VoxSize;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
    
    spm_jobman('run',matlabbatch)
    
	%_____________________SMOOTH_____________________
    
    
    %%%%%   spm8 %%%%%%%% 
% % 	zatlabbatch{1,7}.spm.spatial.smooth.data = cellstr(data_w);
%     zatlabbatch{1,7}.spm.spatial.smooth.data = data_w;
%     zatlabbatch{1,7}.spm.spatial.smooth.fwhm = [5 5 5];
% 	zatlabbatch{1,7}.spm.spatial.smooth.dtype = 0;
%     zatlabbatch{1,7}.spm.spatial.smooth.im = 0;
%     zatlabbatch{1,7}.spm.spatial.smooth.prefix = 's';
    
    %%%%%   spm12 %%%%%%%% 
    clear matlabbatch;
    
    matlabbatch{1}.spm.spatial.smooth.data = cellstr(data_w);
    matlabbatch{1}.spm.spatial.smooth.fwhm = smooth_fwh;
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';

	spm_jobman('run',matlabbatch)

%end
end