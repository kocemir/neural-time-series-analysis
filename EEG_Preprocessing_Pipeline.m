

%% THIS PIECE OF CODE IS TO PREPARE DATA TO FURTHER ANALYSES
  
% Re-referencing is a little bit confusing. Therefore, please refer to
% EEGLAB>Re-referecing tutorial
%% SET ALL FOLDERS AND PATHS

clear all;clc;
addpath(genpath('C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\EEG Data\Codes')) % add analysis code to path
addpath(genpath('C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\eeglab2021.0'))
data_path = 'C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\EEG Data\eeg_data'; % name of the raw EEG file
analysis_path = 'C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\EEG Data\Analysis_Files\';

%% LOAD RAW EEG DATA ( VHDR FILES in BVA ), REFERENCE AND FILTER IT, EPOCH AND REMOVE BASELINE


sbj_no = [13,15,16];  % Update this part if 
subject_count = length(sbj_no);
subject = struct();  % Initialize a struct to store .set file (still raw but, '.set' is EEGLAB's storing format )
[ALLEEG,EEG,CURRENTSET,~] = eeglab;   % Initialize EEGLAB

for sb=1:subject_count

    
filename = ['priVWME3_0',num2str(sbj_no(sb)),'.vhdr']; % name root of the raw EEG file
[EEG,~] = pop_loadbv(data_path, filename); % load the raw EEG file but do not store history knowledge
% EEG = eeg_checkset( EEG ); % check the dataset consistency    
%[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG,EEG,sb);  % Store current and all EEG data


save_filename =  [num2str(sbj_no(sb)),'.set']; % name root of the raw EEG file 
EEG = pop_saveset( EEG, 'filename',save_filename,'filepath',analysis_path); % save dataset   

   
% REFERENCE--> Retain TP10 back and reference TP10 and TP9 to all channels
EEG=pop_chanedit(EEG, 'append',35,'changefield',{36,'labels','TP10'},'changefield',{36,'theta','108'},'changefield',{36,'radius','0.627'},'convert',{'topo2all'});
% [ALLEEG ,EEG] = eeg_store(ALLEEG, EEG, sb);
% EEG = eeg_checkset( EEG );
EEG=pop_chanedit(EEG, 'setref',{'1:31','TP10'});
% [ALLEEG ,EEG] = eeg_store(ALLEEG, EEG, sb);
% EEG = eeg_checkset( EEG );
EEG = pop_reref( EEG, [],'refloc',struct('labels',{'TP10'},'sph_radius',{1},'sph_theta',{-108},'sph_phi',{-22.86},'theta',{108},'radius',{0.627},'X',{-0.28475},'Y',{-0.87636},'Z',{-0.38848},'type',{''},'ref',{''},'urchan',{[]},'datachan',{0},'sph_theta_besa',{112.86},'sph_phi_besa',{-18}),'exclude',[32:35] );
% [ALLEEG ,EEG,CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
% EEG = eeg_checkset( EEG );
EEG = pop_reref( EEG, [10 36] ,'exclude',[32:35] );
save_filename =  [num2str(sbj_no(sb)),'_referenced.set']; % name root of the raw EEG file 
EEG = pop_saveset( EEG, 'filename',save_filename,'filepath',analysis_path); % save dataset  


% FILTERING
save_filename =  [num2str(sbj_no(sb)),'_referenced.set']; % name root of the raw EEG file 
EEG = pop_eegfiltnew(EEG, 'locutoff',0.5,'hicutoff',40,'plotfreqz',0);
EEG = pop_saveset( EEG, 'filename',save_filename,'filepath',analysis_path); 


% EPOCHING
EEG = pop_epoch( EEG, {  'S 31'  'S 32'  'S 33'  'S 34'  'S 36'  'S 37'  'S 38'  'S 39'  }, [-1.5  2], 'epochinfo', 'yes');
EEG = pop_rmbase( EEG, [-200 0] ,[]); 
save_filename =  [num2str(sbj_no(sb)),'_referenced_epoched.set']; % name root of the raw EEG file
EEG = pop_saveset( EEG, 'filename',save_filename,'filepath',analysis_path); 

subject.EEG{sb}=EEG;

end

