clear all
sbj_no=[13,15,16]; % sbj number to analyze
subject_count = length(sbj_no);
addpath(genpath('C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\EEG Data\Codes')) % add analysis code to path
addpath(genpath('C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\eeglab2021.0'))
data_path = 'C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\EEG Data\eeg_data'; % name of the raw EEG file
analysis_path_marked = 'C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\EEG Data\Artefact_Files\Marked_Data';
analysis_path = 'C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\EEG Data\Analysis_Files\';

subject= struct();
for ss=1:subject_count
   
    fname = [num2str(sbj_no(ss)),'_MarkingComplete.set'];
    subject.EEG_marked{ss} = pop_loadset('filename',fname,'filepath',analysis_path_marked);
    
    fname = [num2str(sbj_no(ss)),'_referenced_epoched.set'];
    subject.EEG{ss} = pop_loadset('filename',fname,'filepath',analysis_path);
    subject.EEG{ss}.reject.rejmanual=subject.EEG_marked{ss}.reject.rejmanual;
    
end
%%




critical_event_markers=[31:34 36:39];

for sb=1:subject_count

subject.event_table{sb} = struct2table(subject.EEG{sb}.event);
triggers = (subject.event_table{sb}.type);
triggers= cell2mat(triggers);
triggers = str2num(triggers(:,3:end));

ind=ismember(triggers,critical_event_markers');
indices = find(ind==1);

subject.EEG{sb}.eventCodes = (triggers(indices))';


% sanity check 1: check if all event codes belong to the critical event markers
if any(~ismember(subject.EEG{sb}.eventCodes,critical_event_markers))
    error('a non-critical event marker was used in at least one of the epochs!')
end

% sanity check 2: check if the epoch number equals to trial number 
if length(subject.EEG{sb}.eventCodes)~=subject.EEG{sb}.trials
    error('Number of epochs do not match number of trials!')
end
end

%%

conditions=[31:34 36:39]; % 31,33,36,38: left. 32,34,37,39: right
for sb =1:subject_count
   
    
EEG=subject.EEG{sb};
 
condition_count=length(conditions)/2;

    
    
countCond=0; % reset a condition counter 
for icon_bilat = 1:2:length(conditions) 
    
    % successive numbers belong to the same condition, just different hemifieds.
    countCond=countCond+1; % update the condition counter
    cond=conditions(icon_bilat);
    
    % average across trials of this particular condition, this contains all
    % channel information
     subject.participant{sb}.left{countCond}.condata= EEG.data(:,:,EEG.eventCodes==cond & ~EEG.reject.rejmanual);
     subject.participant{sb}.right{countCond}.condata= EEG.data(:,:,EEG.eventCodes==cond+1& ~EEG.reject.rejmanual);
end



end

%%
chans_left = {'PO3';'P3';'P7';'PO7';'O1'};
chans_right = {'PO4';'P4';'P8';'PO8';'O2'};

  all_channels= struct2table(EEG.chanlocs).labels;
  
  for cc=1:size(chans_left,1)
   left_channels(1,cc)= string(chans_left{cc,1});
  
    right_channels(1,cc) = string(chans_right{cc,1});
  end
  for cc=1:size(all_channels,1)
      
   
    temp_channels(1,cc)= string(all_channels{cc,1});
      
  end

  left_index = zeros(size(chans_left,1),1);
  right_index = zeros(size(chans_right,1),1);
  
  for cc=1:size(chans_left,1)
      [~,left_index(cc,1)]= intersect(temp_channels, left_channels(1,cc),'stable');
      [~,right_index(cc,1)] = intersect(temp_channels,right_channels(1,cc),'stable');
  end
  
  %%
 
for sb=1:subject_count
   
    
    for cn =1:(length(conditions)/2)
        
         L_ipsi= subject.participant{sb}.left{cn}.condata(left_index,:,:);
         
         R_ipsi =subject.participant{sb}.right{cn}.condata(right_index,:,:);
      
         subject.participant{sb}.ipsi{cn}.L_ipsi=L_ipsi;
         subject.participant{sb}.ipsi{cn}.R_ipsi=R_ipsi;
         
         subject.participant{sb}.ipsi{cn}.IPSI_ALL = cat(3,L_ipsi,R_ipsi);
         
         L_contra =subject.participant{sb}.left{cn}.condata(right_index,:,:);% zeros(length(chans_right),subject.EEG{sb}.pnts,size(subject.participant{sb}.left{cn}.condata,3)); 
         R_contra= subject.participant{sb}.right{cn}.condata(left_index,:,:);
        
         subject.participant{sb}.contra{cn}.R_contra=R_contra;
         subject.participant{sb}.contra{cn}.L_contra=L_contra;
         
          subject.participant{sb}.contra{cn}.CONTRA_ALL = cat(3,L_contra,R_contra);
           
    end
    
    
end
 %%

 
min_freq =  4; % minimum frequency
max_freq = 50; % maximum frequency
num_frex = 20; % number of frequencies
min_cyc = 3;   % minimum number of cycles
max_cyc = 10;  % maximum number of cycles



for sb=1:subject_count


subject.time{sb} = -1:1/subject.EEG{sb}.srate:1; 
subject.frex{sb} = logspace(log10(min_freq),log10(max_freq),num_frex); 

subject.stds{sb} = logspace(log10(min_cyc),log10(max_cyc),num_frex)./(2*pi*subject.frex{sb});


end
%%


for sb=1:subject_count
    
    subject.n_wavelet{sb} = length(subject.time{sb});
    subject.half_of_wavelet_size{sb} = (subject.n_wavelet{sb}-1)/2;
    
    for cn =1:countCond
        
             
        % define convolution parameters

subject.n_data{sb}.IPSI{cn}.data_length=  subject.EEG{sb}.pnts*size(subject.participant{1,sb}.ipsi{1,cn}.IPSI_ALL,3);
subject.n_convolution{sb}.IPSI{cn}.convolution_length = subject.n_wavelet{sb}+subject.n_data{sb}.IPSI{cn}.data_length-1;
subject.n_conv_pow2{sb}.IPSI{cn}.pow2_length   = pow2(nextpow2(subject.n_convolution{sb}.IPSI{cn}.convolution_length ));


data = subject.participant{1,sb}.ipsi{1,cn}.IPSI_ALL;
data_to_consider_ipsi = reshape(data,[size(data,1),size(data,2)*size(data,3)]);
subject.eeg_fft{sb}.IPSI{cn}.eegfft =fft(data_to_consider_ipsi,subject.n_conv_pow2{sb}.IPSI{cn}.pow2_length ,2);
 
subject.n_data{sb}.CONTRA{cn}.data_length =  subject.EEG{sb}.pnts*size(subject.participant{1,sb}.contra{1,cn}.CONTRA_ALL,3);
subject.n_convolution{sb}.CONTRA{cn}.convolution_length = subject.n_wavelet{sb}+subject.n_data{sb}.CONTRA{cn}.data_length-1;
subject.n_conv_pow2{sb}.CONTRA{cn}.pow2_length= pow2(nextpow2(subject.n_convolution{sb}.CONTRA{cn}.convolution_length));

 data = subject.participant{1,sb}.contra{1,cn}.CONTRA_ALL;
 data_to_consider_contra = reshape(data,[size(data,1),size(data,2)*size(data,3)]);
 subject.eeg_fft{sb}.CONTRA{cn}.eegfft =fft(data_to_consider_contra,subject.n_conv_pow2{sb}.CONTRA{cn}.pow2_length ,2);

        
    end
    
    
    
    
 
end
 
 %%
 pChan_pSubject_All =zeros(5,1750,20,4,3);
for sb = 1:subject_count
     
    baseidx = dsearchn(subject.EEG{sb}.times',[-500 -200]');
   
    % get the particular frequency and standard deviation for this freq
  
            
    for cn =1:countCond
        
        for fi=1:num_frex
            
   
    f = subject.frex{sb}(fi);
    s = subject.stds{sb}(fi);
    
    % make the wavelet of this frequency 
    wavelet = exp(1i*2*pi*f.*(subject.time{sb})).* exp(-(subject.time{sb}).^2./(2*s^2));
            
    %IPSI
    wavelet_fft_ipsi = fft(sqrt(1/(s*sqrt(pi))) * wavelet , subject.n_conv_pow2{sb}.IPSI{cn}.pow2_length );
    fft_matrix_ipsi = wavelet_fft_ipsi.*subject.eeg_fft{1,sb}.IPSI{1,cn}.eegfft;
    
    
 
    eegconv = ifft(fft_matrix_ipsi,subject.n_conv_pow2{sb}.IPSI{cn}.pow2_length,2);
    eegconv = eegconv(:,1:subject.n_convolution{sb}.IPSI{cn}.convolution_length );
    eegconv = eegconv(:,subject.half_of_wavelet_size{sb}+1:end-subject.half_of_wavelet_size{sb});
   
    [chan,data_len,cond_count] = size (subject.participant{1,sb}.ipsi{1,cn}.IPSI_ALL);
    eegconv=reshape(eegconv,[chan,data_len,cond_count]);
   
    temppower = mean(abs(eegconv).^2,3);
    
    subject.participant{1,sb}.ipsi{1,cn}.eegpower(:,:,cn,fi) = 10*log10(temppower./mean(temppower(:,baseidx(1):baseidx(2)),2));
    
    % CONTRA
    
    wavelet_fft_contra = fft(sqrt(1/(s*sqrt(pi))) * wavelet , subject.n_conv_pow2{sb}.CONTRA{cn}.pow2_length );
    fft_matrix_contra = wavelet_fft_contra.*subject.eeg_fft{1,sb}.CONTRA{1,cn}.eegfft;
    
    
 
    eegconv = ifft(fft_matrix_contra,subject.n_conv_pow2{sb}.CONTRA{cn}.pow2_length,2);
    eegconv = eegconv(:,1:subject.n_convolution{sb}.CONTRA{cn}.convolution_length );
    eegconv = eegconv(:,subject.half_of_wavelet_size{sb}+1:end-subject.half_of_wavelet_size{sb});
   
    [chan,data_len,cond_count] = size (subject.participant{1,sb}.contra{1,cn}.CONTRA_ALL);
    eegconv=reshape(eegconv,[chan,data_len,cond_count]);
   
    temppower = mean(abs(eegconv).^2,3);
    
    subject.participant{1,sb}.contra{1,cn}.eegpower(:,:,cn,fi) = 10*log10(temppower./mean(temppower(:,baseidx(1):baseidx(2)),2));
    
    
    
 
            
        end
        
      
       subject.participant{1,sb}.contra_ipsi{1,cn}= squeeze(subject.participant{1,sb}.contra{1,cn}.eegpower(:,:,cn,:))-squeeze(subject.participant{1,sb}.ipsi{1,cn}.eegpower(:,:,cn,:));
        
        pChan_pSubject_All(:,:,:,cn,sb) = subject.participant{1,sb}.contra_ipsi{1,cn};
    end
    

       
    
end
 
 


%%



%%
%%
legatns = ["High-Prec Load 1","High-Prec Load 2","Low-Prec Load 1","Low-Prec Load 2"];


    
   for cn=1:countCond
    
       pos= cn;

subplot(2,2,pos)
contourf(EEG.times,subject.frex{sb},(squeeze( mean(mean(pChan_pSubject_All(:,:,:,cn,:),1),5)))',80,'linecolor','none')
set(gca,'clim',[-1 1],'xlim',[500 1000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
legend(legatns(cn))
tt = ['Logarithmic frequency scaling for All Subjects for Condition  ',num2str(cn)];
title(tt)

   end

 %% STATISTICAL ANALYSIS
 
 
 
% define analysis specs
toi =[400 1000] ;%[400 1000]; % times of interest
foi = [8, 12];

toi_indices = dsearchn(subject.EEG{1}.times',toi'); 
toi_indices = toi_indices(1):12:toi_indices(2);
foi_indices = dsearchn(subject.frex{1}',foi');
foi_indices = foi_indices(1):foi_indices(2);


power_pSbj_pCond = squeeze(mean(mean(mean(pChan_pSubject_All(:,toi_indices,foi_indices,:,:),1),2),3));

%%

power_pCond = mean(power_pSbj_pCond,2); % average across subjects
sem_pCond = sqrt(var(power_pSbj_pCond,0,2)./subject_count); % standard error of the mean

sem_pCondDiff(1:2) = std(power_pSbj_pCond(2,:)-power_pSbj_pCond(1,:))./sqrt(subject_count);
sem_pCondDiff(3:4) = std(power_pSbj_pCond(4,:)-power_pSbj_pCond(3,:))./sqrt(subject_count);

% plot properties
cond_colors = [53/255, 97/255, 143/255; 53/255, 97/255, 143/255; ...
    235/255, 84/255, 199/255;235/255, 84/255, 199/255]; % rgb values
cond_labels = {'High Prec. Load 1','High Prec. Load 2','Low Prec. Load 1','Low Prec. Load 2' };
cond_label_groups = {'High Prec.','Low Prec' };
legend_names = {'Load 1','Load 2','Load 1','Load 2'};



% plot
figure
for iCond=1:countCond
    hb = bar(iCond,power_pCond(iCond)); % bar plot
    hold on
    set(hb(1),'facecolor',cond_colors(iCond,:)); % bar color
    % make load 1 bars half transparent and dashed
    if iCond==1 || iCond==3
        set(hb(1),'LineStyle',':','LineWidth',2);
        hb.FaceAlpha = 0.5;
    else
        set(hb(1),'LineWidth',2);
    end
end

% add legend names (needs to be done before error bars)
legend(legend_names) 

% plot error bars
for iCond=1:countCond
    errorbar(iCond,power_pCond(iCond),sem_pCondDiff(iCond),'k.', 'HandleVisibility','off','LineWidth',1.5);
end


% more properties
set(gca,'XTick',[1.5 3.5]) % a tick in the middle of consecutive bars
set(gca,'XTickLabel',cond_label_groups,'FontSize',12) % condition labels
set(gca, 'YDir','reverse') % negative up
xlabel('Precision requirement','FontSize',14,'FontWeight','bold') % x-label
ylabel('Alpha (8-12 Hz) Power (dB)','FontSize',14,'FontWeight','bold') % y-label



