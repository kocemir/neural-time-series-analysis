%% Artefact Rejection Indexe

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
     subject.participant{sb}.right{countCond}.condata= EEG.data(:,:,EEG.eventCodes==cond+1 & ~EEG.reject.rejmanual);
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
         
         L_contra =subject.participant{sb}.left{cn}.condata(right_index,:,:);% zeros(length(chans_right),subject.EEG{sb}.pnts,size(subject.participant{sb}.left{cn}.condata,3)); 
         R_contra= subject.participant{sb}.right{cn}.condata(left_index,:,:);
        
         subject.participant{sb}.contra{cn}.R_contra=R_contra;
         subject.participant{sb}.contra{cn}.L_contra=L_contra;
         
         
         
        
    end
    
    
end
 %%
 
x_axis_limit = [-200 1000]; % in ms
y_axis_limit = [-20 20];
    
    for sb =1:subject_count
        
        
        for cn=1:(length(conditions)/2)
            
            LR_ipsi = cat(3, subject.participant{sb}.ipsi{cn}.L_ipsi,...
                      subject.participant{sb}.ipsi{cn}.R_ipsi);
                  
           
            LR_contra = cat(3, subject.participant{sb}.contra{cn}.R_contra,... 
            subject.participant{sb}.contra{cn}.L_contra);
             
            subject.participant{sb}.contra_ipsi{cn}.difference = squeeze(mean(LR_contra,3)-mean(LR_ipsi,3));
        
        figure
        subplot(2,4,1)
        plot(EEG.times,squeeze(mean(subject.participant{sb}.ipsi{cn}.L_ipsi,3)),'LineWidth',2) % 
        set(gca,'xlim',x_axis_limit,'ylim',y_axis_limit,'ydir','reverse')
        tt = [' Subject ', num2str(sbj_no(sb)), ' Condition  ' num2str(cn), ' Left Ipsi'];
        title(tt);
        legend(left_channels);
        
        hold on

plot(get(gca,'xlim'),[0 0],'k') % horizontal reference line
plot([0 0],get(gca,'ylim'),'k:') % vertical reference line (dashed)
        
 hold off


        subplot(2,4,2)
        plot(EEG.times,squeeze(mean(subject.participant{sb}.ipsi{cn}.R_ipsi,3)),'LineWidth',2) % 
        set(gca,'xlim',x_axis_limit,'ylim',y_axis_limit,'ydir','reverse')
        tt = [' Subject ', num2str(sbj_no(sb)), 'Condition ' num2str(cn), '  Right Ipsi'];
        title(tt);
        legend(right_channels);
        
        hold on

plot(get(gca,'xlim'),[0 0],'k') % horizontal reference line
plot([0 0],get(gca,'ylim'),'k:') % vertical reference line (dashed)
        hold off
        
        subplot(2,4,3)
        plot(EEG.times,squeeze(mean(LR_ipsi,3)),'LineWidth',2) % 
        set(gca,'xlim',x_axis_limit,'ylim',y_axis_limit,'ydir','reverse')
        tt = [' Subject ', num2str(sbj_no(sb)), ' Condition ' num2str(cn), ' Average Ipsi'];
        title(tt);
           hold on

plot(get(gca,'xlim'),[0 0],'k') % horizontal reference line
plot([0 0],get(gca,'ylim'),'k:') % vertical reference line (dashed)
        hold off
        
        subplot(2,4,4)
        plot(EEG.times,squeeze(mean((mean(LR_ipsi,3)),1)),'LineWidth',2) % 
        set(gca,'xlim',x_axis_limit,'ylim',y_axis_limit,'ydir','reverse')
        tt = [' Subject ', num2str(sbj_no(sb)), ' Condition ' num2str(cn), ' Average Ipsi for All Channels'];
        title(tt);
             
        hold on

plot(get(gca,'xlim'),[0 0],'k') % horizontal reference line
plot([0 0],get(gca,'ylim'),'k:') % vertical reference line (dashed)
        hold off
        
        
        
        subplot(2,4,5)
        plot(EEG.times,squeeze(mean(subject.participant{sb}.contra{cn}.L_contra,3)),'LineWidth',2) % 
        set(gca,'xlim',x_axis_limit,'ylim',y_axis_limit,'ydir','reverse')
        tt = [' Subject ', num2str(sbj_no(sb)), ' Condition ' num2str(cn), ' Left Contra'];
        title(tt);
        legend(right_channels);
        
        
              hold on

plot(get(gca,'xlim'),[0 0],'k') % horizontal reference line
plot([0 0],get(gca,'ylim'),'k:') % vertical reference line (dashed)
        hold off
        
        subplot(2,4,6)
        plot(EEG.times,squeeze(mean(subject.participant{sb}.contra{cn}.R_contra,3)),'LineWidth',2) % 
        set(gca,'xlim',x_axis_limit,'ylim',y_axis_limit,'ydir','reverse')
        tt = [' Subject ', num2str(sbj_no(sb)), 'Condition ' num2str(cn), '  Right Contra'];
        title(tt);
        legend(left_channels)
        
              hold on

plot(get(gca,'xlim'),[0 0],'k') % horizontal reference line
plot([0 0],get(gca,'ylim'),'k:') % vertical reference line (dashed)
        hold off
        
        
        subplot(2,4,7)
        plot(EEG.times,squeeze(mean(LR_contra,3)),'LineWidth',2) % 
        set(gca,'xlim',x_axis_limit,'ylim',y_axis_limit,'ydir','reverse')
        tt = [' Subject ', num2str(sbj_no(sb)), 'Condition ' num2str(cn), ' Average Contra'];
        title(tt);
        
              hold on

plot(get(gca,'xlim'),[0 0],'k') % horizontal reference line
plot([0 0],get(gca,'ylim'),'k:') % vertical reference line (dashed)
        hold off
        
         
        subplot(2,4,8)
        plot(EEG.times,squeeze(mean((mean(LR_contra,3)),1)),'LineWidth',2) % 
        set(gca,'xlim',x_axis_limit,'ylim',y_axis_limit,'ydir','reverse')
        tt = [' Subject ', num2str(sbj_no(sb)), ' Condition ' num2str(cn), ' Average Contra for All Channels'];
        title(tt);
        
              hold on

plot(get(gca,'xlim'),[0 0],'k') % horizontal reference line
plot([0 0],get(gca,'ylim'),'k:') % vertical reference line (dashed)
        hold off
        
        
        sgtitle([' Plot for Subjet ' num2str(sbj_no(sb))]);
       
        end
        
         end
%% Contra-Ipsi Difference

y_axis_limit=[-5 5];

pChan_pSubject_pCond_All = zeros(5,1750,4,3);
for sb= 1:subject_count
     

    
   for cn = 1:countCond
       subplot(1,3,sb);
     
       plot(subject.EEG{sb}.times,squeeze(mean(subject.participant{1,sb}.contra_ipsi{1,cn}.difference,1)));
       set(gca,'xlim',x_axis_limit,'ylim',y_axis_limit,'ydir','reverse')
     
  
      tt=[' Subject ',  num2str(sbj_no(sb)), ' Contra-Ipsi Difference'];
      title(tt);
       hold on
      
       pChan_pSubject_pCond_All(:,:,cn,sb)= subject.participant{1,sb}.contra_ipsi{1,cn}.difference;
       
   end
    
    legend(["High Precision 1","High Precision 2","Low Precision 1","Low Precision 2"])
    hold off  
end

%% All Subject per Condition

       y_axis_limit = [-3 3];
       plot(subject.EEG{sb}.times,squeeze(mean((mean(pChan_pSubject_pCond_All(:,:,:,:),1)),4)));
       set(gca,'xlim',x_axis_limit,'ylim',y_axis_limit,'ydir','reverse')
       legend(["High Precision 1","High Precision 2","Low Precision 1","Low Precision 2"])
  
       tt=" ALL SUBJECT AVERAGE PER CONDITION ";
       title(tt);
    
 