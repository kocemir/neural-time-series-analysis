%% Psy 510 Spring 2021 Lecture code
% Computing time-frequency power in one channel in one dataset across a
% range of frequencies
%
% Based on Mike Cohen's code that accompanies the book, titled 
% "Analyzing Neural Time Series Data" (MIT Press)
%
% Eren Günseli
% Sabancı University, 2021
%%
clear; clc;

sbj_no=[13,15,16]; % sbj number to analyze
subject_count=length(sbj_no);
addpath(genpath('C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\EEG Data\Codes')) % add analysis code to path
addpath(genpath('C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\eeglab2021.0'))
data_path = 'C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\EEG Data\eeg_data'; % name of the raw EEG file
analysis_path = 'C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\EEG Data\Analysis_Files\';



subject= struct();
for ss=1:subject_count
   
    fname = [num2str(sbj_no(ss)),'_referenced_epoched.set'];
    subject.EEG{ss} = pop_loadset('filename',fname,'filepath',analysis_path);
end


%%
critical_event_markers = [31:34 36:39]; % event markers used for epoching

for sb=1:subject_count
 
EEG=subject.EEG{sb};
    
EEG.eventCodes=nan(1,length(EEG.epoch)); % preallocate
for iEpoch = 1 : length(EEG.epoch) % loop through epochs
    found_critical_event=0;
    while found_critical_event==0
        for event_type=1:length(EEG.epoch(iEpoch).eventtype) % loop through all events in this epoch
            for event_marker=critical_event_markers % loop through event markers of interest
                s = EEG.epoch(iEpoch).eventtype{event_type}(3:end); % get characters 3 to end (i.e. the digits (first three characters are irrelevant)
                s = str2double(s);             % covert string to double
                if any(s==critical_event_markers) % is it a critical event marker?
                    EEG.eventCodes(1,iEpoch) = s;      % save the event code structure
                    found_critical_event=1; % we found the critical event in this epoch. let's move on to the next epoch.
                end
            end
        end
    end
end

% sanity check 1: check if all event codes belong to the critical event markers
if any(~ismember(EEG.eventCodes,critical_event_markers))
    error('a non-critical event marker was used in at least one of the epochs!')
end

% sanity check 2: check if the epoch number equals to trial number 
if length(EEG.eventCodes)~=EEG.trials
    error('Number of epochs do not match number of trials!')
end

subject.EEG{sb}=EEG;
end
%%

% parameters
conditions=[31:34 36:39]; % 31,33,36,38: left. 32,34,37,39: right

% preallocate
nbconds_bilat=length(conditions)/2;


for sb =1:subject_count
    
    EEG=subject.EEG{sb};
% start condition loop
   countCond=0; % reset a condition counter 
for icon_bilat = 1:2:length(conditions) % we will go in steps of two because 
    % successive numbers belong to the same condition, just different hemifieds.
    % eg. low-precision set size one left vs right See lines 15-23
    countCond=countCond+1; % update the condition counter
    cond=conditions(icon_bilat);
    
    ind = EEG.eventCodes==cond | EEG.eventCodes==cond+1;
    index= find(ind==1);
    subject.condition{sb}.indeks{countCond} = index;
    
   
end
end




% %%
% 
% % definitions, selections...
% chan2use = {"Fz","F4","F3"}; % channels to analyze
% chanlocs_table = struct2table(EEG.chanlocs);

%ind = [2,3,28];

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
%


%%
for sb=1:subject_count


% define convolution parameters
subject.n_wavelet{sb}            = length(subject.time{sb});
subject.n_data{sb}               = subject.EEG{sb}.pnts*subject.EEG{sb}.trials;
subject.n_convolution{sb}        = subject.n_wavelet{sb}+subject.n_data{sb}-1;

% find the closest number larger than n_convolution that is a factor of 2
% (to decrease computation time, see Cohen's book, section 13.9)
subject.n_conv_pow2{sb}          = pow2(nextpow2(subject.n_convolution{sb}));

% calculate hald the wavelet size to be used for trimming later
subject.half_of_wavelet_size{sb} = (subject.n_wavelet{sb}-1)/2;



% I USED 17th channel data only

data_to_consider=subject.EEG{sb}.data(:,:,:);
data_eeg = reshape(data_to_consider,[size(data_to_consider,1), subject.EEG{sb}.pnts*subject.EEG{sb}.trials]);

subject.eeg_fft{sb}=fft(data_eeg,subject.n_conv_pow2{sb},2);

end

%%
% initialize

for sb = 1:subject_count


eegpower = zeros(size(data_to_consider,1),subject.EEG{sb}.pnts,countCond,num_frex); % frequencies X time X trials

% find time point indices of the baseline boundaries (to be used to
% baseline normalization below)
baseidx = dsearchn(subject.EEG{sb}.times',[-500 -200]');


% loop through frequencies and compute synchronization

for fi=1:num_frex
    
    % get the particular frequency and standard deviation for this freq
    f = subject.frex{sb}(fi);
    s = subject.stds{sb}(fi);
    
    % make the wavelet of this frequency 
    wavelet = exp(1i*2*pi*f.*(subject.time{sb})).* exp(-(subject.time{sb}).^2./(2*s^2));
 
    %% Check if wavelets taper to zero at the edges
%     plot_wavelets = 1; % Do you want to plot your wavelets?
%     if plot_wavelets == 1
%         if fi==1
%             figure
%         end
%         subplot(5,5,fi)
%         plot((subject.time{sb}),real(wavelet)) % plots the real component
%         close all
%     end
%     
    
  
    wavelet_fft = fft(sqrt(1/(s*sqrt(pi))) * wavelet , subject.n_conv_pow2{sb});
    
    
    fft_matrix = wavelet_fft.*subject.eeg_fft{sb};
    eegconv = ifft(fft_matrix,subject.n_conv_pow2{sb},2);
    eegconv = eegconv(:,1:subject.n_convolution{sb});
    eegconv = eegconv(:,subject.half_of_wavelet_size{sb}+1:end-subject.half_of_wavelet_size{sb});
   
    eegconv=reshape(eegconv,[size(data_to_consider,1),subject.EEG{sb}.pnts,subject.EEG{sb}.trials]);
  
    
    for cn =1:countCond
         
     temppower = mean(abs(eegconv(:,:,subject.condition{1,sb}.indeks{1,cn})).^2,3); 
     eegpower(:,:,cn,fi) = 10*log10(temppower./mean(temppower(:,baseidx(1):baseidx(2)),2));
        
    end
   
  

end
 
subject.eeg_power{sb}=eegpower;

end

%%
% 
% for sb =1:subject_count
%     
%     for cn = 1:nbconds_bilat
%         
%         ind_array=subject.condition{1,sb}.indeks{1,cn};
%        
%         condition_eegpower = subject.eeg_power{sb}(:,:,ind_array,:);
%         subject.participant{sb}.condition_power{cn}= (squeeze(mean(condition_eegpower,3)));
%         
%     end
%     
%     
% end

%%

subject_all = nan(size(data_to_consider,1),EEG.pnts,countCond,num_frex,subject_count);
for sb=1:subject_count
   
    figure
    
   for cn=1:nbconds_bilat
    
       pos= cn;

subplot(2,2,pos)
contourf(subject.EEG{sb}.times,subject.frex{sb},squeeze(subject.eeg_power{sb}(22,:,cn,:))',120,'linecolor','none')
set(gca,'clim',[-3 3],'xlim',[-200 1000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
tt = ['Logarithmic frequency scaling for Condition',num2str(cn)];
title(tt)


 subject_all(:,:,:,:,sb) = subject.eeg_power{sb};

   end
   

 sgtitle(['Time-Frequency Plot for Subjet ' num2str(sbj_no(sb))]);
end
%% Average all subjects per condition 

chan2look = [  12  18 ];
legatns = ["High-Prec Load 1","High-Prec Load 2","Low-Prec Load 1","Low-Prec Load 2"];
for cn=1:countCond

    subplot(2,2,cn)
 contourf(subject.EEG{sb}.times,subject.frex{sb},(squeeze(mean(mean(subject_all(chan2look,:,cn,:,:),5),1)))',120,'linecolor','none')
set(gca,'clim',[-3 3],'xlim',[-500 1300],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
legend(legatns(cn))
tt = ['Logarithmic frequency scaling for All Subjects for Condition  ',num2str(cn)];
title(tt)


end

%% Statistical Analysis


% define analysis specs
toi = [400 1000]; % times of interest
foi = [8.88081098856547, 11.5855548835346 ];

toi_indices = dsearchn(subject.EEG{1}.times',toi'); 
toi_indices = toi_indices(1):toi_indices(2);
foi_indices = dsearchn(subject.frex{1}',foi');
foi_indices = foi_indices(1):foi_indices(2);


power_pSbj_pCond = squeeze(mean(mean(mean(subject_all(chan2look,toi_indices,:,foi_indices,:),1),2),4));

%%

power_pCond = mean(power_pSbj_pCond,2); % average across subjects
sem_pCond = sqrt(var(power_pSbj_pCond,0,2)./subject_count); % standard error of the mean


% standard error of the mean (difference of load 2 - load 1)
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


