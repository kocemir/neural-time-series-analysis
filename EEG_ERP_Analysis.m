%% LOAD PREVIOUSLY PRE-PROCESSED .SET DATA
clear
subject_count=3;
sbj_no=[13,15,16]; % sbj number to analyze

addpath(genpath('C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\EEG Data\Codes')) % add analysis code to path
addpath(genpath('C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\eeglab2021.0'))
data_path = 'C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\EEG Data'; % name of the raw EEG file
analysis_path = 'C:\Users\kocemirhan\Desktop\SPRING 2020-2021\PSY510\EEG Data\Analysis_Files\';

subject= struct();
for ss=1:subject_count
   
    fname = [num2str(sbj_no(ss)),'_referenced_epoched.set'];
    subject.EEG{ss} = pop_loadset('filename',fname,'filepath',analysis_path);
end
%%  Store event codes for each

tic
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
toc
%% HERE, CONDITIONWISE DATA ORGANIZATION

conditions=[31:34 36:39]; % 31,33,36,38: left. 32,34,37,39: right

for sb =1:subject_count
   
    
EEG=subject.EEG{sb};
 
condition_count=length(conditions)/2;
data_pChan_pTime_pCond=nan(EEG.nbchan,EEG.pnts,condition_count);
    
    
countCond=0; % reset a condition counter 
for icon_bilat = 1:2:length(conditions) 
    
    % successive numbers belong to the same condition, just different hemifieds.
    countCond=countCond+1; % update the condition counter
    cond=conditions(icon_bilat);
    
    % average across trials of this particular condition, this contains all
    % channel information
    subject.data_pChan_pTime_pCond{sb}(:,:,countCond)= mean(EEG.data(:,:,EEG.eventCodes==cond | EEG.eventCodes==cond+1),3); 
end

end
%% PLOT RANDOM TRIAL DATA

% parameters
x_axis_limit = [-200 1200]; % in ms
num_trials2plot = 12; % how many trials to plot
channels=struct2table(subject.EEG{1}.chanlocs).labels;
which_channel_to_plot = 'Cz'; % write the name of the channel
% and find the index (channel number) of that label
channel_index = strmatch(which_channel_to_plot,channels);

for sb=1:subject_count

    EEG=subject.EEG{sb};
    figure
    set(gcf,'Name',[ 'Subject_' num2str(sbj_no(sb)) '_' num2str(num_trials2plot) ' random trials from channel ' EEG.chanlocs(channel_index).labels ])

for i=1:num_trials2plot
    
    % figure out how many subplots we need
    subplot(ceil(num_trials2plot/ceil(sqrt(num_trials2plot))),ceil(sqrt(num_trials2plot)),i)
    % pick a random trial (using randsample, which is in the stats toolbox)
    random_trial_to_plot = randsample(EEG.trials,1);
  
    % plot trial and specify x-axis and title
    plot(EEG.times,squeeze(EEG.data(channel_index,:,random_trial_to_plot)));
    set(gca,'xlim',x_axis_limit,'ytick',[])
    title([ 'Trial ' num2str(random_trial_to_plot) ])
end

end
%% HERE, AVERAGE ALL THE TRIALS AND SEE HOW YOU SMOOTHEN THE SIGNAL AND REDUCE SIGNAL

for sb =1:subject_count
    
EEG= subject.EEG{sb};

figure
% plot all trials
plot(EEG.times,squeeze(EEG.data(channel_index,:,:)));

hold on
% plot ERP (simply the average time-domain signal)
plot(EEG.times,squeeze(mean(EEG.data(channel_index,:,:),3)),'g','linew',5)
set(gca,'xlim',[-200 1200],'ylim',[-60 60])
set(gcf,'Name',[ 'Subject ' num2str(sbj_no(sb)) '  and channel  '  EEG.chanlocs(channel_index).labels ])
title(['Subject ' num2str(sbj_no(sb)) ' All Data ']);



% now plot only the ERP
figure
plot(EEG.times,squeeze(mean(EEG.data(channel_index,:,:),3))) % Note the "3" as second input to "mean"; this takes the average of the 3rd dimension.

% plot lines indicating baseline activity and stim onset
hold on
plot(get(gca,'xlim'),[0 0],'k')
plot([0 0],get(gca,'ylim'),'k:')

% add axis labels and title
xlabel('Time (ms)')
ylabel('\muV') % note that matlab interprets "\mu" as the Greek character for micro
title(['Subject' num2str(sbj_no(sb)) '__ERP (average of ' num2str(EEG.trials) ' trials) from electrode ' EEG.chanlocs(channel_index).labels ])

% plot upside down, following ERP convention
%set(gca) 

% below is some advanced but flexible code to change the x-axis label
set(gca,'xlim',[-200 1200])
xticklabel=cellstr(get(gca,'xticklabel'));
xticklabel{str2double(xticklabel)==0}='stim';
set(gca,'xticklabel',xticklabel)

hold off
end 
%% PLOT AVERAGE FOR EACH CONDITION

% parameters



for sb = 1: subject_count
    
which_channel_to_plot = 'Cz'; % in string
x_axis_limit = [-200 1200]; % in ms
y_axis_limit = [-15 15];

channels=struct2table(subject.EEG{1}.chanlocs).labels;
which_channel_to_plot = 'Cz'; % write the name of the channel
channel_index = strmatch(which_channel_to_plot,channels); % Find the index of this channel
 
subject.ERP_pTime_pCond{sb}=squeeze(mean(subject.data_pChan_pTime_pCond{1,sb}(channel_index,:,:),1));

% plot
figure
plot(EEG.times,subject.ERP_pTime_pCond{sb},'LineWidth',2) % Note the "3" as second input to "mean"; this takes the average of the 3rd dimension.
set(gca,'xlim',x_axis_limit,'ylim',y_axis_limit)
hold on

plot(get(gca,'xlim'),[0 0],'k') % horizontal reference line
plot([0 0],get(gca,'ylim'),'k:') % vertical reference line (dashed)
xticklabel=cellstr(get(gca,'xticklabel')); % read the x tick labels
xticklabel{str2double(xticklabel)==0}='Memory'; % convert 0 to 'Memory'
set(gca,'xticklabel',xticklabel)
set(gcf,'Name',[ 'Subject ' num2str(sbj_no(sb)) '  and channel  '  EEG.chanlocs(channel_index).labels ])
title([ 'Channel ', EEG.chanlocs(channel_index).labels ,' for Subject ', num2str(sbj_no(sb))])

xlabel('Time (ms)')
ylabel('\muV') % note that matlab interprets "\mu" as the Greek character for micro
legend('High-Prec Load 1','High-Prec Load 2','Low-Prec Load 1','Low-Prec Load 2','Location','southeast')

end



% At this stage, we have created ERP_pTime_pCond 3D matrixes to store  averaged trials for each different condition

%% MULTIPLE SUBJECT ANALYES

sizes = [size(subject.ERP_pTime_pCond{1}),subject_count];
all_subject_ERP_pTime_pCond= nan(sizes);


for sb=1:subject_count
    
    all_subject_ERP_pTime_pCond(:,:,sb)=subject.ERP_pTime_pCond{sb};
    
end
avg_all_subject_ERP_pTime_pCond = squeeze(mean(all_subject_ERP_pTime_pCond,3));

figure
plot(EEG.times,avg_all_subject_ERP_pTime_pCond,'LineWidth',3) % Note the "3" as second input to "mean"; this takes the average of the 3rd dimension.
set(gca,'xlim',x_axis_limit,'ytick',[])
hold on



plot(get(gca,'xlim'),[0 0],'k') % horizontal reference line
plot([0 0],get(gca,'ylim'),'k:') % vertical reference line (dashed)
xticklabel=cellstr(get(gca,'xticklabel')); % read the x tick labels
xticklabel{str2double(xticklabel)==0}='Memory'; % convert 0 to 'Memory'
set(gca,'xticklabel',xticklabel)
title([ 'Channel ' which_channel_to_plot ' average of all subject ']);

xlabel('Time (ms)')
ylabel('\muV') % note that matlab interprets "\mu" as the Greek character for micro
legend('High-Prec Load 1','High-Prec Load 2','Low-Prec Load 1','Low-Prec Load 2','Location','southeast')

%% TOPOGRAPHIC PLOTS

channel_table = struct2table(EEG.chanlocs);
channel_table = channel_table(1:30,:);  %% Do not take into account the last 4 channels
EEG.chanlocs_new = table2struct(channel_table);

times2plot =400:900;
for sb=1:subject_count

    topoLim = [-6 6];
    figure

for cc= 1:countCond
    
    subplot(2,2,cc);
    avg_erp=mean(subject.data_pChan_pTime_pCond{sb}(:,times2plot,cc),2); % average within this time range
    topoplot(avg_erp, EEG.chanlocs_new,'electrodes','labels','plotrad',.6,'whitebk','on','maplimits',topoLim,'emarker',{'.','r',5,1});
    set(gcf,'numbertitle','off','name',[num2str(times2plot(1)) '-' num2str(times2plot(end)) ' ms']) % title
    title([' Condition ' num2str(cc)]);
end
sgtitle([' Plot for Subjet ' num2str(sbj_no(sb))]);

end

close all




































%% ADVANCED PLOTTING

% Average of all epochs per channel
for sb=1:subject_count
figure;

pop_timtopo(subject.EEG{sb},[-200  1200], [NaN], 'ERP data and scalp maps' ,'maplimits',[-10 15],'electrodes','labels');

figure;
pop_plottopo(subject.EEG{sb}, [1:34] , 'ERP data by channel',0, 'ydir',1)
end

%% MOVIE
figure; [Movie,Colormap] = eegmovie(subject.data_pChan_pTime_pCond{sb}(:,400:1000,1), EEG.srate, EEG.chanlocs, 'vert', 0, 'startsec', -0.1, 'topoplotopt', {'numcontour' 0});
seemovie(Movie,-2,Colormap);

%% MOVIE HEADPLOT

headplotparams1 = { 'meshfile', 'mheadnew.mat' , 'transform', [-0.657149,-7.7344,7.702,0.0661127,0.0233321,-1.54609,92.4703,87.6027,82.8731] };



% set up the spline file
headplot('setup', EEG.chanlocs, 'STUDY_headplot.spl', headplotparams1{:}); close
figure('color', 'w'); [Movie,Colormap] = eegmovie( subject.data_pChan_pTime_pCond{sb}(:,400:900,1), EEG.srate, EEG.chanlocs_new, 'framenum', 'on', 'vert', 0, 'startsec', -0.1, 'mode', '3d', 'headplotopt', { headplotparams1{:}, 'material', 'metal'}, 'camerapath', [80 2 30 0]); 
seemovie(Movie,-5,Colormap);










