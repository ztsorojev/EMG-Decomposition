% EEM255 - Neuroengineering
% Project Phase II
%
% Author: Zaurbek Tsorojev
% ID: 805029443
%****************************************************
close all;
warning('off','all')
tic;
%**************************************************************
%% Input file
%**************************************************************
% Plot the EMG data as well as the spectrum to find the cutoff frequency for the filter.

M = csvread('EMG_example_2_fs_2k.csv'); %read in csv file
%M = csvread('EMG_example_1_90s_fs_2k.csv'); %CHANGE THE THRESHOLD (line 123) TO 5*RMS IF USING THIS DATABASE
%M = csvread('EMG_example_20s_2000Hz.csv');

time= M(:,1); % first column is the time series

L = length(time);
fs= (time(2)-time(1))^-1; % calculate the sample frequecy

    %%*************************************
    % IF YOU CHANGE DATABASE: add this code and put z = 0
    %%*************************************
    %{
    fs= 2000; % the sample frequency
    time= 0:1/fs:(length(M(:,1))-1)/fs; 
    time= time'; % transforming time into a column vector to stay consistent 
                 % with the definition from the first database
    L = length(time);
    %}
    
f = fs/2*linspace(0,1,(L/2)+1);

%%********************
z = 1; %PUT TO 0 IF USING OTHER DATABASE (where 1st column is not time)!!!
       %Represents column for time series;
%%********************

channel_number= size(M,2)-z; % num of channels in the database

for i=1:channel_number,
    
    figure('Color',[1 1 1]);plot(time,M(:,i+z)); %plot each channel
    str= sprintf('Channel %d',i);
    xlabel('seconds');title(str);xlim([time(1) time(end)]); % label and title each plots

    %Calculate the Fourier transform to find the frequency for the passband filter
    Y = fft(M(:,i+z));
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1); %Single-sided spectrum
    P1(2:end-1) = 2*P1(2:end-1);

    %Plot the spectrum
    figure('Color',[1 1 1]);
    plot(f, P1);
    str= sprintf('Single-Sided Amplitude Spectrum of Channel %d',i);
    xlabel('f (Hz)');title(str); % label and title each plots

end
channel_select= 1; % select channel for testing. channel_select<= channel_number
test_input= M(:,channel_select+z); % test_input will go through all the individual sections

%**************************************************************
%% Filter Signal
%**************************************************************
% Differentiate the signal to improve the SNR and apply a high-pass to remove the noise. 

% Noise shaping filter - 1st derivative: improves the SNR
test_input(1:end-1) = diff(test_input);

%Butterworth highpass filter to remove undesired frequencies
[b,a] = butter(15,150/(fs/2),'high');
y = filter(b,a,test_input); %Apply the filter

figure('Color',[1 1 1]);
plot(time,y);
str= sprintf('Filtered signal of Channel %d', channel_select);
title(str);
xlabel('Ttime (s)');
ylabel('Voltage');
xlim([0 time(end)])

% Plot spectrum of filtered signal to see the effect of filtering
Y = fft(y);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure('Color',[1 1 1]);
plot(f, P1);
str= sprintf('Single-Sided Amplitude Spectrum of Filtered Channel %d',channel_select);
xlabel('f (Hz)');title(str);

%**************************************************************
%% Detect Spikes
%**************************************************************
% Detect all the peaks above the threshold.


% I choose a window of 15 ms for spike alignment and set it as the minimum
% distance between two peaks in findpeaks()
W = 0.015;
window = W*fs; %number of samples in the window


% FIND THE THRESHOLD = 3*RMS(baseline noise)
%***********
% METHOD: calculate the RMS of the signal on a 0.2 second window that you
% shift until reaching half the time span of the signal (we assume it is
% enough time to be able to extract the noise). The minimum of all RMS
% values will correspond to the RMS of the baseline noise.
Window = fs/5; % 0.2 second corresponds to fs points
Rms = zeros(1,round(time(end)/2));

for i=0:length(Rms)-1
    Rms(i+1) = rms(y(i*Window+1:(i+1)*Window));
end

RMS = min(Rms);
threshold = 4*RMS*ones(1, length(time)); 


%Find the peaks of the signal above the threshold and set minimum distance
%between peaks (to avoid counting the same two times)
[pks_pos,lcts_pos] = findpeaks(y, time', 'MinPeakHeight', threshold(1),'MinPeakDistance', W);
[pks_nneg,lcts_neg] = findpeaks(-y, time', 'MinPeakHeight', threshold(1), 'MinPeakDistance', W);

pks_neg = (-1).*pks_nneg; %pks_nneg contains the oppositive values of the negative peaks 
                          %as we took -y in the findpeak() function
                          

% Concat all the peaks locations into one vector and same for positions
lcts = [lcts_pos  lcts_neg];
pks = [pks_pos'  pks_neg'];


%Plot all the selected peaks
figure('Color',[1 1 1]);
plot(time,y, time, threshold, time, -threshold); 
hold on
plot(lcts_pos,pks_pos,'o', lcts_neg, pks_neg,'o','MarkerSize',5)
str= sprintf('Spike Detection of Channel %d',channel_select);
xlabel('Time (s)'); ylabel('Voltage'); title(str);xlim([time(1) time(end)]); 
str= sprintf('Detect Spikes of Channel %d',channel_select);
title(str);
xlabel('Time (s)');
ylabel('Voltage');


%**************************************************************
%% Align Spikes
%**************************************************************
% Plot all the spikes in the same time window.


%******
%Find indexes that will be used later in the for loop for the spike
%alignment
%******
  % Calculate index of t = W
  t_temp=abs(time-W); %now the 0 will be the index of t=W so we can use the function find()
  index_W=find(t_temp==(min(t_temp)),1); %find() will find the 0 element (the new min)

  % Calculate index of t = t(end) - W 
  t_temp2=abs(time-(time(end) - W)); 
  index_Wend=find(t_temp2==(min(t_temp2)),1); 
  
  %time interval for the spikes alignment window
  t_spike = 0:1/fs:W;
  samples = t_spike*fs; %conversion of the time vector to the number of samples
  
% Preallocation of matrix that will contain all the spikes (all the positve
% first, then all the negative)
X = zeros(W*fs+1, length(pks));

% Plot the spike alignement
hold on;
figure('Color',[1 1 1]);
for i=1:length(pks)
    
    % When the spikes are more than W/2 seconds of the border 0s and
    % time(end)
    if(lcts(i)>W/2 && lcts(i)<time(end)-W/2)
        %interval = lcts(i)-W/2:1/fs:lcts(i)+W/2;
        
        %Find index of t = lcts(i)+W/2
        t_temp=abs(time - (lcts(i)+W/2));
        index_p=find(t_temp==(min(t_temp)),1);
        
        %Find index of t = lcts(i)-W/2
        t_temp2=abs(time - (lcts(i)-W/2));
        index_a=find(t_temp2==(min(t_temp2)),1); 
        
        %Add the spikes to X and plot it
        X(:,i) = y(index_a:index_p);
        plot(samples, X(:,i));
        hold on
        
    % When the peak is within W/2 of the border time(end), it will not be
    % centered in the spike alignment plot
    elseif(lcts(i)>time(end)-W/2)    
        X(:,i) = y(index_Wend:end);
        plot(samples, X(:,i));
        hold on 
        
    % When the peak is within W/2 of the beginning time, it will not be
    % centered in the spike alignment plot
    else 
        X(:,i) = y(1:index_W);
        plot(samples, X(:,i));
        hold on
    end
end
str= sprintf('Spike Alignment of Channel %d',channel_select);
xlabel('Samples');
ylabel('Voltage');
title(str);
hold off; 

%**************************************************************
%% Extract Features
%**************************************************************
% Apply PCA method to extract the features of the spikes


[coeff,score,variance,~,~,Xm] = pca(X'); %the symbol ~ is used to skip unneeded outputs
%Note: we transpose X because for pca(), rows = spikes and columns = time samples

% Plot the first two PC scores (which are the most significant components)
figure('Color',[1 1 1]);
scatter3(score(:,1),score(:,2), score(:,3),15)
title('Scatterplot of the First Two PC Scores');
xlabel('PC Score 1');
ylabel('PC Score 2');
zlabel('PC Score 3');
box on

%**************************************************************
%% Cluster Spikes
%**************************************************************
% Partition the spikes into clusters to find the number of neurons
% responsible for the signal

K = 3; %number of clusters

[idx,C] = kmeans(score,K); %Clustering:idx contains

% Plot the kmean clustering in 3D
figure('Color',[1 1 1]);
scatter3(score(:,1),score(:,2), score(:,3),15,idx)
title('Kmean Clustering');
xlabel('PC Score 1');
ylabel('PC Score 2');
zlabel('PC Score 3');
box on

% Plot the kmean clustering, with the centroid in 2D
figure('Color',[1 1 1]);
plot(score(idx==1,1),score(idx==1,2),'r.','MarkerSize',12)
hold on
plot(score(idx==2,1),score(idx==2,2),'b.','MarkerSize',12)
hold on
plot(score(idx==3,1),score(idx==3,2),'g.','MarkerSize',12)
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3)
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')

title('Kmean Clustering');
xlabel('PC Score 1');
ylabel('PC Score 2');
hold off


%**************************************************************
%% Classify Spikes
%**************************************************************
% Plot the spikes alignement where each spike is attributed to a cluster

% Split X based on spikes of each cluster
Xa = X(:,idx==1)';
Xb = X(:,idx==2)';
Xc = X(:,idx==3)';

% Plot the clustered spike alignement
hold on;
figure('Color',[1 1 1]);
for i=1:length(Xa(:,1))
    plot(samples, Xa(i,:),'r')
    hold on 
end
for i=1:length(Xb(:,1))
    plot(samples, Xb(i,:),'b')
    hold on    
end
for i=1:length(Xc(:,1))
    plot(samples, Xc(i,:),'g')
    hold on    
end
title('Spike Clustering');
xlabel('Samples');
ylabel('Voltage');
hold off;

%**************************************************************
%% Analysis
%**************************************************************
% Plot the templates of the action potentials based on the clustering
% Plot the selected spikes on the original signal based on clustering


% Preallocate the vector containing the templates
templates = zeros(K, length(coeff(:,1)));

% Create the templates using the outputs of pca and kmeans
% Xm: is the mean, C: the centroid of the cluster, coeff: the principal component
% coefficients
for i = 1:K
    templates(i,:) = Xm + C(i,1)*coeff(:,1)';
end

% Plot the templates
hold on;
figure('Color',[1 1 1]);
plot(samples, templates(1,:),'r', 'LineWidth', 2);
hold on
plot(samples, templates(2,:),'b', 'LineWidth', 2);
hold on
plot(samples, templates(3,:),'g', 'LineWidth', 2);
str= sprintf('Spike Templates of Channel %d',channel_select);
xlabel('Samples');
ylabel('Voltage');
title(str);
hold off;

%Split the peaks vector into its clusters
pks1 = pks(idx==1);
pks2 = pks(idx==2);
pks3 = pks(idx==3);


% Plot the firing pattern 
figure('Color',[1 1 1]);
subplot(3,1,1)
hold on
plot(time,y); 
plot(lcts(idx==1),pks(idx==1),'ro','MarkerSize',3)
hold off
str= sprintf('Firing Pattern of Channel %d',channel_select);
xlim([0 time(end)])
ylabel('Voltage');
title(str);
box on

subplot(3,1,2)
hold on
plot(time,y);
plot(lcts(idx==2),pks(idx==2),'bo','MarkerSize',3)
hold off
xlim([0 time(end)])
ylabel('Voltage');
box on

subplot(3,1,3)
hold on
plot(time,y);
plot(lcts(idx==3),pks(idx==3),'go','MarkerSize',3)
hold off
xlabel('Samples');
xlim([0 time(end)])
ylabel('Voltage');
box on


toc;
