%% Analyse V0.1.10
%           spike2 and Tracker (Manu) PhD SetUp Ben for Data Aron Isabella Inga
%
% /edit: INPUT & OUTPUT now by loading_XXX.m
% INPUT:    2 files that got exportet from spike2.smr to .mat UNBINNED
%           1 file from Manus Tracker (BeeTracker V1.2 (- instead of _)
%           ->values.mat
%           log.mat - if 0 new compute, if 1 take saved workspace.mat
% OUTPUT:   plots(not saved) and the workspacxeFULL.mat at the very end

%% The MIT License (MIT)
% 
% Copyright (c) 2016 Benjamin H Paffhausen
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% spikes_down are smoothed by 10 not 100 (1 frame_size instead of the usual
% 10)
%       edit: rAverageX=smooth(spikes_down,5); so rAverage of 5 Bins
% ISI multiply by 10 AND spikes_down_neu (clear/bumpy) to get diviation to 1 

%%  choose workspace
%loading_1408191255_020
%loading_1511241359_part_Ch1u2u3u4;
%loading_c1508261641_001_ALL_part_Ch1u2u3u4  % ab hälfte, ende auch weg und löcher wegen queenGroup; Ch3&4 leer?
%loading_c1508271559_001_ALL_part_Ch1u2u3u4
%load('2016-03-03_1511241359_part_Ch1u2u3u4_workspace.mat')
% 
% clear all
% close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RE-sort periferBees || inactive ||
%{
x_copy=x;
x=x_copy;

plot(x(1000:2000,5),y(1000:2000,5),'.')

x=x_copy;

marker=nan(11000,2);
for j=5 %2:12
    for i=1000:10590            %length(x(:,1))-10

        if abs(x(i+1,j)-x(i,j))>5 | isnan(x(i+1,j))
            for g=1:20
                [c, index]=min(abs(x(i+g,2:12)-x(i,j)));
                wuerfel(1,g)=c;
                wuerfel(2,g)=index+1;
            end
            breaker=1;
            for k=1:20
                if wuerfel(1,k)<3 && breaker==1
                    temp=0;
                    temp=x(i+k:end,j);
                    x(i+k:end,j)=x(i+k:end,wuerfel(2,k));
                    x(i+k:end,wuerfel(2,k))=temp;
                    breaker=k+10;
                    marker(i,1)=k;
                    marker(i,2)=j;
                    
                end
            end
        end
    end
end


    figure()
    hold on
    plot(x_copy(1000:10590,2:12)+1,'ro')
    plot(x_copy(1000:10590,5)+1,'go')
    plot(x(1000:10590,5),'.')
    plot(x(1000:10590,1),'y*')
    hold off
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% start
tic
figure_position1 = [2840 1080 1000 1080];   % oberes fenster [2840 1080 1000 1080]
figure_position2 = [2840 0 1000 1080];      % unteres fenster [2840 0 1000 1080]
figure_positiondouble = [2840 0 1000 2160];      % unteres fenster [2840 0 1000 1080]
distance_speed_threshold = 50;              % grenze ab der perifere bienen springen/ keine laufgeschw. [50]
rAverageSpeed_size = 50;                    % average' breite [50]
angle_resolution = 36;                      % polarPlot bin'ing [36]
turning_angle_smooth = 10;                  % bin'ing for direction/angle, 360-0° [10]
binning_dist=100;
binning_speed=100;
track_bins=50;
step_behavior=100;
density_plot_scaler=10;
distance_for_contact=50;
%rAverage1 = smooth(spikes1_down_neu,5);  %%############ new, more boxy, more raw
%rAverage2 = smooth(spikes2_down_neu,5);
rAverage1 = smooth(spikes1_down,5);  %%########### smoothed beforehand by 10
rAverage2 = smooth(spikes2_down,5);
rAverage3 = smooth(spikes3_down,5);
rAverage4 = smooth(spikes4_down,5);

% rAverage1 = spikes1_down(2:end)-spikes1_down(1:end-1);
% rAverage1=[rAverage1; 0];
% rAverage1=rAverage1-min(rAverage1);
% rAverage2 = spikes2_down(2:end)-spikes2_down(1:end-1);
% rAverage2=[rAverage2; 0];
% rAverage2=rAverage2-min(rAverage2);

% isi1=isi1.*spikes1_down_neu*10;
% isi2=isi2.*spikes2_down_neu*10;
% isi3=isi3.*spikes3_down_neu*10;
% isi4=isi4.*spikes4_down_neu*10;
% 
% isi1(isi1>1)=1;
% isi2(isi2>1)=1;
% isi3(isi3>1)=1;
% isi4(isi4>1)=1;

isi1_keep=isi1;
isi2_keep=isi2;
isi3_keep=isi3;
isi4_keep=isi4;

rAverage1_keep=rAverage1;
rAverage2_keep=rAverage2;
rAverage3_keep=rAverage3;
rAverage4_keep=rAverage4;
 
% isi1=(nanvar(isi1)*1.5)+nanmean(isi1);
% isi2=(nanvar(isi2)*1.5)+nanmean(isi2);

isi_grenze=20;
isi1corr=isi1.*spikes1_down_neu*10;         % ISI um die aktuelle spikerate korrigiert 
isi1corr(isi1corr>isi_grenze)= isi_grenze;
isi1corr(isi1corr==0)=NaN;
isi2corr=isi2.*spikes2_down_neu*10;
isi2corr(isi2corr>isi_grenze)= isi_grenze;
isi2corr(isi2corr==0)=NaN;
isi3corr=isi3.*spikes3_down_neu*10;
isi3corr(isi3corr>isi_grenze)= isi_grenze;
isi3corr(isi3corr==0)=NaN;
isi4corr=isi4.*spikes4_down_neu*10;
isi4corr(isi4corr>isi_grenze)= isi_grenze;
isi4corr(isi4corr==0)=NaN;

% rAverage1(rAverage1<2)=2;
% rAverage1(rAverage1>2)=10;
%rAverage1 = var1;
%rAverage1=(log10(((10+rAverage1))))-1;
%rAverage2=(log10(((10+rAverage2))))-1;
%rAverage1(rAverage1>1.5*var(rAverage1)+mean(rAverage1))=1.5*var(rAverage1)+mean(rAverage1);
%rAverage2(rAverage2>1.5*var(rAverage2)+mean(rAverage2))=1.5*var(rAverage2)+mean(rAverage2);

isi1=isi1corr;
isi2=isi2corr;
isi3=isi3corr;
isi4=isi4corr;

% 
% rAverage1=smooth(rAverage1,100);
% rAverage2=smooth(rAverage2,100);
% rAverage3=smooth(rAverage3,100);
% rAverage4=smooth(rAverage4,100);

max1=max(rAverage1);
max2=max(rAverage2);
if max2<max1                % the bigger rAverage max stays for both colorbars
    max2=max1;
end

%% spikeshape overlay Unit 1 2 3 4  [1]

figure('OuterPosition',figure_position1)
subplot(2,2,1)
plot(1)
%{
valuesHK1=values1';
plot(valuesHK1(:,1:100:end))
%xlim([0 30]);
%ylim([-0.2 .2]);
set(gca,'YDir', 'reverse')
str = {'spikes in Unit1: ',length(values1),'spikes displayed:',floor(length(values1)/100)};
annotation('textbox',[.2 .6 .3 .3],'String',str,'FitBoxToText','on');
subplot(2,2,2)
valuesHK2=values2';
plot(valuesHK2(:,1:100:end))
% xlim([0 30]);
% ylim([-0.2 .2]);
set(gca,'YDir', 'reverse')
str = {'spikes in Unit2: ',length(values2),'spikes displayed:',floor(length(values2)/100)};
annotation('textbox',[.7 .6 .3 .3],'String',str,'FitBoxToText','on');
subplot(2,2,3)
valuesHK3=values3';
plot(valuesHK3(:,1:100:end))
% xlim([0 30]);
% ylim([-0.2 .2]);
set(gca,'YDir', 'reverse')
str = {'spikes in Unit3: ',length(values3),'spikes displayed:',floor(length(values3)/100)};
annotation('textbox',[.2 0 .3 .3],'String',str,'FitBoxToText','on');
subplot(2,2,4)
valuesHK4=values4';
plot(valuesHK4(:,1:100:end))
% xlim([0 30]);
% ylim([-0.2 .2]);
set(gca,'YDir', 'reverse')
str = {'spikes in Unit4: ',length(values4),'spikes displayed:',floor(length(values4)/100)};
annotation('textbox',[.7 0 .3 .3],'String',str,'FitBoxToText','on');
%}


%% rAverage combinatoric 4x [2]
figure('OuterPosition',figure_position1)
subplot(2,2,1)
scatter(rAverage2,rAverage3,[],rAverage1,'filled')
title('rAverage2,rAverage3,[],rAverage1')
subplot(2,2,2)
scatter(rAverage3,rAverage4,[],rAverage2,'filled')
title('rAverage3,rAverage4,[],rAverage2')
subplot(2,2,3)
scatter(rAverage4,rAverage1,[],rAverage3,'filled')
title('rAverage4,rAverage1,[],rAverage3')
subplot(2,2,4)
scatter(rAverage1,rAverage2,[],rAverage4,'filled')% spikerate over time rAverage1
title('rAverage1,rAverage2,[],rAverage4')


%% rAverage 1-4 timePlot [3]
figure('OuterPosition',figure_position1)
subplot(2,2,1)
plot(rAverage1*10)
title('spikerate over time rAverage1 in Hz')
subplot(2,2,2)
plot(rAverage2*10)
title('spikerate over time rAverage2 in Hz')
subplot(2,2,3)
plot(rAverage3*10)
title('spikerate over time rAverage3 in Hz')
subplot(2,2,4)
plot(rAverage4*10)
title('spikerate over time rAverage4 in Hz')


%% rAverage combinatoric 4x rAverage1_delta [4 & 5]
rAverage1_low=smooth(rAverage1,100);
rAverage2_low=smooth(rAverage2,100);
rAverage3_low=smooth(rAverage3,100);
rAverage4_low=smooth(rAverage4,100);
rAverage1_delta=rAverage1_low-rAverage1;
rAverage2_delta=rAverage2_low-rAverage2;
rAverage3_delta=rAverage3_low-rAverage3;
rAverage4_delta=rAverage4_low-rAverage4;
figure('OuterPosition',figure_position2)
subplot(2,2,1)
scatter(rAverage2_delta,rAverage3_delta,[],rAverage1_delta,'filled')
title('rAverage2,rAverage3,[],rAverage1 HIGHPASS')
subplot(2,2,2)
scatter(rAverage3_delta,rAverage4_delta,[],rAverage2_delta,'filled')
title('rAverage3,rAverage4,[],rAverage2 HIGHPASS')
subplot(2,2,3)
scatter(rAverage4_delta,rAverage1_delta,[],rAverage3_delta,'filled')
title('rAverage4,rAverage1,[],rAverage3 HIGHPASS')
subplot(2,2,4)
scatter(rAverage1_delta,rAverage2_delta,[],rAverage4_delta,'filled')
title('rAverage1,rAverage2,[],rAverage4 HIGHPASS')
figure('OuterPosition',figure_position1)
subplot(2,2,1)
plot(rAverage1_delta)
title('rAverage1 HIGHPASS')
subplot(2,2,2)
plot(rAverage2_delta)
title('rAverage2 HIGHPASS')
subplot(2,2,3)
plot(rAverage3_delta)
title('rAverage3 HIGHPASS')
subplot(2,2,4)
plot(rAverage4_delta)
title('rAverage4 HIGHPASS')


%% plot rAverage 1-4 LOW [6 & 7]
figure('OuterPosition',figure_position1)
plot([rAverage1_low rAverage2_low rAverage3_low rAverage4_low])
figure('OuterPosition',figure_position2)
plot([rAverage1 rAverage2 rAverage3 rAverage4])


%% spikerate over time rAverage2 [8]
figure('OuterPosition',figure_position2)
%scatter(rAverage1(1:end),isi1corr(1:end),10,c(1:end),'filled')
title('spikerate over time rAverage2')


%% trajectorie peripherie bienen / recBee [9]
figure('OuterPosition',figure_position1)
title('recBee trajectory over periferBees')
xlabel('video width in pixel')
ylabel('video hight in pixel')
hold on
plot(x(:,2:12),y(:,2:12),'c.')      % einfarbig
plot(x(:,1),y(:,1))                 % trajektorie Rec Bee
xlim([0 1600]);
ylim([0 1200]);
axis equal
hold off


%% neuro falsecolor over track rAverage1 [10]
figure('OuterPosition',figure_position1)
scatter(x(:,1),y(:,1),[],rAverage1,'filled')
title('trajectory of recBee, neuronal activity false color rAverage1')
colorbar
caxis([0 max1])


%% neuro falsecolor over track rAverage2 [11]
figure('OuterPosition',figure_position2)
scatter(x(:,1),y(:,1),[],rAverage2,'filled')
title('trajectory of recBee, neuronal activity false color rAverage2')
colorbar
caxis([0 max1])


%% unit1/unit2 TIMERESOLUTION [12]
figure('OuterPosition',figure_position2)
a = 10;
c = linspace(1,60,length(rAverage1));
randius=rand([length(x) 1]);
rA1_rA2_c_rand=[rAverage1 rAverage2 c' randius];
rA1_rA2_c_rand=sortrows(rA1_rA2_c_rand,4);
scatter(rA1_rA2_c_rand(:,1),rA1_rA2_c_rand(:,2),a,rA1_rA2_c_rand(:,3),'filled')
axis equal
axis square
title('Unit2 over Unit1, time cource color coded')
xlabel('sp.Frequency Unit1 in spikes/frame')
ylabel('sp.Frequency Unit2 in spikes/frame')
colorbar


%% ISI1 & ISI2 TIMERESOLUTION [13]
figure('OuterPosition',figure_position2)
c = linspace(1,60,length(rAverage1));
randius=rand([length(x) 1]);
isi1_isi2_c_rand=[isi1 isi2 c' randius];
isi1_isi2_c_rand=sortrows(isi1_isi2_c_rand,4);
scatter(isi1_isi2_c_rand(:,1),isi1_isi2_c_rand(:,2),a,isi1_isi2_c_rand(:,3),'filled')
title('ISI2 over ISI1, time cource color coded')
xlabel('ISI in sec')
ylabel('ISI in sec')
colorbar


%% ISI1 VS ISI2 VS rAverage1 VS rAvewrage2  DIAGONAL HIST [14]
figure('OuterPosition',figure_positiondouble)
X=[isi1,isi2,isi3,isi4,rAverage1,rAverage2,rAverage3,rAverage4];
plotmatrix(X)
[S,AX,BigAx,H,HAx] = plotmatrix(X);
ylabel(AX(1),'ISI1')
ylabel(AX(2),'ISI2')
ylabel(AX(3),'ISI3')
ylabel(AX(4),'ISI4')
ylabel(AX(5),'rAverage1')
ylabel(AX(6),'rAverage2')
ylabel(AX(7),'rAverage3')
ylabel(AX(8),'rAverage4')
xlabel(AX(8),'ISI1')
xlabel(AX(16),'ISI2')
xlabel(AX(24),'ISI3')
xlabel(AX(32),'ISI4')
xlabel(AX(40),'rAverage1')
xlabel(AX(48),'rAverage2')
xlabel(AX(56),'rAverage3')
xlabel(AX(64),'rAverage4')

% rAverage5 = smooth(spikes5_down,5);
% rAverage6 = smooth(spikes6_down,5);
% rAverage7 = smooth(spikes7_down,5);
% rAverage8 = smooth(spikes8_down,5);
% rAverage9 = smooth(spikes9_down,5);
% rAverage10 = smooth(spikes10_down,5);
% rAverage11 = smooth(spikes11_down,5);
% rAverage12 = smooth(spikes12_down,5);
% X=[isi1,isi2,isi3,isi4,isi5,isi6,isi7,isi8,isi9,isi10,isi11,isi12,rAverage1,rAverage2,rAverage3,rAverage4,rAverage5,rAverage6,rAverage7,rAverage8,rAverage9,rAverage10,rAverage11,rAverage12];
% plotmatrix(X)
% [S,AX,BigAx,H,HAx] = plotmatrix(X);


%% unit1/unit2  NOW IN 3D [15]
surface_bin=200;
testrAverage1=floor((surface_bin/max(rAverage1))*(rAverage1));
testrAverage2=floor((surface_bin/max(rAverage2))*(rAverage2));
testrAverage1(testrAverage1==0)=1;
testrAverage2(testrAverage2==0)=1;
unit_hist=zeros(surface_bin, surface_bin);
for i=1:length(rAverage1)
    unit_hist(testrAverage1(i),testrAverage2(i))=unit_hist(testrAverage1(i),testrAverage2(i))+1;
end
unit_hist=sqrt(unit_hist);           %   #######################WURZEL GEZOGEN!!!!!###################
figure('OuterPosition',figure_position1)
surface(unit_hist(2:end,2:end),'edgecolor','none')
%view(3)
title('Unit2 over Unit1 - density plot  WURZEL GEZOGEN!!!!!')
xlabel('sp.Frequency rAverage1 WURZEL GEZOGEN!!!!!')
ylabel('sp.Frequency rAverage2 WURZEL GEZOGEN!!!!!')
colorbar
%set(gca, 'ZScale', 'log')


%% unit1/unit2 NOW IN 3D TIMEBLOCKER (5min) 1/2 & 2/1 [16]
figure('OuterPosition',figure_position1)
for j=1:12
subplot(2,12,j)
unit_hist=zeros(max(testrAverage1), max(testrAverage2));
for i=1+(floor(length(rAverage1)/12)*(j-1)):floor(length(rAverage1)/12)*j
    unit_hist(testrAverage1(i),testrAverage2(i))=unit_hist(testrAverage1(i),testrAverage2(i))+1;
end
unit_hist=sqrt(unit_hist);   
surface(unit_hist,'edgecolor','none')
k=.0765*(j-1)+.01;
set(gca,'Position',[k .52 .075 .45])
%set(gca, 'ZScale', 'log')
end
str = {'Unit2 over Unit1'};
annotation('textbox',[.8 .6 .3 .3],'String',str,'FitBoxToText','on','EdgeColor','none');
for j=1:12
subplot(2,12,j+12)
unit_hist=zeros(max(testrAverage1), max(testrAverage2));
for i=1+(floor(length(rAverage1)/12)*(j-1)):floor(length(rAverage1)/12)*j
    unit_hist(testrAverage1(i),testrAverage2(i))=unit_hist(testrAverage1(i),testrAverage2(i))+1;
end
unit_hist=sqrt(unit_hist);   
surface(unit_hist','edgecolor','none')
k=.0765*(j-1)+.01;
set(gca,'Position',[k .02 .075 .45])
%set(gca, 'ZScale', 'log')
end
str = {'Unit1 over Unit2'};
annotation('textbox',[.8 .1 .3 .3],'String',str,'FitBoxToText','on','EdgeColor','none');


%% unit1/unit2 spikeHistogramm  - timeResolution [17]
figure('OuterPosition',figure_position1)
rAverage1_1000=smooth(rAverage1,1000);       % rolling average 5 werte gemittelt (standart)
rAverage2_1000=smooth(rAverage2,1000);
scatterhist(rAverage1_1000,rAverage2_1000,'MarkerSize',1)
axis square
axis equal
title('Unit2 over Unit1 rAverage 2min')
xlabel('sp.Frequency # Unit1 in spikes/frame')
ylabel('sp.Frequency # Unit2 in spikes/frame')
clear rAverage1_1000;
clear rAverage1_1000;


%% unit1/unit2 spikeHistogramm     RolAverage [18]
figure('OuterPosition',figure_position2)
scatterhist(rAverage1,rAverage2,'MarkerSize',1)
axis square
axis equal
title('Unit2 over Unit1 rAverage')
xlabel('sp.Frequency # Unit1 in spikes/frame')
ylabel('sp.Frequency # Unit2 in spikes/frame')


%% isi1/isi2 spikeHistogramm     interspikeIntervall [19]
figure('OuterPosition',figure_position2)
scatterhist(isi1,isi2,'MarkerSize',1)
axis square
axis equal
title('Unit2 over Unit1 ISI')
xlabel('interspikeInterval # Unit1 in sec')
ylabel('interspikeInterval # Unit2 in sec')


%% distance/speed periferBees
% distance > 50 = sprünge - entfernen
speed=zeros(size(x));
for i=1:length(x)
    for j=1:12
        if distance(i,j)<distance_speed_threshold
            speed(i,j)=distance(i,j);
        end
    end
end
% rolling average speed zu rAverageSpeed
rAverageSpeed=zeros(size(x));
for i=1:12
    rAverageSpeed(:,i)=smooth(speed(:,i),rAverageSpeed_size);
end


%% Histogramm RecBee SPEED [20]
figure('OuterPosition',figure_position1)
hist(rAverageSpeed(:,1),1000)
title('walking speed, averaded, histo of RecBee')
xlabel('speed in pixel/frame')
ylabel('#')
xlim([0 10])


%% RecBee SPEED vs spiking [21]
rAverageSpeed_rec=rAverageSpeed(:,1);
%rAverageSpeed_rec(rAverageSpeed_rec>10)=10;
rAverageSpeed_rAverage1 = [rAverageSpeed_rec, rAverage1];
rAverageSpeed_rAverage1_sort=sortrows(rAverageSpeed_rAverage1);
speed_hist_rAverage1=zeros(floor(10*length(x)/binning_speed),binning_speed);
speed_hist_rAverage1(speed_hist_rAverage1==0)=nan;
max_speed=max(rAverageSpeed_rAverage1_sort(:,1));
step_speed=max_speed/binning_speed;
j=1;
k=1;
for i=1:length(x)
    if  rAverageSpeed_rAverage1_sort(i,1)<step_speed*j && rAverageSpeed_rAverage1_sort(i,1)>= step_speed*(j-1)
        speed_hist_rAverage1(k,j)=rAverageSpeed_rAverage1_sort(i,2);
        k=k+1;
    else
        speed_hist_rAverage1(1,j+1)=rAverageSpeed_rAverage1_sort(i,2);
        j=j+1;
        k=2;
    end
end
speed_hist_rAverage1=speed_hist_rAverage1(:,1:end-1);
figure('OuterPosition',figure_position1)
subplot(3,1,1)
speed_hist_rAverage1(speed_hist_rAverage1==0)=nan;
boxplot(speed_hist_rAverage1)
title('walking speed, averaded vs rAverage1 spike activity')
xlabel('speed from lowest to highest')
ylabel('spikes per frame')

speed_hist_rAverage1_amount=NaN(1,binning_speed);
for i=1:size(speed_hist_rAverage1,2)
    speed_hist_rAverage1_amount(i)=sum(isnan(speed_hist_rAverage1(:,i)));
end
subplot(3,1,2)
plot((size(speed_hist_rAverage1,1))-speed_hist_rAverage1_amount);
title('walking speed, ocurance')
xlabel('speed from lowest to highest')
ylabel('#')

rAverageSpeed_rAverage2 = [rAverageSpeed(:,1), rAverage2];
rAverageSpeed_rAverage2_sort=sortrows(rAverageSpeed_rAverage2);
speed_hist_rAverage2=zeros(floor(10*length(x)/binning_speed),binning_speed);
speed_hist_rAverage2(speed_hist_rAverage2==0)=nan;
j=1;
k=1;
for i=1:length(x)
    if  rAverageSpeed_rAverage2_sort(i,1)<step_speed*j && rAverageSpeed_rAverage2_sort(i,1)>= step_speed*(j-1)
        speed_hist_rAverage2(k,j)=rAverageSpeed_rAverage2_sort(i,2);
        k=k+1;
    else
        speed_hist_rAverage2(1,j+1)=rAverageSpeed_rAverage2_sort(i,2);
        j=j+1;
        k=2;
    end
end
speed_hist_rAverage2=speed_hist_rAverage2(:,1:end-1);
subplot(3,1,3)
boxplot(speed_hist_rAverage2)
title('walking speed, averaded vs rAverage2 spike activity')
xlabel('speed from lowest to highest')
ylabel('spikes per frame')
% ylim([0 10])


%% Histogramm periferBees SPEED [22]
rAverageSpeedrow=zeros(length(x)*11,1);
for j=2:12                         % all in one row for histogramm
    for i=1:size(x)
        rAverageSpeedrow(i+(size(x)*(j-2)))=rAverageSpeed(i,j);
    end
end
figure('OuterPosition',figure_position2)
hist(rAverageSpeedrow,1000)
xlim([0 10])
ylim([0 20000])
title('walking speed, averaded, histo of periferBees')
xlabel('speed in pixel/frame')
ylabel('#')


%% turning directions [23]
turn=zeros(length(x),1);
smooth_angle=smooth(angle(:,1),turning_angle_smooth);
for i=1:length(x)-60
    turn(i)=angle(i,1)-angle(i+1,1);
end
direction_hist=hist(turn);      % 2 and 9 turn oder 360; 5-decrease, 6 increase
% (left/right turn, delta = turning distro
disp('left turns:')
disp(direction_hist(5))
disp('right turns:')
disp(direction_hist(6))
disp('left turnsovers:')
disp(sum(direction_hist(1:2)))
disp('right turnsovers:')
disp(sum(direction_hist(9:10)))
% angle over time RecBee
figure('OuterPosition',figure_position1)
%str = {'left: ',direction_hist(5),'right: ',direction_hist(6),'left rollover: ',sum(direction_hist(1:2)),'right rollover: ',sum(direction_hist(9:10))};
%annotation('textbox',[.2 0 .3 .3],'String',str,'FitBoxToText','on','EdgeColor','none');
subplot(2,1,1);
plot(angle(:,1));
xlabel('time in frames')
ylabel('angle of recBee in degree')
title('angular orientation over time of RecBee')
testtext=0;
angle_full=angle(:,1);
for i=1:length(x)-1
if angle_full(i)-angle_full(i+1)>200 
    angle_full(i+1:end)=angle_full(i+1:end)+360;
    testtext=testtext+1;
elseif angle_full(i+1)-angle_full(i)>200
    angle_full(i+1:end)=angle_full(i+1:end)-360;
end
end
subplot(2,1,2);
str = {'overall turns: ',floor((angle_full(end)-angle_full(1))/360)};
annotation('textbox',[.2 0 .3 .3],'String',str,'FitBoxToText','on','EdgeColor','none');
plot(angle_full)
xlabel('time in frames')
ylabel('angle of recBee in degree')
title('angular orientation over time of RecBee cumulative (360° = 0°)')


%% Activity1 per angel [24 & 25]
figure('OuterPosition',figure_position1)
subplot(2,1,1)
angle_resolution_corr=360/angle_resolution;
angleUnits=zeros(floor((length(x)/100)*angle_resolution_corr),angle_resolution);
angleUnits(angleUnits == 0) = NaN;
for j=1:angle_resolution
    h=1;
    for i=1:length(x)
        if angle(i,1) > (j-1)*angle_resolution_corr && angle(i,1) <= j*angle_resolution_corr
            angleUnits(h,j)=rAverage1(i);
            h=h+1;
        end
    end
end
boxplot(angleUnits,'plotstyle','compact','whisker',3)
ylabel('distro of spike Activity Unit1')
xlabel('angle of recBee in degree')
title('amound of spikes per frame of Unit 1 per angle RecBee')

angleUnits_rose_mean=zeros(angle_resolution,1);
for i=1:int64(angle_resolution)
angleUnits_rose_mean(i)=nanmedian(angleUnits(:,i));
end
angleUnits_rose=zeros(max(floor(angleUnits_rose_mean*100)),1);
j=1;
for i=1:int64(angle_resolution)
    for k=1:floor(100*angleUnits_rose_mean(i))
angleUnits_rose(j)=i;
j=j+1;
    end
end
subplot(2,1,2)
rose(-0.01+deg2rad(10*angleUnits_rose),angle_resolution+1);
title('spike distribution of RecBee rAverage1')


% Activity2 per angel
figure('OuterPosition',figure_position2)
subplot(2,1,1)
angle_resolution_corr=360/angle_resolution;
angleUnits=zeros(floor((length(x)/100)*angle_resolution_corr),angle_resolution);
angleUnits(angleUnits == 0) = NaN;
for j=1:angle_resolution
    h=1;
    for i=1:length(x)
        if angle(i,1) > (j-1)*angle_resolution_corr && angle(i,1) <= j*angle_resolution_corr
            angleUnits(h,j)=rAverage2(i);
            h=h+1;
        end
    end
end
boxplot(angleUnits,'plotstyle','compact','whisker',3)
ylabel('distro of spike Activity Unit1')
xlabel('angle of recBee in degree')
title('amound of spikes per frame of Unit 2 per angle RecBee')

angleUnits_rose_mean=zeros(angle_resolution,1);
for i=1:int64(angle_resolution)
angleUnits_rose_mean(i)=nanmedian(angleUnits(:,i));
end
angleUnits_rose=zeros(max(floor(angleUnits_rose_mean*100)),1);
j=1;
for i=1:int64(angle_resolution)
    for k=1:floor(100*angleUnits_rose_mean(i))
angleUnits_rose(j)=i;
j=j+1;
    end
end
subplot(2,1,2)
rose(-0.01+deg2rad(10*angleUnits_rose),angle_resolution);
title('spike distribution of RecBee rAverage2')



% 
% 
% 
% %% angleUnits_rose Hz per angle of RecBee 
% angleUnits_rose_mean=zeros(angle_resolution,1);
% for i=1:int64(angle_resolution)
% angleUnits_rose_mean(i)=nanmean(angleUnits(:,i));
% end
% angleUnits_rose=zeros(max(floor(angleUnits_rose_mean*100)),1);
% j=1;
% for i=1:int64(angle_resolution)
%     for k=1:floor(100*angleUnits_rose_mean(i))
% angleUnits_rose(j)=i;
% j=j+1;
%     end
% end
% figure('OuterPosition',figure_position2) 
% subplot(2,1,1)
% rose(deg2rad(10*angleUnits_rose),angle_resolution);
% title('spike distribution of RecBee rAverage1')
% 
% subplot(2,1,2)
% 


%% angle distro of rec bee  PolarPlot [26]
figure('OuterPosition',figure_position2) 
subplot(2,1,2)
rad_angle=degtorad(angle(:,1));
rose(rad_angle,angle_resolution)
title('angular distribution of RecBee')
subplot(2,1,1)
hist(angle(:,1),angle_resolution)
title('angular distribution of RecBee')


%% activity1 over angle in timeColor [27 & 28]
figure('OuterPosition',figure_position1)
a = 10;
c = linspace(1,60,length(rAverage1));
scatter(rAverage1,angle(:,1),a,c,'filled')
ylim([0 360])
xlabel('spiks per frame of Unit1')     %Label the horizontal axis
ylabel('angle of recBee')         %Label the vertical axis
title('angular distr. of spikeActivity of Unit1, color-timecourse')
colorbar
clear a;

% activity2 over angle in timeColor
figure('OuterPosition',figure_position2)
a = 10;
c = linspace(1,60,length(rAverage2));
scatter(rAverage2,angle(:,1),a,c,'filled')
ylim([0 360])
xlabel('spiks per frame of Unit1')     %Label the horizontal axis
ylabel('angle of recBee')         %Label the vertical axis
title('angular distr. of spikeActivity of Unit2, color-timecourse')
colorbar
clear a;


%% distance to the closest periferBee over spike Activity [29]
%distance_to_main_sort=zeros(length(x),weite);
distance_to_main_sort=distance_to_main;
distance_to_main_sort(distance_to_main == 0) = nan;  % get rid of 0's by jumps/bee changes
for i=1:length(x)
    distance_to_main_sort(i,2:12) = sort(distance_to_main_sort(i,2:12)); % all other then the recBee sorted
end
distance_to_main_sort_rAverage = [distance_to_main_sort(:,2),rAverage1,rAverage2]; % closest distance, unit1, unit2
distance_to_main_sort_rAverage = sortrows(distance_to_main_sort_rAverage);
dist_hist_rAverage1=zeros(floor(10*length(x)/binning_dist),binning_dist);
dist_hist_rAverage1(dist_hist_rAverage1==0)=nan;
dist_hist_rAverage2=zeros(floor(10*length(x)/binning_dist),binning_dist);
dist_hist_rAverage2(dist_hist_rAverage2==0)=nan;
max_dist=max(distance_to_main_sort_rAverage(:,1));
step_dist=max_dist/binning_dist;
j=1;
k=1;
for i=1:length(x)
    if  distance_to_main_sort_rAverage(i,1)<step_dist*j && distance_to_main_sort_rAverage(i,1)>=step_dist*(j-1)
        dist_hist_rAverage1(k,j)=distance_to_main_sort_rAverage(i,2);
        dist_hist_rAverage2(k,j)=distance_to_main_sort_rAverage(i,3);
        k=k+1;
    else
        dist_hist_rAverage1(1,j+1)=distance_to_main_sort_rAverage(i,2);
        dist_hist_rAverage2(1,j+1)=distance_to_main_sort_rAverage(i,3);
        j=j+1;
        k=2;
    end
end
figure('OuterPosition',figure_positiondouble)
subplot(3,1,1)
boxplot(dist_hist_rAverage1)
ylabel('spiks per frame of Unit1')
xlabel('distance of closest periferBee')
title('distance to closest periferBee vs spikes per frame rollAverage, unit1')
%set(gca,'Position',[.53 .53 .44 .44])
dist_hist_rAverage1_amount=NaN(1,binning_dist);
for i=1:binning_dist
    dist_hist_rAverage1_amount(i)=sum(isnan(dist_hist_rAverage1(:,i)));
end
subplot(3,1,2)
plot((size(dist_hist_rAverage1,1))-dist_hist_rAverage1_amount);
ylabel('#')
xlabel('distance of closest periferBee')
title('distance to closest periferBee vs #')
subplot(3,1,3)
boxplot(dist_hist_rAverage2)
ylabel('spiks per frame of Unit1')
xlabel('distance of closest periferBee')
title('distance to closest periferBee vs spikes per frame rollAverage, unit2')
%set(gca,'Position',[.53 .03 .44 .44])


%% plotting the closest bees with neuronal activity as color on tracks [30 & 31]
combiTest=cat(3,distance_to_main,x,y,angle);
for i=1:length(x)
    temp(:,:)=combiTest(i,2:end,:);
    temp= sortrows(temp); % all other then the recBee sorted
    combiTest(i,2:end,:)=temp;
    clear temp;
end
figure('OuterPosition',figure_positiondouble)
subplot(2,1,1)
hold on


combiTest_rA=[combiTest(:,2,2),combiTest(:,2,3),rAverage1];
combiTest_rA=sortrows(combiTest_rA,-3);
scatter(combiTest_rA(:,1),combiTest_rA(:,2),[],combiTest_rA(:,3),'filled')




scatter(combiTest(:,1,2),combiTest(:,1,3),[],'k','filled')
title('trajectory of periferBee CLOSE, neuronal activity false color rAverage1')
colorbar
caxis([0 max1])
hold off
subplot(2,1,2)
hold on
%scatter(combiTest(:,3,2),combiTest(:,3,3),[],rAverage2,'filled')

combiTest_rA=[combiTest(:,2,2),combiTest(:,2,3),rAverage2];
combiTest_rA=sortrows(combiTest_rA,-3);
scatter(combiTest_rA(:,1),combiTest_rA(:,2),[],combiTest_rA(:,3),'filled')

scatter(combiTest(:,1,2),combiTest(:,1,3),[],'k','filled')
title('trajectory of periferBee CLOSE, neuronal activity false color rAverage2')
colorbar
caxis([0 max1])
hold off

figure('OuterPosition',figure_positiondouble)
subplot(2,1,1)
 hold on
% scatter(combiTest(:,7,2),combiTest(:,7,3),[],rAverage1,'filled')
% scatter(combiTest(:,6,2),combiTest(:,6,3),[],rAverage1,'filled')
% scatter(combiTest(:,5,2),combiTest(:,5,3),[],rAverage1,'filled')
% scatter(combiTest(:,4,2),combiTest(:,4,3),[],rAverage1,'filled')
% scatter(combiTest(:,3,2),combiTest(:,3,3),[],rAverage1,'filled')
% scatter(combiTest(:,2,2),combiTest(:,2,3),[],rAverage1,'filled')

combiTest_rA=[  combiTest(:,7,2),combiTest(:,7,3),combiTest(:,6,2),combiTest(:,6,3),...
                combiTest(:,5,2),combiTest(:,5,3),combiTest(:,4,2),combiTest(:,4,3),...
                combiTest(:,3,2),combiTest(:,3,3),combiTest(:,2,2),combiTest(:,2,3),rAverage1];
combiTest_rA=sortrows(combiTest_rA,-13);

scatter(combiTest_rA(:,1),combiTest_rA(:,2),[],combiTest_rA(:,13),'filled')
scatter(combiTest_rA(:,3),combiTest_rA(:,4),[],combiTest_rA(:,13),'filled')
scatter(combiTest_rA(:,5),combiTest_rA(:,6),[],combiTest_rA(:,13),'filled')
scatter(combiTest_rA(:,7),combiTest_rA(:,8),[],combiTest_rA(:,13),'filled')
scatter(combiTest_rA(:,9),combiTest_rA(:,10),[],combiTest_rA(:,13),'filled')
scatter(combiTest_rA(:,11),combiTest_rA(:,12),[],combiTest_rA(:,13),'filled')

scatter(combiTest(:,1,2),combiTest(:,1,3),[],'k','filled')
title('trajectory of periferBees, neuronal activity false color rAverage1')
colorbar
caxis([0 max1])
hold off
subplot(2,1,2)
hold on
% scatter(combiTest(:,7,2),combiTest(:,7,3),[],rAverage2,'filled')
% scatter(combiTest(:,6,2),combiTest(:,6,3),[],rAverage2,'filled')
% scatter(combiTest(:,5,2),combiTest(:,5,3),[],rAverage2,'filled')
% scatter(combiTest(:,4,2),combiTest(:,4,3),[],rAverage2,'filled')
% scatter(combiTest(:,3,2),combiTest(:,3,3),[],rAverage2,'filled')
% scatter(combiTest(:,2,2),combiTest(:,2,3),[],rAverage2,'filled')

combiTest_rA=[  combiTest(:,7,2),combiTest(:,7,3),combiTest(:,6,2),combiTest(:,6,3),...
                combiTest(:,5,2),combiTest(:,5,3),combiTest(:,4,2),combiTest(:,4,3),...
                combiTest(:,3,2),combiTest(:,3,3),combiTest(:,2,2),combiTest(:,2,3),rAverage2];
combiTest_rA=sortrows(combiTest_rA,-13);

scatter(combiTest_rA(:,1),combiTest_rA(:,2),[],combiTest_rA(:,13),'filled')
scatter(combiTest_rA(:,3),combiTest_rA(:,4),[],combiTest_rA(:,13),'filled')
scatter(combiTest_rA(:,5),combiTest_rA(:,6),[],combiTest_rA(:,13),'filled')
scatter(combiTest_rA(:,7),combiTest_rA(:,8),[],combiTest_rA(:,13),'filled')
scatter(combiTest_rA(:,9),combiTest_rA(:,10),[],combiTest_rA(:,13),'filled')
scatter(combiTest_rA(:,11),combiTest_rA(:,12),[],combiTest_rA(:,13),'filled')
scatter(combiTest(:,1,2),combiTest(:,1,3),[],'k','filled')
title('trajectory of periferBees, neuronal activity false color rAverage2')
colorbar
caxis([0 max1])
hold off


%% plotting the closest bees with neuronal activity as color on tracks relative to recBee rA1 & 2 [32 & 33]
% angle and position of periferBee RELATIV to egocentric recBee
x_rel=NaN(size(x));
y_rel=NaN(size(x));
angle_rel=NaN(size(x));
for i=1:length(x)
    if isnan(angle(i,1))
    else
    x_rel(i,:)=x(i,:)-x(i,1);
    y_rel(i,:)=y(i,:)-y(i,1);
    angle_rel(i,:)=mod(angle(i,:)-angle(i,1),360);
    R=rotx(angle(i,1));
    for j=2:12
        temp_a=[1;x_rel(i,j);y_rel(i,j)];
        temp_b=R*temp_a;
        y_rel(i,j)=temp_b(3);
        x_rel(i,j)=temp_b(2);
    end
    end
end

combiTest_rel=cat(3,distance_to_main,x_rel,y_rel);
for i=1:length(x)
    temp(:,:)=combiTest_rel(i,2:12,:);
    temp= sortrows(temp); % all other then the recBee sorted
    combiTest_rel(i,2:12,:)=temp;
    clear temp;
end
%  counter=zeros(length(x),1);            % for scatter3, timeresolved
%  for i=1:length(x)
%      counter(i)=i;
%  end
figure('OuterPosition',figure_positiondouble)
combiTest_rA1=[ combiTest_rel(:,12,2),combiTest_rel(:,12,3),...
                combiTest_rel(:,11,2),combiTest_rel(:,11,3),combiTest_rel(:,10,2),combiTest_rel(:,10,3),...
                combiTest_rel(:,9,2),combiTest_rel(:,9,3),combiTest_rel(:,8,2),combiTest_rel(:,8,3),...
                combiTest_rel(:,7,2),combiTest_rel(:,7,3),combiTest_rel(:,6,2),combiTest_rel(:,6,3),...
                combiTest_rel(:,5,2),combiTest_rel(:,5,3),combiTest_rel(:,4,2),combiTest_rel(:,4,3),...
                combiTest_rel(:,3,2),combiTest_rel(:,3,3),combiTest_rel(:,2,2),combiTest_rel(:,2,3),rAverage1];
combiTest_rA1=sortrows(combiTest_rA1,23);
combiTest_rA2=[ combiTest_rel(:,12,2),combiTest_rel(:,12,3),...
                combiTest_rel(:,11,2),combiTest_rel(:,11,3),combiTest_rel(:,10,2),combiTest_rel(:,10,3),...
                combiTest_rel(:,9,2),combiTest_rel(:,9,3),combiTest_rel(:,8,2),combiTest_rel(:,8,3),...
                combiTest_rel(:,7,2),combiTest_rel(:,7,3),combiTest_rel(:,6,2),combiTest_rel(:,6,3),...
                combiTest_rel(:,5,2),combiTest_rel(:,5,3),combiTest_rel(:,4,2),combiTest_rel(:,4,3),...
                combiTest_rel(:,3,2),combiTest_rel(:,3,3),combiTest_rel(:,2,2),combiTest_rel(:,2,3),rAverage2];
combiTest_rA2=sortrows(combiTest_rA2,23);

% for i=1:11 
% subplot(4,6,i)
% scatter(combiTest_rA1(:,(i*2)-1),combiTest_rA1(:,i*2),[],combiTest_rA1(:,23),'filled')
% %title('trajectory of periferBees, neuronal activity false color rAverage1')
% subplot(4,6,i+12)
% scatter(combiTest_rA2(:,(i*2)-1),combiTest_rA2(:,i*2),[],combiTest_rA2(:,23),'filled')
% end
% title('trajectory of periferBees, neuronal activity false color rAverage1 upper rows, lower rA2')
subplot(2,1,1)
scatter(combiTest_rA1(:,21),combiTest_rA1(:,22),[],combiTest_rA1(:,23),'filled')
 %title('trajectory of periferBees, neuronal activity false color rAverage1')
subplot(2,1,2)
scatter(combiTest_rA2(:,21),combiTest_rA2(:,22),[],combiTest_rA2(:,23),'filled')
title('trajectory of periferBees, neuronal activity false color ISI1CORR upper rows, lower ISI2CORR')

figure('OuterPosition',figure_positiondouble)
combiTest_rA1=[ combiTest_rel(:,12,2),combiTest_rel(:,12,3),...
                combiTest_rel(:,11,2),combiTest_rel(:,11,3),combiTest_rel(:,10,2),combiTest_rel(:,10,3),...
                combiTest_rel(:,9,2),combiTest_rel(:,9,3),combiTest_rel(:,8,2),combiTest_rel(:,8,3),...
                combiTest_rel(:,7,2),combiTest_rel(:,7,3),combiTest_rel(:,6,2),combiTest_rel(:,6,3),...
                combiTest_rel(:,5,2),combiTest_rel(:,5,3),combiTest_rel(:,4,2),combiTest_rel(:,4,3),...
                combiTest_rel(:,3,2),combiTest_rel(:,3,3),combiTest_rel(:,2,2),combiTest_rel(:,2,3),isi1corr];
combiTest_rA1=sortrows(combiTest_rA1,23);
combiTest_rA2=[ combiTest_rel(:,12,2),combiTest_rel(:,12,3),...
                combiTest_rel(:,11,2),combiTest_rel(:,11,3),combiTest_rel(:,10,2),combiTest_rel(:,10,3),...
                combiTest_rel(:,9,2),combiTest_rel(:,9,3),combiTest_rel(:,8,2),combiTest_rel(:,8,3),...
                combiTest_rel(:,7,2),combiTest_rel(:,7,3),combiTest_rel(:,6,2),combiTest_rel(:,6,3),...
                combiTest_rel(:,5,2),combiTest_rel(:,5,3),combiTest_rel(:,4,2),combiTest_rel(:,4,3),...
                combiTest_rel(:,3,2),combiTest_rel(:,3,3),combiTest_rel(:,2,2),combiTest_rel(:,2,3),isi2corr];
combiTest_rA2=sortrows(combiTest_rA2,23);

% for i=1:11 
% subplot(4,6,i)
% scatter(combiTest_rA1(:,(i*2)-1),combiTest_rA1(:,i*2),[],combiTest_rA1(:,23),'filled')
% %title('trajectory of periferBees, neuronal activity false color rAverage1')
% subplot(4,6,i+12)
% scatter(combiTest_rA2(:,(i*2)-1),combiTest_rA2(:,i*2),[],combiTest_rA2(:,23),'filled')
% end
subplot(2,1,1)
scatter(combiTest_rA1(:,21),combiTest_rA1(:,22),[],combiTest_rA1(:,23),'filled')
 %title('trajectory of periferBees, neuronal activity false color rAverage1')
subplot(2,1,2)

scatter(combiTest_rA2(:,21),combiTest_rA2(:,22),[],combiTest_rA2(:,23),'filled')
title('trajectory of periferBees, neuronal activity false color ISI1CORR upper rows, lower ISI2CORR')


%% plotting the closest bees with neuronal activity as color on tracks relative to recBee rA 3 & 4 [34 & 35]
figure('OuterPosition',figure_positiondouble)
combiTest_rA1=[ combiTest_rel(:,12,2),combiTest_rel(:,12,3),...
                combiTest_rel(:,11,2),combiTest_rel(:,11,3),combiTest_rel(:,10,2),combiTest_rel(:,10,3),...
                combiTest_rel(:,9,2),combiTest_rel(:,9,3),combiTest_rel(:,8,2),combiTest_rel(:,8,3),...
                combiTest_rel(:,7,2),combiTest_rel(:,7,3),combiTest_rel(:,6,2),combiTest_rel(:,6,3),...
                combiTest_rel(:,5,2),combiTest_rel(:,5,3),combiTest_rel(:,4,2),combiTest_rel(:,4,3),...
                combiTest_rel(:,3,2),combiTest_rel(:,3,3),combiTest_rel(:,2,2),combiTest_rel(:,2,3),rAverage3];
combiTest_rA1=sortrows(combiTest_rA1,23);
combiTest_rA2=[ combiTest_rel(:,12,2),combiTest_rel(:,12,3),...
                combiTest_rel(:,11,2),combiTest_rel(:,11,3),combiTest_rel(:,10,2),combiTest_rel(:,10,3),...
                combiTest_rel(:,9,2),combiTest_rel(:,9,3),combiTest_rel(:,8,2),combiTest_rel(:,8,3),...
                combiTest_rel(:,7,2),combiTest_rel(:,7,3),combiTest_rel(:,6,2),combiTest_rel(:,6,3),...
                combiTest_rel(:,5,2),combiTest_rel(:,5,3),combiTest_rel(:,4,2),combiTest_rel(:,4,3),...
                combiTest_rel(:,3,2),combiTest_rel(:,3,3),combiTest_rel(:,2,2),combiTest_rel(:,2,3),rAverage4];
combiTest_rA2=sortrows(combiTest_rA2,23);

subplot(2,1,1)
scatter(combiTest_rA1(:,21),combiTest_rA1(:,22),[],combiTest_rA1(:,23),'filled')
subplot(2,1,2)
scatter(combiTest_rA2(:,21),combiTest_rA2(:,22),[],combiTest_rA2(:,23),'filled')
title('trajectory of periferBees, neuronal activity false color rAverage3 upper rows, lower rAverage4')

figure('OuterPosition',figure_positiondouble)
combiTest_rA1=[ combiTest_rel(:,12,2),combiTest_rel(:,12,3),...
                combiTest_rel(:,11,2),combiTest_rel(:,11,3),combiTest_rel(:,10,2),combiTest_rel(:,10,3),...
                combiTest_rel(:,9,2),combiTest_rel(:,9,3),combiTest_rel(:,8,2),combiTest_rel(:,8,3),...
                combiTest_rel(:,7,2),combiTest_rel(:,7,3),combiTest_rel(:,6,2),combiTest_rel(:,6,3),...
                combiTest_rel(:,5,2),combiTest_rel(:,5,3),combiTest_rel(:,4,2),combiTest_rel(:,4,3),...
                combiTest_rel(:,3,2),combiTest_rel(:,3,3),combiTest_rel(:,2,2),combiTest_rel(:,2,3),isi3corr];
combiTest_rA1=sortrows(combiTest_rA1,23);
combiTest_rA2=[ combiTest_rel(:,12,2),combiTest_rel(:,12,3),...
                combiTest_rel(:,11,2),combiTest_rel(:,11,3),combiTest_rel(:,10,2),combiTest_rel(:,10,3),...
                combiTest_rel(:,9,2),combiTest_rel(:,9,3),combiTest_rel(:,8,2),combiTest_rel(:,8,3),...
                combiTest_rel(:,7,2),combiTest_rel(:,7,3),combiTest_rel(:,6,2),combiTest_rel(:,6,3),...
                combiTest_rel(:,5,2),combiTest_rel(:,5,3),combiTest_rel(:,4,2),combiTest_rel(:,4,3),...
                combiTest_rel(:,3,2),combiTest_rel(:,3,3),combiTest_rel(:,2,2),combiTest_rel(:,2,3),isi4corr];
combiTest_rA2=sortrows(combiTest_rA2,23);

subplot(2,1,1)
scatter(combiTest_rA1(:,21),combiTest_rA1(:,22),[],combiTest_rA1(:,23),'filled')
subplot(2,1,2)
scatter(combiTest_rA2(:,21),combiTest_rA2(:,22),[],combiTest_rA2(:,23),'filled')
title('trajectory of periferBees, neuronal activity false color ISI3CORR upper rows, lower ISI4CORR')


% scatter(combiTest_rA(:,3),combiTest_rA(:,4),[],combiTest_rA(:,13),'filled')
% scatter(combiTest_rA(:,5),combiTest_rA(:,6),[],combiTest_rA(:,13),'filled')
% scatter(combiTest_rA(:,7),combiTest_rA(:,8),[],combiTest_rA(:,13),'filled')
% scatter(combiTest_rA(:,9),combiTest_rA(:,10),[],combiTest_rA(:,13),'filled')
% scatter(combiTest_rA(:,11),combiTest_rA(:,12),[],combiTest_rA(:,13),'filled')
% hold on
% scatter(combiTest_rel(:,5,2),combiTest_rel(:,5,3),[],rAverage2,'filled')
% scatter(combiTest_rel(:,4,2),combiTest_rel(:,4,3),[],rAverage2,'filled')
% scatter(combiTest_rel(:,3,2),combiTest_rel(:,3,3),[],rAverage2,'filled')
% scatter(combiTest_rel(:,2,2),combiTest_rel(:,2,3),[],rAverage2,'filled')
% hold off
% title('trajectory of periferBees, neuronal activity false color rAverage2')
% end
%{
figure('OuterPosition',figure_positiondouble)
subplot(2,1,1)
hold on
scatter(combiTest_rel(:,2,2),combiTest_rel(:,2,3),[],rAverage1,'filled')
hold off
title('trajectory of periferBee CLOSEST, neuronal activity false color rAverage1')
subplot(2,1,2)
hold on
scatter(combiTest_rel(:,2,2),combiTest_rel(:,2,3),[],rAverage2,'filled')
hold off
title('trajectory of periferBee CLOSEST, neuronal activity false color rAverage2')
%}


%% distance to closest vs activity boxplots [36]
figure('OuterPosition',figure_positiondouble)
distance_rAverage1234=[distance_to_main(:,2:12) rAverage1 rAverage2 rAverage3 rAverage4];
teiler=10;
for k=1:4
distance_bool=nan(length(x),teiler);
for i = 1:length(x)
distance_rAverage1234(i,1:11)=sortrows(distance_rAverage1234(i,1:11)')';
end
distance_rAverage1234=sortrows(distance_rAverage1234,-(11+k));
maxNOW=max(distance_rAverage1234(:,11+k));
for j=1:teiler
for i=1:length(x)
    if distance_rAverage1234(i,11+k)>(j-1)*(maxNOW/teiler) && distance_rAverage1234(i,11+k)<(j)*(maxNOW/teiler)
distance_bool(i,j) =  distance_rAverage1234(i,1);
    end
end
end
subplot(2,2,k)
boxplot(distance_bool,'Orientation','horizontal')
title('distance of closest PeriferBee (x-axis) vs. activity (y-axis)')
end


%% PSTH closest contact  - peak delta > 10pxls && closness < 50pxls [37 & 38]
[closest_peak,closest_peak_i] = findpeaks(-1*smooth(distance_to_main_sort(:,2),50)); % find peaks, negativ so close = max
%[closest_peak,closest_peak_i] = findpeaks(-1*(distance_to_main_sort(:,2))); % find peaks, negativ so close = max
contact=zeros(length(x),1);
PSTH_contact_unit1=zeros(floor(length(closest_peak_i)/10),50);
PSTH_contact_unit1(PSTH_contact_unit1 == 0) = nan;  % get rid of 0's
PSTH_contact_unit2=PSTH_contact_unit1;
PSTH_contact_rAverage1=PSTH_contact_unit1;
PSTH_contact_rAverage2=PSTH_contact_unit1;
j=1;
for i=5:length(closest_peak_i)-10
      % old if conditionclosest_peak(i) - closest_peak(i-1) > 10 && 
    if    (closest_peak(i)+distance_for_contact)>0 % if a peak is closer then the last one it is a real max
        contact(closest_peak_i(i)) = closest_peak(i);  % all 0's but max
        PSTH_contact_rAverage1(j,:)=rAverage1(((closest_peak_i(i)-25):(closest_peak_i(i)+24)));
        PSTH_contact_rAverage1(j,:)=PSTH_contact_rAverage1(j,:)-PSTH_contact_rAverage1(j,25);
        PSTH_contact_rAverage2(j,:)=rAverage2(((closest_peak_i(i)-25):(closest_peak_i(i)+24)));
        PSTH_contact_rAverage2(j,:)=PSTH_contact_rAverage2(j,:)-PSTH_contact_rAverage2(j,25);
        j=j+1;
    end
end
figure('OuterPosition',figure_positiondouble) % plot closness over time, cycle contact(isch)
hold on
plot(-1*(distance_to_main_sort(:,2)))
plot(contact,'o','MarkerSize',12)
ylim([-50 0])
title('plot closness over time, cycle contact(isch)')
hold off
figure('OuterPosition',figure_positiondouble)
subplot(4,1,1)
boxplot(PSTH_contact_rAverage1)                           % PSTH of rAverage1 @ closeness
str = {'contacts: ',size(PSTH_contact_rAverage1,2)};
annotation('textbox',[.5 .3 .2 .2],'String',str,'FitBoxToText','on','EdgeColor','none');
title('PSTH of rAverage1 @ closeness +/- 25 bin')
%set(gca,'Position',[.53 .53 .44 .44])
subplot(4,1,2)
boxplot(PSTH_contact_rAverage2)                           % PSTH of rAverage2 @ closeness
title('PSTH of rAverage2 @ closeness +/- 25 bin')
%set(gca,'Position',[.53 .03 .44 .44])
subplot(4,1,3)
plot(PSTH_contact_rAverage1')
title('rAverage1 normalized @ closeness +/- 25 bin')
subplot(4,1,4)
plot(PSTH_contact_rAverage2')
title('rAverage2 normalized @ closeness +/- 25 bin')


%% neuro activity before and after contact [39]
figure('OuterPosition',figure_positiondouble)
j=0;
PSTH_contact_rAverage1_bin=zeros(size(PSTH_contact_rAverage1,1),2);
for i=1:size(PSTH_contact_rAverage1,1)
    PSTH_contact_rAverage1_bin(i,1)=0.05*sum(PSTH_contact_rAverage1(i,1:20));
    PSTH_contact_rAverage1_bin(i,2)=0.05*sum(PSTH_contact_rAverage1(i,31:50));
end
subplot(2,1,2)
boxplot(PSTH_contact_rAverage1_bin)
title('20 bins befor and 20 bins after contact summed boxplot rAverage1')

for i=1:size(PSTH_contact_rAverage1,2)   
    if PSTH_contact_rAverage1_bin(i,1)<PSTH_contact_rAverage1_bin(i,2)
        PSTH_contact_rAverage1_bin(i,:)=PSTH_contact_rAverage1_bin(i,:)+1;
        j=j+1;
    end
end
subplot(2,1,1)
plot(PSTH_contact_rAverage1_bin')
str = {'up: ',j,'down: ',length(PSTH_contact_rAverage1_bin(~isnan(PSTH_contact_rAverage1_bin(:,1))))-j};
annotation('textbox',[.5 .4 .3 .3],'String',str,'FitBoxToText','on','EdgeColor','none');
title('20 bins befor and 20 bins after contact summed, sortet & counted ascending/descending rAverage1')


%% behavior analysis, contact+-step 
j=0;
for i=1:length(x)
    if contact(i) ~= 0
        j=j+1;
    end
end
contact_counter=j;
behavior_contact=zeros(contact_counter,step_behavior,10);
behavior_contact_norm=behavior_contact;
j=1;
for i=1+(step_behavior/2):length(x)
    if contact(i) ~= 0
        behavior_contact(j,:,1)=x(i-(step_behavior/2):i+(step_behavior/2)-1,1);
        behavior_contact(j,:,2)=y(i-(step_behavior/2):i+(step_behavior/2)-1,1);
        behavior_contact(j,:,3)=angle(i-(step_behavior/2):i+(step_behavior/2)-1,1);
        behavior_contact(j,:,4)=rAverageSpeed(i-(step_behavior/2):i+(step_behavior/2)-1,1);
        behavior_contact(j,:,5)=x(i-(step_behavior/2):i+(step_behavior/2)-1,2);
        behavior_contact(j,:,6)=y(i-(step_behavior/2):i+(step_behavior/2)-1,2);
        behavior_contact(j,:,7)=angle(i-(step_behavior/2):i+(step_behavior/2)-1,2);
        behavior_contact(j,:,8)=rAverageSpeed(i-(step_behavior/2):i+(step_behavior/2)-1,2);
        behavior_contact(j,:,9)=rAverage1(i-(step_behavior/2):i+(step_behavior/2)-1);
        behavior_contact(j,:,10)=rAverage2(i-(step_behavior/2):i+(step_behavior/2)-1);
        j=j+1;
    end
end
for i=1:contact_counter
    for j=1:10
        behavior_contact_norm(i,:,j)=behavior_contact(i,:,j)-behavior_contact(i,(step_behavior/2),j);
    end
end

behavior_contact_mean=zeros(length(behavior_contact(:,1,1)),2,length(behavior_contact(1,1,:)));
behavior_contact_mean_norm=behavior_contact_mean;
for i=1:length(behavior_contact(:,1,1))
    for j=1:10
        behavior_contact_mean(i,1,j)=mean(behavior_contact(i,1:(step_behavior/2),j));
        behavior_contact_mean(i,2,j)=mean(behavior_contact(i,(step_behavior/2)-1:step_behavior,j));
    end
end
for i=1:length(behavior_contact(:,1,1))
    for j=1:10
        behavior_contact_mean_norm(i,1,j)=mean(behavior_contact_norm(i,1:(step_behavior/2),j));
        behavior_contact_mean_norm(i,2,j)=mean(behavior_contact_norm(i,(step_behavior/2)-1:step_behavior,j));
    end
end


%% behavior analysis, contact+-step [40]
figure('OuterPosition',figure_positiondouble)
str = {'each contact(close) +- step, x y angle speed REC, x y angle speed perifer colse, rAverage1 & rAverage2'};
annotation('textbox',[.1 0 .9 .99],'String',str,'FitBoxToText','on','EdgeColor','none');
j=1;
for i=1:10
    subplot(5,4,j)
    plot(behavior_contact(:,:,i)')
    subplot(5,4,j+1)
    boxplot(behavior_contact_mean(:,:,i))
    j=j+2;
end


%% behavior analysis, contact+-step  NORMALIZED [41]
figure('OuterPosition',figure_positiondouble)
str = {'each contact(close) +- step, x y angle speed REC, x y angle speed perifer colse, rAverage1 & rAverage2 NORM'};
annotation('textbox',[.1 0 .9 .99],'String',str,'FitBoxToText','on','EdgeColor','none');
j=1;
for i=1:10
    subplot(5,4,j)
    plot(behavior_contact_norm(:,:,i)')
    subplot(5,4,j+1)
    boxplot(behavior_contact_mean_norm(:,:,i))
    j=j+2;
end


%% behavior analysis, contact+-step  NORMALIZED UP & DOWN rAverage1 [42]
behavior_contact_mean_norm_rauf1=zeros(contact_counter,step_behavior,10);
behavior_contact_mean_norm_rauf1(behavior_contact_mean_norm_rauf1==0)=nan;
behavior_contact_mean_norm_runter1=zeros(contact_counter,step_behavior,10);
behavior_contact_mean_norm_runter1(behavior_contact_mean_norm_runter1==0)=nan;
j=1;
k=1;
for i=1:length(behavior_contact(:,1,1))
    if behavior_contact_mean_norm(i,1,9) < behavior_contact_mean_norm(i,2,9)
    behavior_contact_mean_norm_rauf1(j,:,:)=behavior_contact_norm(i,:,:);
    j=j+1;
    else
    behavior_contact_mean_norm_runter1(k,:,:)=behavior_contact_norm(i,:,:);
    k=k+1;
    end
end
figure('OuterPosition',figure_positiondouble)
str = {'each contact(close) +- step, x y angle speed REC, x y angle speed perifer colse, rAverage1 & rAverage2 sortet to rAverage1 left down right up'};
annotation('textbox',[.1 0 .9 .99],'String',str,'FitBoxToText','on','EdgeColor','none');
rauf_counter1=length(behavior_contact_mean_norm_rauf1(~isnan(behavior_contact_mean_norm_rauf1(:,1,10))));
runter_counter1=length(behavior_contact_mean_norm_runter1(~isnan(behavior_contact_mean_norm_runter1(:,1,10))));
str = {'rAverage1 up:',rauf_counter1,'down:',runter_counter1};
annotation('textbox',[.1 .05 .9 .92],'String',str,'FitBoxToText','on','EdgeColor','none');
j=1;
for i=1:10
    subplot(5,4,j)
    plot(behavior_contact_mean_norm_rauf1(:,:,i)')
     subplot(5,4,j+1)
        plot(behavior_contact_mean_norm_runter1(:,:,i)')
    j=j+2;
end


%% behavior analysis, contact+-step  NORMALIZED UP & DOWN rAverage2 [43]
behavior_contact_mean_norm_rauf2=zeros(contact_counter,step_behavior,10);
behavior_contact_mean_norm_rauf2(behavior_contact_mean_norm_rauf2==0)=nan;
behavior_contact_mean_norm_runter2=zeros(contact_counter,step_behavior,10);
behavior_contact_mean_norm_runter2(behavior_contact_mean_norm_runter2==0)=nan;
j=1;
k=1;
for i=1:length(behavior_contact(:,1,1))
    if behavior_contact_mean_norm(i,1,10) < behavior_contact_mean_norm(i,2,10)
    behavior_contact_mean_norm_rauf2(j,:,:)=behavior_contact_norm(i,:,:);
    j=j+1;
    else
    behavior_contact_mean_norm_runter2(k,:,:)=behavior_contact_norm(i,:,:);
    k=k+1;
    end
end
figure('OuterPosition',figure_positiondouble)
str = {'each contact(close) +- step, x y angle speed REC, x y angle speed perifer colse, rAverage1 & rAverage2 sortet to rAverage2 left down right up'};
annotation('textbox',[.1 0 .9 .99],'String',str,'FitBoxToText','on','EdgeColor','none');
rauf_counter2=length(behavior_contact_mean_norm_rauf2(~isnan(behavior_contact_mean_norm_rauf2(:,1,10))));
runter_counter2=length(behavior_contact_mean_norm_runter2(~isnan(behavior_contact_mean_norm_runter2(:,1,10))));
str = {'rAverage2 up:',rauf_counter2,'down:',runter_counter2};
annotation('textbox',[.1 .05 .9 .92],'String',str,'FitBoxToText','on','EdgeColor','none');
j=1;
for i=1:10
    subplot(5,4,j)
    plot(behavior_contact_mean_norm_rauf2(:,:,i)')
     subplot(5,4,j+1)
        plot(behavior_contact_mean_norm_runter2(:,:,i)')
    j=j+2;
end
%min(property)


%% histo rAverage1 vs distance closest - density plot         [44]         distance perifer rA1
figure('OuterPosition',figure_position1)
property=((surface_bin/max(distance_to_main_sort(:,2)))*distance_to_main(:,2));
fullrAverage1=floor((surface_bin/max(rAverage1))*(rAverage1));
 property_step=0.5;
% fullrAverage1=ceil(density_plot_scaler*(rAverage1));
 fullrAverage1(fullrAverage1==0)=1;
 fullproperty=ceil(property_step*(property));
 fullproperty(fullproperty==0)=1;
 fullproperty(isnan(fullproperty))=1;
% fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(rAverage1)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,1)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('rAverage1 vs distance closest - density plot')
xlabel('sp.Frequency rAverage1')
ylabel('property')
set(gca,'Position',[.02 .51 .19 .44])
% surface_bin=200;
% testrAverage1=floor((surface_bin/max(rAverage1))*(rAverage1));
% testrAverage2=floor((surface_bin/max(rAverage2))*(rAverage2));
% testrAverage1(testrAverage1==0)=1;
% testrAverage2(testrAverage2==0)=1;
% unit_hist=zeros(surface_bin, surface_bin);


%% histo rAverage2 vs distance closest - density plot                      distance perifer rA2
fullrAverage1=ceil(density_plot_scaler*(rAverage2));
fullrAverage1(fullrAverage1==0)=1;
fullproperty=ceil(property_step*(property));
fullproperty(fullproperty==0)=1;
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(rAverage1)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,6)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('rAverage2 vs distance closest - density plot')
xlabel('sp.Frequency rAverage1')
ylabel('property')
set(gca,'Position',[.02 .02 .19 .44])


%% histo rAverage1 vs x  - density plot                                    x coordinate recBee rA1
property=x(:,1);
property_step=0.5;
fullrAverage1=ceil(density_plot_scaler*(rAverage1));
fullrAverage1(fullrAverage1==0)=1;
fullproperty=ceil(property_step*(property));
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(rAverage1)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,2)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('rAverage1 vs x coordinat - density plot')
xlabel('sp.Frequency rAverage1')
ylabel('property')
set(gca,'Position',[.22 .51 .19 .44])


%% histo rAverage2 vs x  - density plot                                    x coordinate recBee rA2
fullrAverage1=ceil(density_plot_scaler*(rAverage2));
fullrAverage1(fullrAverage1==0)=1;
fullproperty=ceil(property_step*(property));
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(rAverage1)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,7)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('rAverage2 vs x coordinat - density plot')
xlabel('sp.Frequency rAverage2')
ylabel('property')
set(gca,'Position',[.22 .02 .19 .44])


%% histo rAverage1 vs y - density plot                                     y coordinate recBee rA1
property=y(:,1);
property_step=0.5;
fullrAverage1=ceil(density_plot_scaler*(rAverage1));
fullrAverage1(fullrAverage1==0)=1;
fullproperty=ceil(property_step*(property));
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(rAverage1)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,3)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('rAverage1 vs y coordinat - density plot')
xlabel('sp.Frequency rAverage1')
ylabel('property')
set(gca,'Position',[.42 .51 .19 .44])


%% histo rAverage2 vs y - density plot                                     y coordinate recBee rA2
fullrAverage1=ceil(density_plot_scaler*(rAverage2));
fullrAverage1(fullrAverage1==0)=1;
fullproperty=ceil(property_step*(property));
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(rAverage1)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,8)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('rAverage2 vs y coordinat - density plot')
xlabel('sp.Frequency rAverage2')
ylabel('property')
set(gca,'Position',[.42 .02 .19 .44])


%% histo rAverage1 vs angle - density plot                                 angle recBee rA1
property=angle(:,1);
property_step=0.5;
fullrAverage1=ceil(density_plot_scaler*(rAverage1));
fullrAverage1(fullrAverage1==0)=1;
fullproperty=ceil(property_step*(property));
fullproperty(fullproperty<0)=0;
fullproperty=fullproperty+1;
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(rAverage1)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,4)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('rAverage1 vs angle recBee - density plot')
xlabel('sp.Frequency rAverage1')
ylabel('property')
set(gca,'Position',[.62 .51 .19 .44])


%% histo rAverage2 vs angle - density plot                                 angle recBee rA2
fullrAverage1=ceil(density_plot_scaler*(rAverage2));
fullrAverage1(fullrAverage1==0)=1;
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(rAverage1)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,9)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('rAverage2 vs angle RecBee - density plot')
xlabel('sp.Frequency rAverage1')
ylabel('property')
set(gca,'Position',[.62 .02 .19 .44])


%% histo rAverage1 vs distance closest - density plot           [45]       distance perifer isi1
figure('OuterPosition',figure_position2)
property=distance_to_main_sort(:,2);
property_step=0.5;
density_plot_scaler=50;
fullrAverage1=ceil(density_plot_scaler*(isi1));
fullrAverage1(fullrAverage1==0)=1;
%fullrAverage1(fullrAverage1>100)=100;

fullproperty=ceil(property_step*(property));
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(isi1)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,1)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('isi1 vs distance closest - density plot')
%xlabel('isi in sec')
ylabel('property')
set(gca,'Position',[.02 .51 .19 .44])


%% histo rAverage2 vs distance closest - density plot                      distance perifer isi2
fullrAverage1=ceil(density_plot_scaler*(isi2));
fullrAverage1(fullrAverage1==0)=1;
fullproperty=ceil(property_step*(property));
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(isi2)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,6)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('isi2 vs distance closest - density plot')
%xlabel('sp.Frequency rAverage1')
ylabel('property')
set(gca,'Position',[.02 .02 .19 .44])


%% histo rAverage1 vs x  - density plot                                    x coordinate recBee isi1
property=x(:,1);
property_step=0.5;
fullrAverage1=ceil(density_plot_scaler*(isi1));
fullrAverage1(fullrAverage1==0)=1;
fullproperty=ceil(property_step*(property));
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(isi1)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,2)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('isi1 vs x coordinat - density plot')
%xlabel('sp.Frequency rAverage1')
ylabel('property')
set(gca,'Position',[.22 .51 .19 .44])


%% histo rAverage2 vs x  - density plot                                    x coordinate recBee isi2
fullrAverage1=ceil(density_plot_scaler*(isi2));
fullrAverage1(fullrAverage1==0)=1;
fullproperty=ceil(property_step*(property));
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(isi2)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,7)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('isi2 vs x coordinat - density plot')
%xlabel('sp.Frequency rAverage1')
ylabel('property')
set(gca,'Position',[.22 .02 .19 .44])


%% histo rAverage1 vs y - density plot                                     y coordinate recBee isi1
property=y(:,1);
property_step=0.5;
fullrAverage1=ceil(density_plot_scaler*(isi1));
fullrAverage1(fullrAverage1==0)=1;
fullproperty=ceil(property_step*(property));
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(isi1)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,3)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('isi1 vs y coordinat - density plot')
%xlabel('sp.Frequency rAverage1')
ylabel('property')
set(gca,'Position',[.42 .51 .19 .44])


%% histo rAverage2 vs y - density plot                                     y coordinate recBee isi2
fullrAverage1=ceil(density_plot_scaler*(isi2));
fullrAverage1(fullrAverage1==0)=1;
fullproperty=ceil(property_step*(property));
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(isi2)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,8)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('isi2 vs y coordinat - density plot')
%xlabel('sp.Frequency rAverage1')
ylabel('property')
set(gca,'Position',[.42 .02 .19 .44])


%% histo rAverage1 vs angle - density plot                                 000rAverage1 VS isi2
property=rAverage1;
property_step=10;
fullrAverage1=ceil(density_plot_scaler*(isi1));
fullrAverage1(fullrAverage1==0)=1;
fullproperty=ceil(property_step*(property));
fullproperty(fullproperty<0)=0;
fullproperty=fullproperty+1;
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(isi1)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,5)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('isi1 vs rAverage1 - density plot')
%xlabel('sp.Frequency rAverage1')
ylabel('property')
set(gca,'Position',[.82 .51 .16 .44])


%% histo rAverage2 vs angle - density plot                                 000rAverage1 VS isi2
fullrAverage1=ceil(density_plot_scaler*(isi2));
fullrAverage1(fullrAverage1==0)=1;
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(isi2)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,10)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('isi2 vs rAverage1 - density plot')
%xlabel('sp.Frequency rAverage1')
ylabel('property')
set(gca,'Position',[.82 .02 .16 .44])


%% histo rAverage1 vs angle - density plot                                 angle recBee isi1
property=angle(:,1);
property_step=0.5;
fullrAverage1=ceil(density_plot_scaler*(isi1));
fullrAverage1(fullrAverage1==0)=1;
fullproperty=ceil(property_step*(property));
fullproperty(fullproperty<0)=0;
fullproperty=fullproperty+1;
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(isi1)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,4)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('isi1 vs angle recBee - density plot')
%xlabel('sp.Frequency rAverage1')
ylabel('property')
set(gca,'Position',[.62 .51 .19 .44])


%% histo rAverage2 vs angle - density plot                                 angle recBee isi2
fullrAverage1=ceil(density_plot_scaler*(isi2));
fullrAverage1(fullrAverage1==0)=1;
fullproperty(isnan(fullproperty))=1;
fullrAverage1(isnan(fullrAverage1))=1;
hist_block=zeros(max(fullrAverage1), max(fullproperty));
for i=1:length(isi2)
    hist_block(fullrAverage1(i),fullproperty(i))=hist_block(fullrAverage1(i),fullproperty(i))+1;
end
hist_block(1,1)=0;
subplot(2,5,9)
surface(sqrt(hist_block(2:end,2:end)'),'edgecolor','none') 
title('isi2 vs angle RecBee - density plot')
%xlabel('sp.Frequency rAverage1')
ylabel('property')
set(gca,'Position',[.62 .02 .19 .44])


%% rAverage change up or down quickly! plots [46 & 47 & 48 & 49]
rAverage1_up=rAverage1;
for i=10:length(rAverage1)
if rAverage1(i)-rAverage1(i-1) > .2    % if spikes go up
rAverage1_up(i)=rAverage1_up(i);
else
    rAverage1_up(i)=nan;
end
end
rAverage1_down=rAverage1;
for i=10:length(rAverage1)
if rAverage1(i-1)-rAverage1(i) > .2    % if spikes go down
rAverage1_down(i)=rAverage1_down(i);
else
    rAverage1_down(i)=nan;
end
end

dist_test_up=distance_to_main_sort(:,2);
dist_test_down=distance_to_main_sort(:,2);
dist_test_neither=zeros(length(x),1);
for i=1:length(x)
if isnan(rAverage1_up(i))
    dist_test_up(i)=0;
end
end
for i=1:length(x)
if isnan(rAverage1_down(i))
    dist_test_down(i)=0;
end
end
for i=1:length(x)
if dist_test_down(i) == 0 && dist_test_up(i) == 0
    dist_test_neither(i)=distance_to_main_sort(i,2);
end
end
figure('OuterPosition',figure_positiondouble)
subplot(2,1,1)
bin_distance_to_main_thisplot=30;
zeroline=zeros(bin_distance_to_main_thisplot,1);
hist([dist_test_up dist_test_down],bin_distance_to_main_thisplot)
upper_lim=(mean(hist([dist_test_up dist_test_down],bin_distance_to_main_thisplot)));
ylim([0 upper_lim(2)])
title('spike rate rise (blue) or falls(red) over distence to main SORT histo')
subplot(2,1,2)
hold on
spike_change_hist1=hist([dist_test_up dist_test_down],bin_distance_to_main_thisplot);
plot((spike_change_hist1(:,1)-spike_change_hist1(:,2)));
plot(zeroline)
hold off
title('spike rate rise (positiv) or falls (negativ) over distence to main SORT histo')


dist_test_up=x(:,1);
dist_test_down=x(:,1);
dist_test_neither=zeros(length(x),1);
for i=1:length(x)
if isnan(rAverage1_up(i))
    dist_test_up(i)=0;
end
end
for i=1:length(x)
if isnan(rAverage1_down(i))
    dist_test_down(i)=0;
end
end
for i=1:length(x)
if dist_test_down(i) == 0 && dist_test_up(i) == 0
    dist_test_neither(i)=x(i,1);
end
end
figure('OuterPosition',figure_positiondouble)
subplot(2,1,1)
bin_x_thisplot=300;
zeroline=zeros(bin_x_thisplot,1);
hist([dist_test_up dist_test_down],bin_x_thisplot)
upper_lim=(mean(hist([dist_test_up dist_test_down],bin_x_thisplot)));
ylim([0 upper_lim(1)])
title('spike rate rise (blue) or falls(red) over x-position recBee histo')
subplot(2,1,2)
hold on
spike_change_hist1=hist([dist_test_up dist_test_down],bin_x_thisplot);
plot((spike_change_hist1(:,1)-spike_change_hist1(:,2)));
plot(zeroline)
hold off
title('spike rate rise (positiv) or falls (negativ) over x-position recBee histo')


dist_test_up=x_rel(:,2);
dist_test_down=x_rel(:,2);
dist_test_neither=zeros(length(x),1);
for i=1:length(x)
if isnan(rAverage1_up(i))
    dist_test_up(i)=0;
end
end
for i=1:length(x)
if isnan(rAverage1_down(i))
    dist_test_down(i)=0;
end
end
for i=1:length(x)
if dist_test_down(i) == 0 && dist_test_up(i) == 0
    dist_test_neither(i)=x_rel(i,2);
end
end
figure('OuterPosition',figure_positiondouble)
subplot(2,1,1)
bin_x_rel_thisplot=30;
zeroline=zeros(bin_x_rel_thisplot,1);
hist([dist_test_up dist_test_down],bin_x_rel_thisplot)
upper_lim=(mean(hist([dist_test_up dist_test_down],bin_x_rel_thisplot)));
ylim([0 upper_lim(2)])
title('spike rate rise (blue) or falls(red) over x-position perifer RELATIVE histo')
subplot(2,1,2)
hold on
spike_change_hist1=hist([dist_test_up dist_test_down],bin_x_rel_thisplot);
plot((spike_change_hist1(:,1)-spike_change_hist1(:,2)));
plot(zeroline)
hold off
title('spike rate rise (positiv) or falls (negativ) over x-position perifer RELATIVE histo')


dist_test_up=angle(:,1);
dist_test_down=angle(:,1);
dist_test_neither=zeros(length(x),1);
for i=1:length(x)
if isnan(rAverage1_up(i))
    dist_test_up(i)=0;
end
end
for i=1:length(x)
if isnan(rAverage1_down(i))
    dist_test_down(i)=0;
end
end
for i=1:length(x)
if dist_test_down(i) == 0 && dist_test_up(i) == 0
    dist_test_neither(i)=angle(i,2);
end
end
figure('OuterPosition',figure_positiondouble)
subplot(2,1,1)
bin_angle_thisplot=30;
zeroline=zeros(bin_angle_thisplot,1);
hist([dist_test_up dist_test_down],bin_angle_thisplot);
upper_lim=(mean(hist([dist_test_up dist_test_down],bin_angle_thisplot)));
ylim([0 upper_lim(2)])
title('spike rate rise (blue) or falls(red) over angle of RecBee histo')
subplot(2,1,2)
hold on
spike_change_hist1=hist([dist_test_up dist_test_down],bin_angle_thisplot);
plot((spike_change_hist1(:,1)-spike_change_hist1(:,2)));
plot(zeroline)
hold off
title('spike rate rise (positiv) or falls (negativ) over angle of RecBee histo')


%% Activity1 per angel_rel [50 & 51]
figure('OuterPosition',figure_position1)
winkel_perifer=mod(rad2deg(cart2pol(combiTest_rel(:,2,2),combiTest_rel(:,2,3))),360);
angle_resolution_corr=360/angle_resolution;
angleUnits_rel=NaN(floor((length(x)/100)*angle_resolution_corr),angle_resolution);
for j=1:angle_resolution
    h=1;
    for i=1:length(x)
        if winkel_perifer(i) > (j-1)*angle_resolution_corr && winkel_perifer(i) <= j*angle_resolution_corr
            angleUnits_rel(h,j)=rAverage1(i);
            h=h+1;
        end
    end
end
subplot(2,1,1)
boxplot(angleUnits_rel,'plotstyle','compact','whisker',3)
ylabel('distro of spike Activity Unit1')
xlabel('angle of closest periferBee in degree')
title('amount of spikes per frame of Unit 1 per angle_REL periferBee')
subplot(2,1,2)
angleUnits_rose_mean=zeros(angle_resolution,1);
for i=1:int64(angle_resolution)
angleUnits_rose_mean(i)=nanmedian(angleUnits_rel(:,i));
end
angleUnits_rose=zeros(max(floor(angleUnits_rose_mean*100)),1);
j=1;
for i=1:int64(angle_resolution)
    for k=1:floor(100*angleUnits_rose_mean(i))
angleUnits_rose(j)=i;
j=j+1;
    end
end
subplot(2,1,2)
p=rose(-0.01+deg2rad(10*angleUnits_rose),angle_resolution);
h = findall(gca,'type','line');
h(h == p) = [];
delete(h);

% Activity2 per angel
figure('OuterPosition',figure_position2)
angle_resolution_corr=360/angle_resolution;
angleUnits_rel=NaN(floor((length(x)/100)*angle_resolution_corr),angle_resolution);
for j=1:angle_resolution
    h=1;
    for i=1:length(x)
        if winkel_perifer(i) > (j-1)*angle_resolution_corr && winkel_perifer(i) <= j*angle_resolution_corr
            angleUnits_rel(h,j)=rAverage2(i);
            h=h+1;
        end
    end
end
subplot(2,1,1)
boxplot(angleUnits_rel,'plotstyle','compact','whisker',3)
ylabel('distro of spike Activity Unit1')
xlabel('angle of closest periferBee in degree')
title('amount of spikes per frame of Unit 2 per angle_REL periferBee')
subplot(2,1,2)
angleUnits_rose_mean=zeros(angle_resolution,1);
for i=1:int64(angle_resolution)
angleUnits_rose_mean(i)=nanmedian(angleUnits_rel(:,i));
end
angleUnits_rose=zeros(max(floor(angleUnits_rose_mean*100)),1);
j=1;
for i=1:int64(angle_resolution)
    for k=1:floor(100*angleUnits_rose_mean(i))
angleUnits_rose(j)=i;
j=j+1;
    end
end
subplot(2,1,2)
p=rose(-0.01+deg2rad(10*angleUnits_rose),angle_resolution);
h = findall(gca,'type','line');
h(h == p) = [];
delete(h);





figure('OuterPosition',figure_position2)

scatter(winkel_perifer,combiTest_rel(:,2,1),300,rAverage1,'.')






%% headdircetion implementation   | Track of recBee with headdirection as blue line [52]
figure('OuterPosition',figure_positiondouble)
arrows=nan(length(x),5);    % x recbee y recbee trash(z-axis) y arrowhead x arrowhead
arrows(:,1)=x(:,1);
arrows(:,2)=y(:,1);
hold on
%arrows(:,4)=y(:,1)+1;       % one up so 0°, gets turned by angle(i,1) later

for i=1:length(x)-1
    if isnan(angle(i,1))
    else
R=rotx(angle(i,1));
arrows(i,3:5)=R*[0;10;0];      % building the arrow starting from [0 0 0] %% 'their-x' axis first, empty
arrows(i,4)=arrows(i,4)+arrows(i,2); % adding the arrow to the current pos of RecBee  
arrows(i,5)=arrows(i,5)+arrows(i,1); % y and x are twisted
arrows(i,3)=arrows(i,2);
arrows(i,2)=arrows(i,5);
plot(arrows(i,1:2),arrows(i,3:4),'Color',[((rAverage1(i)/max2)),0,(1-(rAverage1(i)/max2))])
    end
end
hold off
title('headdirection on /as track of Recbee')


%% headdircetion implementation   | Track of periferBees with headdirection as blue line [53]
figure('OuterPosition',figure_positiondouble)
arrows=nan(length(x),5);    % x recbee y recbee trash(z-axis) y arrowhead x arrowhead
arrows(:,1)=combiTest(:,2,2);
arrows(:,2)=combiTest(:,2,3);
hold on
%arrows(:,4)=y(:,1)+1;       % one up so 0°, gets turned by angle(i,1) later
for i=1:length(x)
    if isnan(combiTest(i,2,4))
    else
R=rotx(combiTest(i,2,4));
arrows(i,3:5)=R*[0;10;0];      % building the arrow starting from [0 0 0] %% 'their-x' axis first, empty
arrows(i,4)=arrows(i,4)+arrows(i,2); % adding the arrow to the current pos of RecBee  
arrows(i,5)=arrows(i,5)+arrows(i,1); % y and x are twisted
arrows(i,3)=arrows(i,2);
arrows(i,2)=arrows(i,5);
plot(arrows(i,1:2),arrows(i,3:4))
    end
end
scatter(combiTest(:,2,2),combiTest(:,2,3),[],log10(rAverage1),'filled')
hold off
title('headdirection on track of periferBees with sp.Activity rAverage1Recbee')


%% rAverage1 hist over time [54]
figure('OuterPosition',figure_positiondouble)
subplot(2,2,1)
clear histoarray;
k=1;
bins=100;
under_even=floor(-1+length(isi1)/300)*300;
histoarray=nan(under_even/300,bins);
for i=1:300:under_even
histoarray(k,:)=hist(rAverage1(i:i+499),bins);
k=k+1;
end
histoarray=log10(histoarray);
histoarray(histoarray==-Inf)=0;
contourf(histoarray,200,'EdgeColor','none');
title('rAverage1 hist over time ROOTED')


%% rAverage2 hist over time
subplot(2,2,2)
clear histoarray;
k=1;
under_even=floor(-1+length(isi1)/300)*300;
histoarray=nan(under_even/300,bins);
for i=1:300:under_even
histoarray(k,:)=hist(rAverage2(i:i+499),bins);
k=k+1;
end
histoarray=log10(histoarray);
histoarray(histoarray==-Inf)=0;
contourf(histoarray,200,'EdgeColor','none');
title('rAverage2 hist over time ROOTED')


%% isi1 hist over time
subplot(2,2,3)
clear histoarray;
k=1;
under_even=floor(-1+length(isi1)/300)*300;
histoarray=nan(under_even/300,bins);
for i=1:300:under_even
histoarray(k,:)=hist(isi1(i:i+499),bins);
k=k+1;
end
histoarray=log10(histoarray);
histoarray(histoarray==-Inf)=0;
contourf(histoarray,20,'EdgeColor','none');
title('isi1 hist over time ROOTED')


%% isi2 hist over time
subplot(2,2,4)
clear histoarray;
k=1;
under_even=floor(-1+length(isi1)/300)*300;
histoarray=nan(under_even/300,bins);
for i=1:300:under_even
histoarray(k,:)=hist(isi2(i:i+499),bins);
k=k+1;
end
histoarray=log10(histoarray);
histoarray(histoarray==-Inf)=0;
contourf(histoarray,20,'EdgeColor','none');
title('isi2 hist over time ROOTED')


%% rAverage3/4 hist over time [55]
figure('OuterPosition',figure_positiondouble)
subplot(2,2,1)
clear histoarray;
k=1;
bins=100;
under_even=floor(-1+length(isi3)/300)*300;
histoarray=nan(under_even/300,bins);
for i=1:300:under_even
histoarray(k,:)=hist(rAverage3(i:i+499),bins);
k=k+1;
end
histoarray=log10(histoarray);
histoarray(histoarray==-Inf)=0;
contourf(histoarray,200,'EdgeColor','none');
title('rAverage3 hist over time ROOTED')


%% rAverage2 hist over time
subplot(2,2,2)
clear histoarray;
k=1;
under_even=floor(-1+length(isi3)/300)*300;
histoarray=nan(under_even/300,bins);
for i=1:300:under_even
histoarray(k,:)=hist(rAverage4(i:i+499),bins);
k=k+1;
end
histoarray=log10(histoarray);
histoarray(histoarray==-Inf)=0;
contourf(histoarray,200,'EdgeColor','none');
title('rAverage4 hist over time ROOTED')


%% isi1 hist over time
subplot(2,2,3)
clear histoarray;
k=1;
under_even=floor(-1+length(isi1)/300)*300;
histoarray=nan(under_even/300,bins);
for i=1:300:under_even
histoarray(k,:)=hist(isi3(i:i+499),bins);
k=k+1;
end
histoarray=log10(histoarray);
histoarray(histoarray==-Inf)=0;
contourf(histoarray,20,'EdgeColor','none');
title('isi3 hist over time ROOTED')


%% isi2 hist over time
subplot(2,2,4)
clear histoarray;
k=1;
under_even=floor(-1+length(isi1)/300)*300;
histoarray=nan(under_even/300,bins);
for i=1:300:under_even
histoarray(k,:)=hist(isi4(i:i+499),bins);
k=k+1;
end
histoarray=log10(histoarray);
histoarray(histoarray==-Inf)=0;
contourf(histoarray,20,'EdgeColor','none');
title('isi4 hist over time ROOTED')


%% isi1.*rAverage1*10 1-4 timePlot [56]   ||construction
figure('OuterPosition',figure_position1)
subplot(2,2,1)
plot(isi1corr)
title('isi1corr histogramm')
subplot(2,2,2)
plot(isi2corr)
title('isi2corr histogramm')
subplot(2,2,3)
plot(isi3corr)
title('isi3corr histogramm')
subplot(2,2,4)
plot(isi4corr)
title('isi4corr histogramm')


%% hist(isi1.*rAverage1*10,1000) 1-4  [57] || construction
figure('OuterPosition',figure_position2)
subplot(2,2,1)
hist(isi1corr,1000)
%ylim([0 150])
%xlim([0 20])
title('isi1corr histogramm')
subplot(2,2,2)
hist(isi2corr,1000)
%ylim([0 150])
%xlim([0 5])
title('isi2corr histogramm')
subplot(2,2,3)
hist(isi3corr,1000)
%ylim([0 150])
%xlim([0 5])
title('isi3corr histogramm')
subplot(2,2,4)
hist(isi4corr,1000)
%ylim([0 150])
%xlim([0 5])
title('isi4corr histogramm')
% ISI1=isi1.*rAverage1*10;
% below=ISI1(ISI1<1);
% abow=ISI1(ISI1>1);
% length(below)-length(abow)
% ISI1=isi2.*rAverage2*10;
% below=ISI1(ISI1<1);
% abow=ISI1(ISI1>1);
% length(below)-length(abow)
% ISI1=isi3.*rAverage3*10;
% below=ISI1(ISI1<1);
% abow=ISI1(ISI1>1);
% length(below)-length(abow)
% ISI1=isi4.*rAverage4*10;
% below=ISI1(ISI1<1);
% abow=ISI1(ISI1>1);
% length(below)-length(abow)


%% walking speed vs rAverage1 (without wrong high speeds) [58]
figure('OuterPosition',figure_position1)
subplot(2,1,1)
timeRow=linspace(1,length(x),length(x));
rA_speed_rA1=[rAverageSpeed_rec rAverage1 timeRow'];
rA_speed_rA1=sortrows(rA_speed_rA1,1);
clear rA_speed_rA1_binned
binning_speed=100;
for i=1:binning_speed
    rA_speed_rA1_binned(:,:,i)=rA_speed_rA1(1+(i-1)*floor(length(rA_speed_rA1)/binning_speed):i*floor(length(rA_speed_rA1)/binning_speed),:);
end
rA_speed_rA1_binned_speed=squeeze(rA_speed_rA1_binned(:,2,:));
boxplot(rA_speed_rA1_binned_speed)
title('rA1 vs Speed of RecBee')
rA_speed_rA1(:,3)=(rA_speed_rA1(:,3)/1000)+50;
subplot(2,1,2)
plot(rA_speed_rA1)%(:,1:2))


%% put the original back in place revers sort the figures
isi1=isi1_keep;
isi2=isi2_keep;
isi3=isi3_keep;
isi4=isi4_keep;
rAverage1=rAverage1_keep;
rAverage2=rAverage2_keep;
rAverage3=rAverage3_keep;
rAverage4=rAverage4_keep;
j=gcf;
j=j.Number;
for i=1:j
figure(j)
j=j-1;
end
%% save the workspace || inactive ||
%save('workspaceFULL_5rA.mat');
analyse_time=toc;
disp(['    ',num2str(floor(analyse_time)), ' sec'])
%%################################################  END OF LTS


%% Save all The Pics (doent work good!!) || inactive ||
% 
% for i=1:40
%     bueld=getframe;
%     imwrite(bueld.cdata,[num2str(i),'bild.png'])
%     close
% end


%% trash || inactive ||
%{ 

h = animatedline;
axis([0,4*pi,-1,1])

x = linspace(0,4*pi,1000);
y = sin(x);
for k = 1:length(x)
    addpoints(h,x(k),y(k));
    drawnow
end


t = 0:.01:2*pi;
x = cos(2*t).*(cos(t).^2);
y = sin(2*t).*(sin(t).^2);
comet(isi1,rAverage1);

% contacts plotten (sollten nah sein)

% periferBees sorten für weniger artefakte  (sdiehe oben)

hold on
plot(behavior_contact(:,:,1)',behavior_contact(:,:,2)')

plot(behavior_contact(:,:,6)',behavior_contact(:,:,7)')

hold off

t = zeros(1,100);
for n = 1:100
    A = rand(n,n);
    b = rand(n,1);
    tic;
    x = A\b;
    t(n) = toc;
end
plot(t)

% 
% h = get(0,'children');
% for i=1:length(h)
%   saveas(h(i), ['figure' num2str(i)], 'png');
% end
% 
% 
% 

% hist([rAverage1 rAverage2],60)


plot(arrows(10000,1:2),arrows(10000,3:4))
angle(10000,1)
arrows(10000,:)
plot(arrows(1,1:2),arrows(1,3:4))
i=99;
R=rotx(i);
x=[1;-1];
y(i,:)=R*x;
temp_b(2)
arrows(:,3)=arrows(:,2);
arrows(:,2)=arrows(:,5);
R=rotx(7)
R*[0; 10;0]
figure(89)
hold on
plot(x(20000:20040,1),y(20000:20040,1),'.r')
plot(x(20000:20040,2),y(20000:20040,2),'.b')

plot(x_rel(20000:20040,1),y_rel(20000:20040,1),'.g')
plot(x_rel(20000:20040,2),y_rel(20000:20040,2),'.y')
hold off

plot(angle_rel(20000:20040,2))

hold on
for j=2:12
scatter(y_rel(:,j),x_rel(:,j),'MarkerFaceColor',[j/12,j/12,0],'Marker','o','MarkerEdgeColor','none')
end
hold off

i=99;
R=rotx(i);
x=[1;-1];
y(i,:)=R*x;
temp_b(2)

plot(y(i,2),y(i,3),'MarkerFaceColor',[i/360,i/360,0/360],'Marker','o','MarkerEdgeColor','none')
%}
