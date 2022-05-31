function ace=ace_analysisSpindles()
% by Q all stage
% Please cite 
% Naji, M., Krishnan, G. P., McDevitt, E. A., Bazhenov, M., & Mednick, S. C. (2019). 
% Coupling of autonomic and central events during sleep benefits declarative memory consolidation. Neurobiology of learning and memory, 157, 139-150.

% Chen, P. C., Whitehurst, L. N., Naji, M., & Mednick, S. C. (2020). 
% Autonomic/central coupling benefits working memory in healthy young adults. Neurobiology of Learning and Memory, 173, 107267.

% requires 3 files: sleep edf file, single column text file for sleep
% stages (one for every 30 sec), and mat file of K_u_b_i_o_s output of
% ECG R peaks 

% --------------parameters
segmin=3; % at least 3 minutes for a stable stage
hrbwin=20; % window size around HR burst
avbin=5;  % bin size for averaging EEG powers
m=1.25; % threshold for HR burst detection mean(RR)-m*std(RR)
oc=0; %1 / 0 outlier control on eeg powers
swa=[0.1 4]; alpha=[8 13]; theta=[4 7]; sigma=[12 16]; % EEG bands
lf=[0.04 0.15]; hf=[0.15 0.4]; % LF, HF bands
swa1=[0.5 1];
swa2=[1 4];
RR_upperlim = 2;
RR_lowerlim = 0.55; 
%%%

[FileName,PathName] = uigetfile([pwd '/*.edf'],'Select the sleep edf file');
disp('Reading EDF file ... ');
[X,Channels,fs]=edf2mat_samefss([PathName FileName]);
% [Selection,ok] = listdlg('ListString',Channels,'Name','EEG electrode selection');
for i=1:length(Channels)
cn(i)=Channels{1,i}(1);
end
cn1=ismember(cn,['F' 'C' 'O' 'P' 'f' 'c' 'o' 'p' 'T']);
Selection=find(cn1==1);

XEEG=X(Selection,:);
Xothers=X;Xothers(Selection,:)=[]; clear X;
%
disp('filtering...');
[ba2,aa2]=butter(4,swa(1)/(fs/2),'high');
[ba1,aa1]=butter(4,swa(2)/(fs/2),'low');
X_delta=filtfilt(ba2,aa2,XEEG'); X_delta=filtfilt(ba1,aa1,X_delta); 
X_deltaPWR=abs(hilbert(X_delta)); X_delta=X_delta'; X_deltaPWR=X_deltaPWR';

[ba2,aa2]=butter(4,swa1(1)/(fs/2),'high');
[ba1,aa1]=butter(4,swa1(2)/(fs/2),'low');
X_delta1=filtfilt(ba2,aa2,XEEG'); X_delta1=filtfilt(ba1,aa1,X_delta1); 
X_delta1PWR=abs(hilbert(X_delta1)); X_delta1=X_delta1'; X_delta1PWR=X_delta1PWR';

[ba2,aa2]=butter(4,swa2(1)/(fs/2),'high');
[ba1,aa1]=butter(4,swa2(2)/(fs/2),'low');
X_delta2=filtfilt(ba2,aa2,XEEG'); X_delta2=filtfilt(ba1,aa1,X_delta2); 
X_delta2PWR=abs(hilbert(X_delta2)); X_delta2=X_delta2'; X_delta2PWR=X_delta2PWR';


[ba4,aa4]=butter(4,alpha(1)/(fs/2),'high');
[ba3,aa3]=butter(4,alpha(2)/(fs/2),'low'); 
X_alpha=filtfilt(ba4,aa4,XEEG'); X_alpha=filtfilt(ba3,aa3,X_alpha); 
X_alphaPWR=abs(hilbert(X_alpha)); X_alpha=X_alpha'; X_alphaPWR=X_alphaPWR';

[ba6,aa6]=butter(4,theta(1)/(fs/2),'high');
[ba5,aa5]=butter(4,theta(2)/(fs/2),'low');
X_theta=filtfilt(ba6,aa6,XEEG'); X_theta=filtfilt(ba5,aa5,X_theta);
X_thetaPWR=abs(hilbert(X_theta)); X_theta=X_theta'; X_thetaPWR=X_thetaPWR';

[ba8,aa8]=butter(4,sigma(1)/(fs/2),'high');
[ba7,aa7]=butter(4,sigma(2)/(fs/2),'low'); 
X_sigma=filtfilt(ba8,aa8,XEEG'); X_sigma=filtfilt(ba7,aa7,X_sigma); 
X_sigmaPWR=abs(hilbert(X_sigma)); X_sigma=X_sigma'; X_sigmaPWR=X_sigmaPWR';
%
% eval(['x' Channels{1, 6}  '= X(6,:);']);


% read RR intervals file
tempallfilesMat = dir(fullfile([PathName '*.mat']));
alln_mat={tempallfilesMat.name};
ii=zeros(1,length(alln_mat));
for i=1:length(alln_mat)
ii(i)=sum(alln_mat{1,i}(1:length(FileName(1:end-4)))==FileName(1:end-4));
end
rf=find(ii==length(FileName(1:end-4)));
s=strcat(PathName, alln_mat(rf));
load(s{1,1});
if ~exist('RES')
RES=Res;
end
clear Res
fgr=fieldnames(RES.CNT.rate);
fs2=eval(['RES.CNT.rate' '.' fgr{1,1}]);
RR_tot_ind=round((RES.HRV.Data.T_RR-RES.CNT.Offset)*fs2);
RR_tot=(RR_tot_ind(2:end)-RR_tot_ind(1:end-1))./fs2;
RR_tot_ind(1)=[];
RR_tot_time=RR_tot_ind./fs2;
RRts=spline(RR_tot_ind./fs2,RR_tot,1/fs2:1/fs2:RR_tot_ind(end)/fs2);
% RR filter
[blf2,alf2]=butter(4,lf(2)/(fs/2),'low');
[blf1,alf1]=butter(4,lf(1)/(fs/2),'high');
RRlf=filtfilt(blf2,alf2,RRts);  RRlf=filtfilt(blf1,alf1,RRlf);

[bhf2,ahf2]=butter(4,hf(1)/(fs/2),'high');
[bhf1,ahf1]=butter(4,hf(2)/(fs/2),'low');
RRhf=filtfilt(bhf2,ahf2,RRts);  RRhf=filtfilt(bhf1,ahf1,RRhf);
hfpw=abs(hilbert(RRhf));
%
disp('RR loaded')
if fs2~=fs
    disp(['Warning! is sampling rate ' fs ' or ' fs2 '?']);
end
% read sleep score file
tempallfilesTxt = dir(fullfile([PathName '*.txt']));
alln_txt={tempallfilesTxt.name}
FileName
ii=zeros(1,length(alln_txt));
for i=1:length(alln_txt)
ii(i)=sum(alln_txt{1,i}(1:length(FileName(1:end-4)))==FileName(1:end-4));
end
rf=find(ii==length(FileName(1:end-4)));
s=strcat(PathName, alln_txt(rf));
mrkr=load(s{1,1});
mrkr(mrkr==-1)=7;
disp('scores loaded')
t=find((mrkr(2:end)-mrkr(1:end-1))~=0);
smp=[0;t*30*fs;size(mrkr,1)*30*fs]';
bnd=zeros(length(smp)-1,2);
for i=1:length(smp)-1
    bnd(i,:)=[smp(i)+1 smp(i+1)]; % beginning and end of each bout
end
Stage=mrkr([1;t+1]); % bout sleep stage
% assignin('base','bnd',bnd);

duration=(bnd(:,2)-bnd(:,1)+1)/(fs*60); 
stg2segs=find(Stage==2 & bnd(:,2)<length(RRts));
stg3segs=find(Stage==3 & bnd(:,2)<length(RRts));

ch_stg2 = cell(1,length(Selection));
ch_spindles_bnd = cell(1,length(Selection));
ch_spindles_quartile = cell(1,length(Selection));
ch_spindles_dur = cell(1,length(Selection));
if ~isempty(stg2segs)
for ch=1:length(Selection)%size(stg2segs,2)
    for sg = 1:length(stg2segs)
    ch_stg2{sg,ch}=XEEG(ch,bnd(stg2segs(sg),1):bnd(stg2segs(sg),2));%if ~isempty(find(X(f3ch,bnd(stg2segs(i),1):bnd(stg2segs(i),2))>380)) rj=[rj i]; end
    end
    [~,ch_spindle_duration,ch_spindle_bnd] = Mywamsley_spindle_detection(ch_stg2(:,ch),XEEG(ch,:),fs);
    ch_spindles_bnd{1,ch} = ch_spindle_bnd;
    ch_spindles_dur{1,ch} = ch_spindle_duration;
end
end

assignin('base','ch_stg2',ch_stg2);
assignin('base','stg2segs',stg2segs);
assignin('base','Stage',Stage);
assignin('base','Stage',Stage);
assignin('base','ch_spindles_bnd',ch_spindles_bnd);
assignin('base','ch_spindles_dur',ch_spindles_dur);

% bouts longer than segmin at each stage
in7=find(Stage==7 & duration>=segmin & bnd(:,2)<length(RRts)); 
in0=find(Stage==0 & duration>=segmin & bnd(:,2)<length(RRts));
in1=find(Stage==1 & duration>=segmin & bnd(:,2)<length(RRts));
in2=find(Stage==2 & duration>=segmin & bnd(:,2)<length(RRts));
in3=find(Stage==3 & duration>=segmin & bnd(:,2)<length(RRts));
in5=find(Stage==5 & duration>=segmin & bnd(:,2)<length(RRts));
% assignin('base','RRts',RRts);

o=find(mrkr>0 & mrkr<7); e=find(mrkr(end:-1:1)>0 & mrkr(end:-1:1)<7);
sleep_beg=(o(1)-1)*30*fs+1;
sleep_end=(length(mrkr)-e(1)+1)*30*fs;
q = 4;
ql=(sleep_end-sleep_beg+1)/q;

ii = 1;q1_s=(ii-1)*ql+sleep_beg; q1_e=ii*ql+sleep_beg-1;
ii = 2;q2_s=(ii-1)*ql+sleep_beg; q2_e=ii*ql+sleep_beg-1;
ii = 3;q3_s=(ii-1)*ql+sleep_beg; q3_e=ii*ql+sleep_beg-1;
ii = 4;q4_s=(ii-1)*ql+sleep_beg; q4_e=ii*ql+sleep_beg-1;

disp('bouts extracted');
%
disp('HR burst analysis...')
[sbj_hrb_ind7,sbj_hrb7,sbj_hrbHFPW7]=myhrbwindowedd(in7,bnd,fs,hrbwin,RRts,hfpw,m,RR_lowerlim, RR_upperlim); sbj_hrb_HFpwr_bins7=myhrbbinnedd_EEG(sbj_hrbHFPW7,avbin,fs,oc);
sbj_hrb_RR_bins7=myhrbbinnedd_EEG(sbj_hrb7,avbin,fs,oc);
sbj_hrb_EEG7=myhrbwindowedd_EEG(XEEG,hrbwin,fs,oc,sbj_hrb_ind7); 
sbj_hrb_deltaPWR7=myhrbwindowedd_EEG(X_deltaPWR,hrbwin,fs,oc,sbj_hrb_ind7); sbj_hrb_deltapwr_bins7=myhrbbinnedd_EEG(sbj_hrb_deltaPWR7,avbin,fs,oc);
sbj_hrb_alphaPWR7=myhrbwindowedd_EEG(X_alphaPWR,hrbwin,fs,oc,sbj_hrb_ind7); sbj_hrb_alphapwr_bins7=myhrbbinnedd_EEG(sbj_hrb_alphaPWR7,avbin,fs,oc);
sbj_hrb_thetaPWR7=myhrbwindowedd_EEG(X_thetaPWR,hrbwin,fs,oc,sbj_hrb_ind7); sbj_hrb_thetapwr_bins7=myhrbbinnedd_EEG(sbj_hrb_thetaPWR7,avbin,fs,oc);
sbj_hrb_sigmaPWR7=myhrbwindowedd_EEG(X_sigmaPWR,hrbwin,fs,oc,sbj_hrb_ind7); sbj_hrb_sigmapwr_bins7=myhrbbinnedd_EEG(sbj_hrb_sigmaPWR7,avbin,fs,oc);
sbj_hrb_delta1PWR7=myhrbwindowedd_EEG(X_delta1PWR,hrbwin,fs,oc,sbj_hrb_ind7); sbj_hrb_delta1pwr_bins7=myhrbbinnedd_EEG(sbj_hrb_delta1PWR7,avbin,fs,oc);
sbj_hrb_delta2PWR7=myhrbwindowedd_EEG(X_delta2PWR,hrbwin,fs,oc,sbj_hrb_ind7); sbj_hrb_delta2pwr_bins7=myhrbbinnedd_EEG(sbj_hrb_delta2PWR7,avbin,fs,oc);

[sbj_hrb_ind0,sbj_hrb0,sbj_hrbHFPW0]=myhrbwindowedd(in0,bnd,fs,hrbwin,RRts,hfpw,m,RR_lowerlim, RR_upperlim); sbj_hrb_HFpwr_bins0=myhrbbinnedd_EEG(sbj_hrbHFPW0,avbin,fs,oc);
sbj_hrb_RR_bins0=myhrbbinnedd_EEG(sbj_hrb0,avbin,fs,oc);
sbj_hrb_EEG0=myhrbwindowedd_EEG(XEEG,hrbwin,fs,oc,sbj_hrb_ind0);
sbj_hrb_deltaPWR0=myhrbwindowedd_EEG(X_deltaPWR,hrbwin,fs,oc,sbj_hrb_ind0); sbj_hrb_deltapwr_bins0=myhrbbinnedd_EEG(sbj_hrb_deltaPWR0,avbin,fs,oc);
sbj_hrb_alphaPWR0=myhrbwindowedd_EEG(X_alphaPWR,hrbwin,fs,oc,sbj_hrb_ind0); sbj_hrb_alphapwr_bins0=myhrbbinnedd_EEG(sbj_hrb_alphaPWR0,avbin,fs,oc);
sbj_hrb_thetaPWR0=myhrbwindowedd_EEG(X_thetaPWR,hrbwin,fs,oc,sbj_hrb_ind0); sbj_hrb_thetapwr_bins0=myhrbbinnedd_EEG(sbj_hrb_thetaPWR0,avbin,fs,oc);
sbj_hrb_sigmaPWR0=myhrbwindowedd_EEG(X_sigmaPWR,hrbwin,fs,oc,sbj_hrb_ind0); sbj_hrb_sigmapwr_bins0=myhrbbinnedd_EEG(sbj_hrb_sigmaPWR0,avbin,fs,oc);
sbj_hrb_delta1PWR0=myhrbwindowedd_EEG(X_delta1PWR,hrbwin,fs,oc,sbj_hrb_ind0); sbj_hrb_delta1pwr_bins0=myhrbbinnedd_EEG(sbj_hrb_delta1PWR0,avbin,fs,oc);
sbj_hrb_delta2PWR0=myhrbwindowedd_EEG(X_delta2PWR,hrbwin,fs,oc,sbj_hrb_ind0); sbj_hrb_delta2pwr_bins0=myhrbbinnedd_EEG(sbj_hrb_delta2PWR0,avbin,fs,oc);

[sbj_hrb_ind1,sbj_hrb1,sbj_hrbHFPW1]=myhrbwindowedd(in1,bnd,fs,hrbwin,RRts,hfpw,m,RR_lowerlim, RR_upperlim); sbj_hrb_HFpwr_bins1=myhrbbinnedd_EEG(sbj_hrbHFPW1,avbin,fs,oc);
sbj_hrb_RR_bins1=myhrbbinnedd_EEG(sbj_hrb1,avbin,fs,oc);
sbj_hrb_EEG1=myhrbwindowedd_EEG(XEEG,hrbwin,fs,oc,sbj_hrb_ind1);
sbj_hrb_deltaPWR1=myhrbwindowedd_EEG(X_deltaPWR,hrbwin,fs,oc,sbj_hrb_ind1); sbj_hrb_deltapwr_bins1=myhrbbinnedd_EEG(sbj_hrb_deltaPWR1,avbin,fs,oc);
sbj_hrb_alphaPWR1=myhrbwindowedd_EEG(X_alphaPWR,hrbwin,fs,oc,sbj_hrb_ind1); sbj_hrb_alphapwr_bins1=myhrbbinnedd_EEG(sbj_hrb_alphaPWR1,avbin,fs,oc);
sbj_hrb_thetaPWR1=myhrbwindowedd_EEG(X_thetaPWR,hrbwin,fs,oc,sbj_hrb_ind1); sbj_hrb_thetapwr_bins1=myhrbbinnedd_EEG(sbj_hrb_thetaPWR1,avbin,fs,oc);
sbj_hrb_sigmaPWR1=myhrbwindowedd_EEG(X_sigmaPWR,hrbwin,fs,oc,sbj_hrb_ind1); sbj_hrb_sigmapwr_bins1=myhrbbinnedd_EEG(sbj_hrb_sigmaPWR1,avbin,fs,oc);
sbj_hrb_delta1PWR1=myhrbwindowedd_EEG(X_delta1PWR,hrbwin,fs,oc,sbj_hrb_ind1); sbj_hrb_delta1pwr_bins1=myhrbbinnedd_EEG(sbj_hrb_delta1PWR1,avbin,fs,oc);
sbj_hrb_delta2PWR1=myhrbwindowedd_EEG(X_delta2PWR,hrbwin,fs,oc,sbj_hrb_ind1); sbj_hrb_delta2pwr_bins1=myhrbbinnedd_EEG(sbj_hrb_delta2PWR1,avbin,fs,oc);

[sbj_hrb_ind2,sbj_hrb2,sbj_hrbHFPW2]=myhrbwindowedd(in2,bnd,fs,hrbwin,RRts,hfpw,m,RR_lowerlim, RR_upperlim); sbj_hrb_HFpwr_bins2=myhrbbinnedd_EEG(sbj_hrbHFPW2,avbin,fs,oc);
sbj_hrb_RR_bins2=myhrbbinnedd_EEG(sbj_hrb2,avbin,fs,oc);
sbj_hrb_EEG2=myhrbwindowedd_EEG(XEEG,hrbwin,fs,oc,sbj_hrb_ind2);
sbj_hrb_deltaPWR2=myhrbwindowedd_EEG(X_deltaPWR,hrbwin,fs,oc,sbj_hrb_ind2); sbj_hrb_deltapwr_bins2=myhrbbinnedd_EEG(sbj_hrb_deltaPWR2,avbin,fs,oc);
sbj_hrb_alphaPWR2=myhrbwindowedd_EEG(X_alphaPWR,hrbwin,fs,oc,sbj_hrb_ind2); sbj_hrb_alphapwr_bins2=myhrbbinnedd_EEG(sbj_hrb_alphaPWR2,avbin,fs,oc);
sbj_hrb_thetaPWR2=myhrbwindowedd_EEG(X_thetaPWR,hrbwin,fs,oc,sbj_hrb_ind2); sbj_hrb_thetapwr_bins2=myhrbbinnedd_EEG(sbj_hrb_thetaPWR2,avbin,fs,oc);
sbj_hrb_sigmaPWR2=myhrbwindowedd_EEG(X_sigmaPWR,hrbwin,fs,oc,sbj_hrb_ind2); sbj_hrb_sigmapwr_bins2=myhrbbinnedd_EEG(sbj_hrb_sigmaPWR2,avbin,fs,oc);
sbj_hrb_delta1PWR2=myhrbwindowedd_EEG(X_delta1PWR,hrbwin,fs,oc,sbj_hrb_ind2); sbj_hrb_delta1pwr_bins2=myhrbbinnedd_EEG(sbj_hrb_delta1PWR2,avbin,fs,oc);
sbj_hrb_delta2PWR2=myhrbwindowedd_EEG(X_delta2PWR,hrbwin,fs,oc,sbj_hrb_ind2); sbj_hrb_delta2pwr_bins2=myhrbbinnedd_EEG(sbj_hrb_delta2PWR2,avbin,fs,oc);
sbj_hrb_spindles2 = myHRB_spindle_intersect(sbj_hrb_ind2, ch_spindles_bnd, hrbwin,fs);
for c = 1:length(sbj_hrb_spindles2)
    sbj_hrb_spindles_densityPerMinute2{1,c} = mean(sbj_hrb_spindles2{1,c},1)*(60/avbin); 
end
assignin('base','sbj_hrb_spindles2',sbj_hrb_spindles2);

[sbj_hrb_ind3,sbj_hrb3,sbj_hrbHFPW3]=myhrbwindowedd(in3,bnd,fs,hrbwin,RRts,hfpw,m,RR_lowerlim, RR_upperlim); sbj_hrb_HFpwr_bins3=myhrbbinnedd_EEG(sbj_hrbHFPW3,avbin,fs,oc);
sbj_hrb_RR_bins3=myhrbbinnedd_EEG(sbj_hrb3,avbin,fs,oc);
sbj_hrb_EEG3=myhrbwindowedd_EEG(XEEG,hrbwin,fs,oc,sbj_hrb_ind3);
sbj_hrb_deltaPWR3=myhrbwindowedd_EEG(X_deltaPWR,hrbwin,fs,oc,sbj_hrb_ind3); sbj_hrb_deltapwr_bins3=myhrbbinnedd_EEG(sbj_hrb_deltaPWR3,avbin,fs,oc);
sbj_hrb_alphaPWR3=myhrbwindowedd_EEG(X_alphaPWR,hrbwin,fs,oc,sbj_hrb_ind3); sbj_hrb_alphapwr_bins3=myhrbbinnedd_EEG(sbj_hrb_alphaPWR3,avbin,fs,oc);
sbj_hrb_thetaPWR3=myhrbwindowedd_EEG(X_thetaPWR,hrbwin,fs,oc,sbj_hrb_ind3); sbj_hrb_thetapwr_bins3=myhrbbinnedd_EEG(sbj_hrb_thetaPWR3,avbin,fs,oc);
sbj_hrb_sigmaPWR3=myhrbwindowedd_EEG(X_sigmaPWR,hrbwin,fs,oc,sbj_hrb_ind3); sbj_hrb_sigmapwr_bins3=myhrbbinnedd_EEG(sbj_hrb_sigmaPWR3,avbin,fs,oc);
sbj_hrb_delta1PWR3=myhrbwindowedd_EEG(X_delta1PWR,hrbwin,fs,oc,sbj_hrb_ind3); sbj_hrb_delta1pwr_bins3=myhrbbinnedd_EEG(sbj_hrb_delta1PWR3,avbin,fs,oc);
sbj_hrb_delta2PWR3=myhrbwindowedd_EEG(X_delta2PWR,hrbwin,fs,oc,sbj_hrb_ind3); sbj_hrb_delta2pwr_bins3=myhrbbinnedd_EEG(sbj_hrb_delta2PWR3,avbin,fs,oc);
sbj_hrb_spindles3 = myHRB_spindle_intersect(sbj_hrb_ind3, ch_spindles_bnd, hrbwin,fs);
for c = 1:length(sbj_hrb_spindles3)
    sbj_hrb_spindles_densityPerMinute3{1,c} = mean(sbj_hrb_spindles3{1,c},1)*(60/avbin); 
end
assignin('base','sbj_hrb_spindles3',sbj_hrb_spindles3);

[sbj_hrb_ind5,sbj_hrb5,sbj_hrbHFPW5]=myhrbwindowedd(in5,bnd,fs,hrbwin,RRts,hfpw,m,RR_lowerlim, RR_upperlim); sbj_hrb_HFpwr_bins5=myhrbbinnedd_EEG(sbj_hrbHFPW5,avbin,fs,oc);
sbj_hrb_RR_bins5=myhrbbinnedd_EEG(sbj_hrb5,avbin,fs,oc);
sbj_hrb_EEG5=myhrbwindowedd_EEG(XEEG,hrbwin,fs,oc,sbj_hrb_ind5);
sbj_hrb_deltaPWR5=myhrbwindowedd_EEG(X_deltaPWR,hrbwin,fs,oc,sbj_hrb_ind5); sbj_hrb_deltapwr_bins5=myhrbbinnedd_EEG(sbj_hrb_deltaPWR5,avbin,fs,oc);
sbj_hrb_alphaPWR5=myhrbwindowedd_EEG(X_alphaPWR,hrbwin,fs,oc,sbj_hrb_ind5); sbj_hrb_alphapwr_bins5=myhrbbinnedd_EEG(sbj_hrb_alphaPWR5,avbin,fs,oc);
sbj_hrb_thetaPWR5=myhrbwindowedd_EEG(X_thetaPWR,hrbwin,fs,oc,sbj_hrb_ind5); sbj_hrb_thetapwr_bins5=myhrbbinnedd_EEG(sbj_hrb_thetaPWR5,avbin,fs,oc);
sbj_hrb_sigmaPWR5=myhrbwindowedd_EEG(X_sigmaPWR,hrbwin,fs,oc,sbj_hrb_ind5); sbj_hrb_sigmapwr_bins5=myhrbbinnedd_EEG(sbj_hrb_sigmaPWR5,avbin,fs,oc);
sbj_hrb_delta1PWR5=myhrbwindowedd_EEG(X_delta1PWR,hrbwin,fs,oc,sbj_hrb_ind5); sbj_hrb_delta1pwr_bins5=myhrbbinnedd_EEG(sbj_hrb_delta1PWR5,avbin,fs,oc);
sbj_hrb_delta2PWR5=myhrbwindowedd_EEG(X_delta2PWR,hrbwin,fs,oc,sbj_hrb_ind5); sbj_hrb_delta2pwr_bins5=myhrbbinnedd_EEG(sbj_hrb_delta2PWR5,avbin,fs,oc);

sbj_hrb_ind0_quartile = zeros(size(sbj_hrb_ind0));
sbj_hrb_ind1_quartile = zeros(size(sbj_hrb_ind1));
sbj_hrb_ind2_quartile = zeros(size(sbj_hrb_ind2));
sbj_hrb_ind3_quartile = zeros(size(sbj_hrb_ind3));
sbj_hrb_ind5_quartile = zeros(size(sbj_hrb_ind5));
sbj_hrb_ind7_quartile = zeros(size(sbj_hrb_ind7));
bnd_q = zeros(size(bnd,1),1);

for b = 1:length(bnd_q)
    if bnd(b,1)>=q1_s && bnd(b,1)<=q1_e
        bnd_q(b) = 1;
    elseif bnd(b,1)>=q2_s && bnd(b,1)<=q2_e
        bnd_q(b) = 2;
    elseif bnd(b,1)>=q3_s && bnd(b,1)<=q3_e
        bnd_q(b) = 3;
    elseif bnd(b,1)>=q4_s && bnd(b,1)<=q4_e
        bnd_q(b) = 4;
    end
end
minutes_q_wake = [0,0,0,0];
minutes_q_stg1 = [0,0,0,0];
minutes_q_stg2 = [0,0,0,0];
minutes_q_sws = [0,0,0,0];
minutes_q_rem = [0,0,0,0];
for qq = 1:4
    minutes_q_wake(qq) = sum(duration(find(Stage==0 & duration>=segmin & bnd_q == qq)));
    minutes_q_stg1(qq) = sum(duration(find(Stage==1 & duration>=segmin & bnd_q == qq)));
    minutes_q_stg2(qq) = sum(duration(find(Stage==2 & duration>=segmin & bnd_q == qq)));
    minutes_q_sws(qq) = sum(duration(find(Stage==3 & duration>=segmin & bnd_q == qq)));
    minutes_q_rem(qq) = sum(duration(find(Stage==5 & duration>=segmin & bnd_q == qq)));
end
ace.minutes_q_wake = minutes_q_wake;
ace.minutes_q_stg1 = minutes_q_stg1;
ace.minutes_q_stg2 = minutes_q_stg2;
ace.minutes_q_sws = minutes_q_sws;
ace.minutes_q_rem = minutes_q_rem;

for ix = 1:length(sbj_hrb_ind0)
    if sbj_hrb_ind0(ix)>=q1_s && sbj_hrb_ind0(ix)<=q1_e
        sbj_hrb_ind0_quartile(ix) = 1;
    elseif sbj_hrb_ind0(ix)>=q2_s && sbj_hrb_ind0(ix)<=q2_e
        sbj_hrb_ind0_quartile(ix) = 2;
    elseif sbj_hrb_ind0(ix)>=q3_s && sbj_hrb_ind0(ix)<=q3_e
        sbj_hrb_ind0_quartile(ix) = 3;
    elseif sbj_hrb_ind0(ix)>=q4_s 
        sbj_hrb_ind0_quartile(ix) = 4;
    elseif sbj_hrb_ind0(ix)<q1_s 
        sbj_hrb_ind0_quartile(ix) = 0;
    end   
end

for ix = 1:length(sbj_hrb_ind1)
    if sbj_hrb_ind1(ix)>=q1_s && sbj_hrb_ind1(ix)<=q1_e
        sbj_hrb_ind1_quartile(ix) = 1;
    elseif sbj_hrb_ind1(ix)>=q2_s && sbj_hrb_ind1(ix)<=q2_e
        sbj_hrb_ind1_quartile(ix) = 2;
    elseif sbj_hrb_ind1(ix)>=q3_s && sbj_hrb_ind1(ix)<=q3_e
        sbj_hrb_ind1_quartile(ix) = 3;
    elseif sbj_hrb_ind1(ix)>=q4_s 
        sbj_hrb_ind1_quartile(ix) = 4;
    elseif sbj_hrb_ind1(ix)<q1_s 
        sbj_hrb_ind1_quartile(ix) = 0;
    end   
end

for ix = 1:length(sbj_hrb_ind2)
    if sbj_hrb_ind2(ix)>=q1_s && sbj_hrb_ind2(ix)<=q1_e
        sbj_hrb_ind2_quartile(ix) = 1;
    elseif sbj_hrb_ind2(ix)>=q2_s && sbj_hrb_ind2(ix)<=q2_e
        sbj_hrb_ind2_quartile(ix) = 2;
    elseif sbj_hrb_ind2(ix)>=q3_s && sbj_hrb_ind2(ix)<=q3_e
        sbj_hrb_ind2_quartile(ix) = 3;
    elseif sbj_hrb_ind2(ix)>=q4_s 
        sbj_hrb_ind2_quartile(ix) = 4;
    elseif sbj_hrb_ind2(ix)<q1_s 
        sbj_hrb_ind2_quartile(ix) = 0;
    end   
end

for ix = 1:length(sbj_hrb_ind3)
    if sbj_hrb_ind3(ix)>=q1_s && sbj_hrb_ind3(ix)<=q1_e
        sbj_hrb_ind3_quartile(ix) = 1;
    elseif sbj_hrb_ind3(ix)>=q2_s && sbj_hrb_ind3(ix)<=q2_e
        sbj_hrb_ind3_quartile(ix) = 2;
    elseif sbj_hrb_ind3(ix)>=q3_s && sbj_hrb_ind3(ix)<=q3_e
        sbj_hrb_ind3_quartile(ix) = 3;
    elseif sbj_hrb_ind3(ix)>=q4_s 
        sbj_hrb_ind3_quartile(ix) = 4;
    elseif sbj_hrb_ind3(ix)<q1_s 
        sbj_hrb_ind3_quartile(ix) = 0;
    end   
end

for ix = 1:length(sbj_hrb_ind5)
    if sbj_hrb_ind5(ix)>=q1_s && sbj_hrb_ind5(ix)<=q1_e
        sbj_hrb_ind5_quartile(ix) = 1;
    elseif sbj_hrb_ind5(ix)>=q2_s && sbj_hrb_ind5(ix)<=q2_e
        sbj_hrb_ind5_quartile(ix) = 2;
    elseif sbj_hrb_ind5(ix)>=q3_s && sbj_hrb_ind5(ix)<=q3_e
        sbj_hrb_ind5_quartile(ix) = 3;
    elseif sbj_hrb_ind5(ix)>=q4_s 
        sbj_hrb_ind5_quartile(ix) = 4;
    elseif sbj_hrb_ind5(ix)<q1_s 
        sbj_hrb_ind5_quartile(ix) = 0;
    end   
end

for ix = 1:length(sbj_hrb_ind7)
    if sbj_hrb_ind7(ix)>=q1_s && sbj_hrb_ind7(ix)<=q1_e
        sbj_hrb_ind7_quartile(ix) = 1;
    elseif sbj_hrb_ind7(ix)>=q2_s && sbj_hrb_ind7(ix)<=q2_e
        sbj_hrb_ind7_quartile(ix) = 2;
    elseif sbj_hrb_ind7(ix)>=q3_s && sbj_hrb_ind7(ix)<=q3_e
        sbj_hrb_ind7_quartile(ix) = 3;
    elseif sbj_hrb_ind7(ix)>=q4_s 
        sbj_hrb_ind7_quartile(ix) = 4;
    elseif sbj_hrb_ind7(ix)<q1_s 
        sbj_hrb_ind7_quartile(ix) = 0;
    end   
end

for ch=1:length(Selection)
    for ix = 1:size(ch_spindles_bnd{1,ch},1)
        if ch_spindles_bnd{1,ch}(ix,1)>=q1_s && ch_spindles_bnd{1,ch}(ix,1)<=q1_e
            ch_spindles_quartile{1,ch}(ix,1) = 1;
        elseif ch_spindles_bnd{1,ch}(ix,1)>=q2_s && ch_spindles_bnd{1,ch}(ix,1)<=q2_e
            ch_spindles_quartile{1,ch}(ix,1) = 2;
        elseif ch_spindles_bnd{1,ch}(ix,1)>=q3_s && ch_spindles_bnd{1,ch}(ix,1)<=q3_e
            ch_spindles_quartile{1,ch}(ix,1) = 3;
        elseif ch_spindles_bnd{1,ch}(ix,1)>=q4_s 
            ch_spindles_quartile{1,ch}(ix,1) = 4;
        elseif ch_spindles_bnd{1,ch}(ix,1)<q1_s 
            ch_spindles_quartile{1,ch}(ix,1) = 0;
        end   
    end
end
assignin('base','ch_spindles_quartile',ch_spindles_quartile);

%%%%
disp('calculating bin data...');
eval('ace.Channels= Channels(1, Selection);');
for i=1:length(Selection)
    eval(['ace.binDelta_stg1_' Channels{1, Selection(i)}  '= sbj_hrb_deltapwr_bins1{i,1};']);
    eval(['ace.binDelta_sws_' Channels{1, Selection(i)}  '= sbj_hrb_deltapwr_bins3{i,1};']);
    eval(['ace.binDelta_stg2_' Channels{1, Selection(i)}  '= sbj_hrb_deltapwr_bins2{i,1};']);
    eval(['ace.binDelta_wake_' Channels{1, Selection(i)}  '= sbj_hrb_deltapwr_bins0{i,1};']);
    eval(['ace.binDelta_rem_' Channels{1, Selection(i)}  '= sbj_hrb_deltapwr_bins5{i,1};']);
    eval(['ace.binDelta_nostage_' Channels{1, Selection(i)}  '= sbj_hrb_deltapwr_bins7{i,1};']);
    
    eval(['ace.binSlowDelta_stg1_' Channels{1, Selection(i)}  '= sbj_hrb_delta1pwr_bins1{i,1};']);
    eval(['ace.binSlowDelta_sws_' Channels{1, Selection(i)}  '= sbj_hrb_delta1pwr_bins3{i,1};']);
    eval(['ace.binSlowDelta_stg2_' Channels{1, Selection(i)}  '= sbj_hrb_delta1pwr_bins2{i,1};']);
    eval(['ace.binSlowDelta_wake_' Channels{1, Selection(i)}  '= sbj_hrb_delta1pwr_bins0{i,1};']);
    eval(['ace.binSlowDelta_rem_' Channels{1, Selection(i)}  '= sbj_hrb_delta1pwr_bins5{i,1};']);
    eval(['ace.binSlowDelta_nostage_' Channels{1, Selection(i)}  '= sbj_hrb_delta1pwr_bins7{i,1};']);
    
    eval(['ace.binFastDelta_stg1_' Channels{1, Selection(i)}  '= sbj_hrb_delta2pwr_bins1{i,1};']);
    eval(['ace.binFastDelta_sws_' Channels{1, Selection(i)}  '= sbj_hrb_delta2pwr_bins3{i,1};']);
    eval(['ace.binFastDelta_stg2_' Channels{1, Selection(i)}  '= sbj_hrb_delta2pwr_bins2{i,1};']);
    eval(['ace.binFastDelta_wake_' Channels{1, Selection(i)}  '= sbj_hrb_delta2pwr_bins0{i,1};']);
    eval(['ace.binFastDelta_rem_' Channels{1, Selection(i)}  '= sbj_hrb_delta2pwr_bins5{i,1};']);
    eval(['ace.binFastDelta_nostage_' Channels{1, Selection(i)}  '= sbj_hrb_delta2pwr_bins7{i,1};']);
    
    eval(['ace.binAlpha_sws_' Channels{1, Selection(i)}  '= sbj_hrb_alphapwr_bins3{i,1};']);
    eval(['ace.binAlpha_stg2_' Channels{1, Selection(i)}  '= sbj_hrb_alphapwr_bins2{i,1};']);
    eval(['ace.binAlpha_stg1_' Channels{1, Selection(i)}  '= sbj_hrb_alphapwr_bins1{i,1};']);
    eval(['ace.binAlpha_wake_' Channels{1, Selection(i)}  '= sbj_hrb_alphapwr_bins0{i,1};']);
    eval(['ace.binAlpha_rem_' Channels{1, Selection(i)}  '= sbj_hrb_alphapwr_bins5{i,1};']);
    eval(['ace.binAlpha_nostage_' Channels{1, Selection(i)}  '= sbj_hrb_alphapwr_bins7{i,1};']);
    
    eval(['ace.binSigma_sws_' Channels{1, Selection(i)}  '= sbj_hrb_sigmapwr_bins3{i,1};']);
    eval(['ace.binSigma_stg2_' Channels{1, Selection(i)}  '= sbj_hrb_sigmapwr_bins2{i,1};']);
    eval(['ace.binSigma_stg1_' Channels{1, Selection(i)}  '= sbj_hrb_sigmapwr_bins1{i,1};']);
    eval(['ace.binSigma_wake_' Channels{1, Selection(i)}  '= sbj_hrb_sigmapwr_bins0{i,1};']);
    eval(['ace.binSigma_rem_' Channels{1, Selection(i)}  '= sbj_hrb_sigmapwr_bins5{i,1};']);
    eval(['ace.binSigma_nostage_' Channels{1, Selection(i)}  '= sbj_hrb_sigmapwr_bins7{i,1};']);
    
    eval(['ace.binTheta_sws_' Channels{1, Selection(i)}  '= sbj_hrb_thetapwr_bins3{i,1};']);
    eval(['ace.binTheta_stg2_' Channels{1, Selection(i)}  '= sbj_hrb_thetapwr_bins2{i,1};']);
    eval(['ace.binTheta_stg1_' Channels{1, Selection(i)}  '= sbj_hrb_thetapwr_bins1{i,1};']);
    eval(['ace.binTheta_wake_' Channels{1, Selection(i)}  '= sbj_hrb_thetapwr_bins0{i,1};']);
    eval(['ace.binTheta_rem_' Channels{1, Selection(i)}  '= sbj_hrb_thetapwr_bins5{i,1};']);
    eval(['ace.binTheta_nostage_' Channels{1, Selection(i)}  '= sbj_hrb_thetapwr_bins7{i,1};']);
    eval(['ace.binSpindleDensity_stg2_' Channels{1, Selection(i)} '= sbj_hrb_spindles_densityPerMinute2{1,i};']);
    eval(['ace.binSpindleCount_stg2_' Channels{1, Selection(i)} '= sbj_hrb_spindles2{1,i};']);
    eval(['ace.binSpindleDensity_sws_' Channels{1, Selection(i)} '= sbj_hrb_spindles_densityPerMinute3{1,i};']);
    eval(['ace.binSpindleCount_sws_' Channels{1, Selection(i)} '= sbj_hrb_spindles3{1,i};']);
end
eval(['ace.binHF_sws'  '= sbj_hrb_HFpwr_bins3;']);
eval(['ace.binHF_stg2' '= sbj_hrb_HFpwr_bins2;']);
eval(['ace.binHF_stg1'   '= sbj_hrb_HFpwr_bins1;']);
eval(['ace.binHF_wake'  '= sbj_hrb_HFpwr_bins0;']);
eval(['ace.binHF_rem'   '= sbj_hrb_HFpwr_bins5;']);
eval(['ace.binHF_nostage'  '= sbj_hrb_HFpwr_bins7;']);

eval(['ace.binRR_sws'  '= sbj_hrb_RR_bins3;']);
eval(['ace.binRR_stg2' '= sbj_hrb_RR_bins2;']);
eval(['ace.binRR_stg1'   '= sbj_hrb_RR_bins1;']);
eval(['ace.binRR_wake'  '= sbj_hrb_RR_bins0;']);
eval(['ace.binRR_rem'   '= sbj_hrb_RR_bins5;']);
eval(['ace.binRR_nostage'  '= sbj_hrb_RR_bins7;']);
%%%%
disp('averaging...');%(find(sbj_hrb_ind2_quartile ==1),:)
for i=1:length(Selection)
    eval(['ace.avDelta_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins1{i,1}, 1:size(sbj_hrb_deltapwr_bins1{i,1},1));']);
    eval(['ace.avDelta_Q1_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==1));']);
    eval(['ace.avDelta_Q2_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==2));']);
    eval(['ace.avDelta_Q3_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==3));']);
    eval(['ace.avDelta_Q4_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==4));']);
    
    eval(['ace.avDelta_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins3{i,1}, 1:size(sbj_hrb_deltapwr_bins3{i,1},1));']);
    eval(['ace.avDelta_Q1_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==1));']);
    eval(['ace.avDelta_Q2_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==2));']);
    eval(['ace.avDelta_Q3_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==3));']);
    eval(['ace.avDelta_Q4_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==4));']);
    
    eval(['ace.avDelta_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins2{i,1}, 1:size(sbj_hrb_deltapwr_bins2{i,1},1));']);
    eval(['ace.avDelta_Q1_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==1));']);
    eval(['ace.avDelta_Q2_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==2));']);
    eval(['ace.avDelta_Q3_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==3));']);
    eval(['ace.avDelta_Q4_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==4));']);
    
    eval(['ace.avDelta_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins0{i,1}, 1:size(sbj_hrb_deltapwr_bins0{i,1},1));']);
    eval(['ace.avDelta_Q1_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==1));']);
    eval(['ace.avDelta_Q2_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==2));']);
    eval(['ace.avDelta_Q3_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==3));']);
    eval(['ace.avDelta_Q4_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==4));']);
    
    eval(['ace.avDelta_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins5{i,1}, 1:size(sbj_hrb_deltapwr_bins5{i,1},1));']);
    eval(['ace.avDelta_Q1_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==1));']);
    eval(['ace.avDelta_Q2_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==2));']);
    eval(['ace.avDelta_Q3_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==3));']);
    eval(['ace.avDelta_Q4_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==4));']);
    
    eval(['ace.avDelta_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins7{i,1}, 1:size(sbj_hrb_deltapwr_bins7{i,1},1));']);
    eval(['ace.avDelta_Q1_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==1));']);
    eval(['ace.avDelta_Q2_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==2));']);
    eval(['ace.avDelta_Q3_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==3));']);
    eval(['ace.avDelta_Q4_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==4));']);
    
    eval(['ace.avSlowDelta_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins1{i,1}, 1:size(sbj_hrb_delta1pwr_bins1{i,1},1));']);
    eval(['ace.avSlowDelta_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins3{i,1}, 1:size(sbj_hrb_delta1pwr_bins3{i,1},1));']);
    eval(['ace.avSlowDelta_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins2{i,1}, 1:size(sbj_hrb_delta1pwr_bins2{i,1},1));']);
    eval(['ace.avSlowDelta_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins0{i,1}, 1:size(sbj_hrb_delta1pwr_bins0{i,1},1));']);
    eval(['ace.avSlowDelta_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins5{i,1}, 1:size(sbj_hrb_delta1pwr_bins5{i,1},1));']);
    eval(['ace.avSlowDelta_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins7{i,1}, 1:size(sbj_hrb_delta1pwr_bins7{i,1},1));']);
    
    eval(['ace.avFastDelta_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins1{i,1}, 1:size(sbj_hrb_delta2pwr_bins1{i,1},1));']);
    eval(['ace.avFastDelta_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins3{i,1}, 1:size(sbj_hrb_delta2pwr_bins3{i,1},1));']);
    eval(['ace.avFastDelta_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins2{i,1}, 1:size(sbj_hrb_delta2pwr_bins2{i,1},1));']);
    eval(['ace.avFastDelta_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins0{i,1}, 1:size(sbj_hrb_delta2pwr_bins0{i,1},1));']);
    eval(['ace.avFastDelta_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins5{i,1}, 1:size(sbj_hrb_delta2pwr_bins5{i,1},1));']);
    eval(['ace.avFastDelta_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins7{i,1}, 1:size(sbj_hrb_delta2pwr_bins7{i,1},1));']);
    
    eval(['ace.avAlpha_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins3{i,1}, 1:size(sbj_hrb_alphapwr_bins3{i,1},1));']);
    eval(['ace.avAlpha_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins2{i,1}, 1:size(sbj_hrb_alphapwr_bins2{i,1},1));']);
    eval(['ace.avAlpha_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins1{i,1}, 1:size(sbj_hrb_alphapwr_bins1{i,1},1));']);
    eval(['ace.avAlpha_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins0{i,1}, 1:size(sbj_hrb_alphapwr_bins0{i,1},1));']);
    eval(['ace.avAlpha_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins5{i,1}, 1:size(sbj_hrb_alphapwr_bins5{i,1},1));']);
    eval(['ace.avAlpha_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins7{i,1}, 1:size(sbj_hrb_alphapwr_bins7{i,1},1));']);
    
    eval(['ace.avSigma_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins3{i,1}, 1:size(sbj_hrb_sigmapwr_bins3{i,1},1));']);
    eval(['ace.avSigma_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins2{i,1}, 1:size(sbj_hrb_sigmapwr_bins2{i,1},1));']);
    eval(['ace.avSigma_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins1{i,1}, 1:size(sbj_hrb_sigmapwr_bins1{i,1},1));']);
    eval(['ace.avSigma_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins0{i,1}, 1:size(sbj_hrb_sigmapwr_bins0{i,1},1));']);
    eval(['ace.avSigma_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins5{i,1}, 1:size(sbj_hrb_sigmapwr_bins5{i,1},1));']);
    eval(['ace.avSigma_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins7{i,1}, 1:size(sbj_hrb_sigmapwr_bins7{i,1},1));']);
    
    eval(['ace.avTheta_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins3{i,1}, 1:size(sbj_hrb_thetapwr_bins3{i,1},1));']);
    eval(['ace.avTheta_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins2{i,1}, 1:size(sbj_hrb_thetapwr_bins2{i,1},1));']);
    eval(['ace.avTheta_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins1{i,1}, 1:size(sbj_hrb_thetapwr_bins1{i,1},1));']);
    eval(['ace.avTheta_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins0{i,1}, 1:size(sbj_hrb_thetapwr_bins0{i,1},1));']);
    eval(['ace.avTheta_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins5{i,1}, 1:size(sbj_hrb_thetapwr_bins5{i,1},1));']);
    eval(['ace.avTheta_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins7{i,1}, 1:size(sbj_hrb_thetapwr_bins7{i,1},1));']);
    
    %%%%
    eval(['ace.avSlowDelta_Q1_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==1));']);
    eval(['ace.avSlowDelta_Q2_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==2));']);
    eval(['ace.avSlowDelta_Q3_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==3));']);
    eval(['ace.avSlowDelta_Q4_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==4));']);
    eval(['ace.avSlowDelta_Q1_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==1));']);
    eval(['ace.avSlowDelta_Q2_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==2));']);
    eval(['ace.avSlowDelta_Q3_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==3));']);
    eval(['ace.avSlowDelta_Q4_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==4));']);
    eval(['ace.avSlowDelta_Q1_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==1));']);
    eval(['ace.avSlowDelta_Q2_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==2));']);
    eval(['ace.avSlowDelta_Q3_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==3));']);
    eval(['ace.avSlowDelta_Q4_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==4));']);
    eval(['ace.avSlowDelta_Q1_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==1));']);
    eval(['ace.avSlowDelta_Q2_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==2));']);
    eval(['ace.avSlowDelta_Q3_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==3));']);
    eval(['ace.avSlowDelta_Q4_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==4));']);
    eval(['ace.avSlowDelta_Q1_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==1));']);
    eval(['ace.avSlowDelta_Q2_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==2));']);
    eval(['ace.avSlowDelta_Q3_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==3));']);
    eval(['ace.avSlowDelta_Q4_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==4));']);
    eval(['ace.avSlowDelta_Q1_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==1));']);
    eval(['ace.avSlowDelta_Q2_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==2));']);
    eval(['ace.avSlowDelta_Q3_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==3));']);
    eval(['ace.avSlowDelta_Q4_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==4));']);
    %%%
    eval(['ace.avFastDelta_Q1_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==1));']);
    eval(['ace.avFastDelta_Q2_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==2));']);
    eval(['ace.avFastDelta_Q3_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==3));']);
    eval(['ace.avFastDelta_Q4_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==4));']);
    eval(['ace.avFastDelta_Q1_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==1));']);
    eval(['ace.avFastDelta_Q2_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==2));']);
    eval(['ace.avFastDelta_Q3_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==3));']);
    eval(['ace.avFastDelta_Q4_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==4));']);
    eval(['ace.avFastDelta_Q1_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==1));']);
    eval(['ace.avFastDelta_Q2_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==2));']);
    eval(['ace.avFastDelta_Q3_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==3));']);
    eval(['ace.avFastDelta_Q4_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==4));']);
    eval(['ace.avFastDelta_Q1_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==1));']);
    eval(['ace.avFastDelta_Q2_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==2));']);
    eval(['ace.avFastDelta_Q3_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==3));']);
    eval(['ace.avFastDelta_Q4_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==4));']);
    eval(['ace.avFastDelta_Q1_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==1));']);
    eval(['ace.avFastDelta_Q2_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==2));']);
    eval(['ace.avFastDelta_Q3_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==3));']);
    eval(['ace.avFastDelta_Q4_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==4));']);
    eval(['ace.avFastDelta_Q1_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==1));']);
    eval(['ace.avFastDelta_Q2_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==2));']);
    eval(['ace.avFastDelta_Q3_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==3));']);
    eval(['ace.avFastDelta_Q4_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==4));']);
    %%%
    eval(['ace.avAlpha_Q1_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==1));']);
    eval(['ace.avAlpha_Q2_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==2));']);
    eval(['ace.avAlpha_Q3_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==3));']);
    eval(['ace.avAlpha_Q4_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==4));']);
    eval(['ace.avAlpha_Q1_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==1));']);
    eval(['ace.avAlpha_Q2_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==2));']);
    eval(['ace.avAlpha_Q3_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==3));']);
    eval(['ace.avAlpha_Q4_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==4));']);
    eval(['ace.avAlpha_Q1_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==1));']);
    eval(['ace.avAlpha_Q2_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==2));']);
    eval(['ace.avAlpha_Q3_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==3));']);
    eval(['ace.avAlpha_Q4_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==4));']);
    eval(['ace.avAlpha_Q1_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==1));']);
    eval(['ace.avAlpha_Q2_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==2));']);
    eval(['ace.avAlpha_Q3_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==3));']);
    eval(['ace.avAlpha_Q4_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==4));']);
    eval(['ace.avAlpha_Q1_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==1));']);
    eval(['ace.avAlpha_Q2_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==2));']);
    eval(['ace.avAlpha_Q3_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==3));']);
    eval(['ace.avAlpha_Q4_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==4));']);
    eval(['ace.avAlpha_Q1_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==1));']);
    eval(['ace.avAlpha_Q2_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==2));']);
    eval(['ace.avAlpha_Q3_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==3));']);
    eval(['ace.avAlpha_Q4_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==4));']);
    %%%
    eval(['ace.avSigma_Q1_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==1));']);
    eval(['ace.avSigma_Q2_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==2));']);
    eval(['ace.avSigma_Q3_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==3));']);
    eval(['ace.avSigma_Q4_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==4));']);
    eval(['ace.avSigma_Q1_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==1));']);
    eval(['ace.avSigma_Q2_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==2));']);
    eval(['ace.avSigma_Q3_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==3));']);
    eval(['ace.avSigma_Q4_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==4));']);
    eval(['ace.avSigma_Q1_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==1));']);
    eval(['ace.avSigma_Q2_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==2));']);
    eval(['ace.avSigma_Q3_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==3));']);
    eval(['ace.avSigma_Q4_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==4));']);
    eval(['ace.avSigma_Q1_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==1));']);
    eval(['ace.avSigma_Q2_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==2));']);
    eval(['ace.avSigma_Q3_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==3));']);
    eval(['ace.avSigma_Q4_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==4));']);
    eval(['ace.avSigma_Q1_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==1));']);
    eval(['ace.avSigma_Q2_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==2));']);
    eval(['ace.avSigma_Q3_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==3));']);
    eval(['ace.avSigma_Q4_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==4));']);
    eval(['ace.avSigma_Q1_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==1));']);
    eval(['ace.avSigma_Q2_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==2));']);
    eval(['ace.avSigma_Q3_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==3));']);
    eval(['ace.avSigma_Q4_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==4));']);
    %%%
    eval(['ace.avTheta_Q1_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==1));']);
    eval(['ace.avTheta_Q2_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==2));']);
    eval(['ace.avTheta_Q3_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==3));']);
    eval(['ace.avTheta_Q4_stg1_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins1{i,1},find(sbj_hrb_ind1_quartile ==4));']);
    eval(['ace.avTheta_Q1_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==1));']);
    eval(['ace.avTheta_Q2_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==2));']);
    eval(['ace.avTheta_Q3_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==3));']);
    eval(['ace.avTheta_Q4_sws_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins3{i,1},find(sbj_hrb_ind3_quartile ==4));']);
    eval(['ace.avTheta_Q1_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==1));']);
    eval(['ace.avTheta_Q2_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==2));']);
    eval(['ace.avTheta_Q3_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==3));']);
    eval(['ace.avTheta_Q4_stg2_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins2{i,1},find(sbj_hrb_ind2_quartile ==4));']);
    eval(['ace.avTheta_Q1_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==1));']);
    eval(['ace.avTheta_Q2_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==2));']);
    eval(['ace.avTheta_Q3_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==3));']);
    eval(['ace.avTheta_Q4_wake_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins0{i,1},find(sbj_hrb_ind0_quartile ==4));']);
    eval(['ace.avTheta_Q1_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==1));']);
    eval(['ace.avTheta_Q2_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==2));']);
    eval(['ace.avTheta_Q3_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==3));']);
    eval(['ace.avTheta_Q4_rem_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins5{i,1},find(sbj_hrb_ind5_quartile ==4));']);
    eval(['ace.avTheta_Q1_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==1));']);
    eval(['ace.avTheta_Q2_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==2));']);
    eval(['ace.avTheta_Q3_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==3));']);
    eval(['ace.avTheta_Q4_nostage_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins7{i,1},find(sbj_hrb_ind7_quartile ==4));']);
    %%%
    
end
eval(['ace.avHF_sws'  '= mean_rj_olyr(sbj_hrb_HFpwr_bins3, 1:size(sbj_hrb_HFpwr_bins3,1));']);
eval(['ace.avHF_stg2' '= mean_rj_olyr(sbj_hrb_HFpwr_bins2, 1:size(sbj_hrb_HFpwr_bins2,1));']);
eval(['ace.avHF_stg1'   '= mean_rj_olyr(sbj_hrb_HFpwr_bins1, 1:size(sbj_hrb_HFpwr_bins1,1));']);
eval(['ace.avHF_wake'  '= mean_rj_olyr(sbj_hrb_HFpwr_bins0, 1:size(sbj_hrb_HFpwr_bins0,1));']);
eval(['ace.avHF_rem'   '= mean_rj_olyr(sbj_hrb_HFpwr_bins5, 1:size(sbj_hrb_HFpwr_bins5,1));']);
eval(['ace.avHF_nostage'  '= mean_rj_olyr(sbj_hrb_HFpwr_bins7, 1:size(sbj_hrb_HFpwr_bins7,1));']);

eval(['ace.RRtimeseries'  '= RRts;']);
eval(['ace.HRB_sws_idx'  '= sbj_hrb_ind3;']);
eval(['ace.HRB_stg2_idx' '= sbj_hrb_ind2;']);
eval(['ace.HRB_stg1_idx'   '= sbj_hrb_ind1;']);
eval(['ace.HRB_wake_idx'   '= sbj_hrb_ind0;']);
eval(['ace.HRB_rem_idx'  '= sbj_hrb_ind5;']);
eval(['ace.HRB_nostage_idx'  '= sbj_hrb_ind7;']);

eval(['ace.HRB_sws_quartile'  '= sbj_hrb_ind3_quartile;']);
eval(['ace.HRB_stg2_quartile' '= sbj_hrb_ind2_quartile;']);
eval(['ace.HRB_stg1_quartile'   '= sbj_hrb_ind1_quartile;']);
eval(['ace.HRB_wake_quartile'   '= sbj_hrb_ind0_quartile;']);
eval(['ace.HRB_rem_quartile'  '= sbj_hrb_ind5_quartile;']);
eval(['ace.HRB_nostage_quartile'  '= sbj_hrb_ind7_quartile;']);

eval(['ace.HRB_sws_density'  '= length(sbj_hrb_ind3)/sum(duration(in3));']);
eval(['ace.HRB_stg2_density' '= length(sbj_hrb_ind2)/sum(duration(in2));']);
eval(['ace.HRB_stg1_density'   '= length(sbj_hrb_ind1)/sum(duration(in1));']);
eval(['ace.HRB_wake_density'   '= length(sbj_hrb_ind0)/sum(duration(in0));']);
eval(['ace.HRB_rem_density'  '= length(sbj_hrb_ind5)/sum(duration(in5));']);
eval(['ace.HRB_nostage_density'  '= length(sbj_hrb_ind7)/sum(duration(in7));']);

eval(['ace.HRB_sws'  '= sbj_hrb3;']);
eval(['ace.HRB_stg2' '= sbj_hrb2;']);
eval(['ace.HRB_stg1'   '= sbj_hrb1;']);
eval(['ace.HRB_wake'   '= sbj_hrb0;']);
eval(['ace.HRB_rem'  '= sbj_hrb5;']);
eval(['ace.HRB_nostage'  '= sbj_hrb7;']);

eval(['ace.HRB_EEG_sws_allCh'  '= sbj_hrb_EEG3;']);
eval(['ace.HRB_EEG_stg2_allCh'  '= sbj_hrb_EEG2;']);
eval(['ace.HRB_EEG_stg1_allCh'  '= sbj_hrb_EEG1;']);
eval(['ace.HRB_EEG_wake_allCh'   '= sbj_hrb_EEG0;']);
eval(['ace.HRB_EEG_rem_allCh'  '= sbj_hrb_EEG5;']);
eval(['ace.HRB_EEG_nostage_allCh'   '= sbj_hrb_EEG7;']);

eval(['ace.HRB_DeltaHilbAmp_sws_allCh'  '= sbj_hrb_deltaPWR3;']);
eval(['ace.HRB_DeltaHilbAmp_stg2_allCh'  '= sbj_hrb_deltaPWR2;']);
eval(['ace.HRB_DeltaHilbAmp_stg1_allCh'  '= sbj_hrb_deltaPWR1;']);
eval(['ace.HRB_DeltaHilbAmp_wake_allCh'   '= sbj_hrb_deltaPWR0;']);
eval(['ace.HRB_DeltaHilbAmp_rem_allCh'  '= sbj_hrb_deltaPWR5;']);
eval(['ace.HRB_DeltaHilbAmp_nostage_allCh'   '= sbj_hrb_deltaPWR7;']);

eval(['ace.HRB_SlowDeltaHilbAmp_sws_allCh'  '= sbj_hrb_delta1PWR3;']);
eval(['ace.HRB_SlowDeltaHilbAmp_stg2_allCh'  '= sbj_hrb_delta1PWR2;']);
eval(['ace.HRB_SlowDeltaHilbAmp_stg1_allCh'  '= sbj_hrb_delta1PWR1;']);
eval(['ace.HRB_SlowDeltaHilbAmp_wake_allCh'   '= sbj_hrb_delta1PWR0;']);
eval(['ace.HRB_SlowDeltaHilbAmp_rem_allCh'  '= sbj_hrb_delta1PWR5;']);
eval(['ace.HRB_SlowDeltaHilbAmp_nostage_allCh'   '= sbj_hrb_delta1PWR7;']);

eval(['ace.HRB_FastDeltaHilbAmp_sws_allCh'  '= sbj_hrb_delta2PWR3;']);
eval(['ace.HRB_FastDeltaHilbAmp_stg2_allCh'  '= sbj_hrb_delta2PWR2;']);
eval(['ace.HRB_FastDeltaHilbAmp_stg1_allCh'  '= sbj_hrb_delta2PWR1;']);
eval(['ace.HRB_FastDeltaHilbAmp_wake_allCh'   '= sbj_hrb_delta2PWR0;']);
eval(['ace.HRB_FastDeltaHilbAmp_rem_allCh'  '= sbj_hrb_delta2PWR5;']);
eval(['ace.HRB_FastDeltaHilbAmp_nostage_allCh'   '= sbj_hrb_delta2PWR7;']);


eval(['ace.HRB_SigmaHilbAmp_sws_allCh'  '= sbj_hrb_sigmaPWR3;']);
eval(['ace.HRB_SigmaHilbAmp_stg2_allCh'  '= sbj_hrb_sigmaPWR2;']);
eval(['ace.HRB_SigmaHilbAmp_stg1_allCh'  '= sbj_hrb_sigmaPWR1;']);
eval(['ace.HRB_SigmaHilbAmp_wake_allCh'   '= sbj_hrb_sigmaPWR0;']);
eval(['ace.HRB_SigmaHilbAmp_rem_allCh'  '= sbj_hrb_sigmaPWR5;']);
eval(['ace.HRB_SigmaHilbAmp_nostage_allCh'   '= sbj_hrb_sigmaPWR7;']);

eval(['ace.HRB_ThetaHilbAmp_sws_allCh'  '= sbj_hrb_thetaPWR3;']);
eval(['ace.HRB_ThetaHilbAmp_stg2_allCh'  '= sbj_hrb_thetaPWR2;']);
eval(['ace.HRB_ThetaHilbAmp_stg1_allCh'  '= sbj_hrb_thetaPWR1;']);
eval(['ace.HRB_ThetaHilbAmp_wake_allCh'   '= sbj_hrb_thetaPWR0;']);
eval(['ace.HRB_ThetaHilbAmp_rem_allCh'  '= sbj_hrb_thetaPWR5;']);
eval(['ace.HRB_ThetaHilbAmp_nostage_allCh'   '= sbj_hrb_thetaPWR7;']);

eval(['ace.HRB_AlphaHilbAmp_sws_allCh'  '= sbj_hrb_alphaPWR3;']);
eval(['ace.HRB_AlphaHilbAmp_stg2_allCh'  '= sbj_hrb_alphaPWR2;']);
eval(['ace.HRB_AlphaHilbAmp_stg1_allCh'  '= sbj_hrb_alphaPWR1;']);
eval(['ace.HRB_AlphaHilbAmp_wake_allCh'   '= sbj_hrb_alphaPWR0;']);
eval(['ace.HRB_AlphaHilbAmp_rem_allCh'  '= sbj_hrb_alphaPWR5;']);
eval(['ace.HRB_AlphaHilbAmp_nostage_allCh'   '= sbj_hrb_alphaPWR7;']);

disp('non-Ace analysis...')
inc0=myblactivityy(in0,bnd,sbj_hrb_ind0,hrbwin,fs);
inc1=myblactivityy(in1,bnd,sbj_hrb_ind1,hrbwin,fs);
inc2=myblactivityy(in2,bnd,sbj_hrb_ind2,hrbwin,fs);
inc3=myblactivityy(in3,bnd,sbj_hrb_ind3,hrbwin,fs);
inc5=myblactivityy(in5,bnd,sbj_hrb_ind5,hrbwin,fs);
inc7=myblactivityy(in7,bnd,sbj_hrb_ind7,hrbwin,fs);

inc0_q1=myblactivityy(in0,bnd,sbj_hrb_ind0(find(sbj_hrb_ind0_quartile==1)),hrbwin,fs);
inc0_q2=myblactivityy(in0,bnd,sbj_hrb_ind0(find(sbj_hrb_ind0_quartile==2)),hrbwin,fs);
inc0_q3=myblactivityy(in0,bnd,sbj_hrb_ind0(find(sbj_hrb_ind0_quartile==3)),hrbwin,fs);
inc0_q4=myblactivityy(in0,bnd,sbj_hrb_ind0(find(sbj_hrb_ind0_quartile==4)),hrbwin,fs);

inc1_q1=myblactivityy(in1,bnd,sbj_hrb_ind1(find(sbj_hrb_ind1_quartile==1)),hrbwin,fs);
inc1_q2=myblactivityy(in1,bnd,sbj_hrb_ind1(find(sbj_hrb_ind1_quartile==2)),hrbwin,fs);
inc1_q3=myblactivityy(in1,bnd,sbj_hrb_ind1(find(sbj_hrb_ind1_quartile==3)),hrbwin,fs);
inc1_q4=myblactivityy(in1,bnd,sbj_hrb_ind1(find(sbj_hrb_ind1_quartile==4)),hrbwin,fs);

inc2_q1=myblactivityy(in2,bnd,sbj_hrb_ind2(find(sbj_hrb_ind2_quartile==1)),hrbwin,fs);
inc2_q2=myblactivityy(in2,bnd,sbj_hrb_ind2(find(sbj_hrb_ind2_quartile==2)),hrbwin,fs);
inc2_q3=myblactivityy(in2,bnd,sbj_hrb_ind2(find(sbj_hrb_ind2_quartile==3)),hrbwin,fs);
inc2_q4=myblactivityy(in2,bnd,sbj_hrb_ind2(find(sbj_hrb_ind2_quartile==4)),hrbwin,fs);

inc3_q1=myblactivityy(in3,bnd,sbj_hrb_ind3(find(sbj_hrb_ind3_quartile==1)),hrbwin,fs);
inc3_q2=myblactivityy(in3,bnd,sbj_hrb_ind3(find(sbj_hrb_ind3_quartile==2)),hrbwin,fs);
inc3_q3=myblactivityy(in3,bnd,sbj_hrb_ind3(find(sbj_hrb_ind3_quartile==3)),hrbwin,fs);
inc3_q4=myblactivityy(in3,bnd,sbj_hrb_ind3(find(sbj_hrb_ind3_quartile==4)),hrbwin,fs);

inc5_q1=myblactivityy(in5,bnd,sbj_hrb_ind5(find(sbj_hrb_ind5_quartile==1)),hrbwin,fs);
inc5_q2=myblactivityy(in5,bnd,sbj_hrb_ind5(find(sbj_hrb_ind5_quartile==2)),hrbwin,fs);
inc5_q3=myblactivityy(in5,bnd,sbj_hrb_ind5(find(sbj_hrb_ind5_quartile==3)),hrbwin,fs);
inc5_q4=myblactivityy(in5,bnd,sbj_hrb_ind5(find(sbj_hrb_ind5_quartile==4)),hrbwin,fs);

inc7_q1=myblactivityy(in7,bnd,sbj_hrb_ind7(find(sbj_hrb_ind7_quartile==1)),hrbwin,fs);
inc7_q2=myblactivityy(in7,bnd,sbj_hrb_ind7(find(sbj_hrb_ind7_quartile==2)),hrbwin,fs);
inc7_q3=myblactivityy(in7,bnd,sbj_hrb_ind7(find(sbj_hrb_ind7_quartile==3)),hrbwin,fs);
inc7_q4=myblactivityy(in7,bnd,sbj_hrb_ind7(find(sbj_hrb_ind7_quartile==4)),hrbwin,fs);

bl_spindleDensity_2 = myblspindle(in2,bnd,sbj_hrb_ind2,hrbwin,fs, ch_spindles_bnd);
bl_spindleDensity_3 = myblspindle(in3,bnd,sbj_hrb_ind3,hrbwin,fs, ch_spindles_bnd);

for ch = 1:length(Selection)
    ch_spindles_bnd_q1{1,ch} = ch_spindles_bnd{1,ch}(find(ch_spindles_quartile{1,ch}==1),:);
    ch_spindles_bnd_q2{1,ch} = ch_spindles_bnd{1,ch}(find(ch_spindles_quartile{1,ch}==2),:);
    ch_spindles_bnd_q3{1,ch} = ch_spindles_bnd{1,ch}(find(ch_spindles_quartile{1,ch}==3),:);
    ch_spindles_bnd_q4{1,ch} = ch_spindles_bnd{1,ch}(find(ch_spindles_quartile{1,ch}==4),:);
end
sbj_hrb_spindles2_q1 = myHRB_spindle_intersect(sbj_hrb_ind2, ch_spindles_bnd_q1, hrbwin,fs);
sbj_hrb_spindles2_q2 = myHRB_spindle_intersect(sbj_hrb_ind2, ch_spindles_bnd_q2, hrbwin,fs);
sbj_hrb_spindles2_q3 = myHRB_spindle_intersect(sbj_hrb_ind2, ch_spindles_bnd_q3, hrbwin,fs);
sbj_hrb_spindles2_q4 = myHRB_spindle_intersect(sbj_hrb_ind2, ch_spindles_bnd_q4, hrbwin,fs);

for c = 1:length(sbj_hrb_spindles2)
    sbj_hrb_spindles_densityPerMinute2_q1{1,c} = mean(sbj_hrb_spindles2_q1{1,c},1)*(60/avbin); 
    sbj_hrb_spindles_densityPerMinute2_q2{1,c} = mean(sbj_hrb_spindles2_q2{1,c},1)*(60/avbin); 
    sbj_hrb_spindles_densityPerMinute2_q3{1,c} = mean(sbj_hrb_spindles2_q3{1,c},1)*(60/avbin); 
    sbj_hrb_spindles_densityPerMinute2_q4{1,c} = mean(sbj_hrb_spindles2_q4{1,c},1)*(60/avbin); 
end

sbj_hrb_spindles3_q1 = myHRB_spindle_intersect(sbj_hrb_ind3, ch_spindles_bnd_q1, hrbwin,fs);
sbj_hrb_spindles3_q2 = myHRB_spindle_intersect(sbj_hrb_ind3, ch_spindles_bnd_q2, hrbwin,fs);
sbj_hrb_spindles3_q3 = myHRB_spindle_intersect(sbj_hrb_ind3, ch_spindles_bnd_q3, hrbwin,fs);
sbj_hrb_spindles3_q4 = myHRB_spindle_intersect(sbj_hrb_ind3, ch_spindles_bnd_q4, hrbwin,fs);

for c = 1:length(sbj_hrb_spindles3)
    sbj_hrb_spindles_densityPerMinute3_q1{1,c} = mean(sbj_hrb_spindles3_q1{1,c},1)*(60/avbin); 
    sbj_hrb_spindles_densityPerMinute3_q2{1,c} = mean(sbj_hrb_spindles3_q2{1,c},1)*(60/avbin); 
    sbj_hrb_spindles_densityPerMinute3_q3{1,c} = mean(sbj_hrb_spindles3_q3{1,c},1)*(60/avbin); 
    sbj_hrb_spindles_densityPerMinute3_q4{1,c} = mean(sbj_hrb_spindles3_q4{1,c},1)*(60/avbin); 
end

bl_spindleDensity_2_q1 = myblspindle(in2,bnd,sbj_hrb_ind2,hrbwin,fs, ch_spindles_bnd_q1);
bl_spindleDensity_2_q2 = myblspindle(in2,bnd,sbj_hrb_ind2,hrbwin,fs, ch_spindles_bnd_q2);
bl_spindleDensity_2_q3 = myblspindle(in2,bnd,sbj_hrb_ind2,hrbwin,fs, ch_spindles_bnd_q3);
bl_spindleDensity_2_q4 = myblspindle(in2,bnd,sbj_hrb_ind2,hrbwin,fs, ch_spindles_bnd_q4);

bl_spindleDensity_3_q1 = myblspindle(in3,bnd,sbj_hrb_ind3,hrbwin,fs, ch_spindles_bnd_q1);
bl_spindleDensity_3_q2 = myblspindle(in3,bnd,sbj_hrb_ind3,hrbwin,fs, ch_spindles_bnd_q2);
bl_spindleDensity_3_q3 = myblspindle(in3,bnd,sbj_hrb_ind3,hrbwin,fs, ch_spindles_bnd_q3);
bl_spindleDensity_3_q4 = myblspindle(in3,bnd,sbj_hrb_ind3,hrbwin,fs, ch_spindles_bnd_q4);


for i=1:length(Selection)
    eval(['ace.blDelta_sws_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc3));']);
    eval(['ace.blDelta_stg2_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc2));']);
    eval(['ace.blDelta_stg1_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc1));']);
    eval(['ace.blDelta_rem_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc5));']);
    eval(['ace.blDelta_nostage_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc7));']);
    eval(['ace.blDelta_wake_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc0));']);
    
    eval(['ace.blSlowDelta_sws_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc3));']);
    eval(['ace.blSlowDelta_stg2_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc2));']);
    eval(['ace.blSlowDelta_stg1_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc1));']);
    eval(['ace.blSlowDelta_rem_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc5));']);
    eval(['ace.blSlowDelta_nostage_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc7));']);
    eval(['ace.blSlowDelta_wake_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc0));']);
    
    eval(['ace.blFastDelta_sws_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc3));']);
    eval(['ace.blFastDelta_stg2_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc2));']);
    eval(['ace.blFastDelta_stg1_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc1));']);
    eval(['ace.blFastDelta_rem_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc5));']);
    eval(['ace.blFastDelta_nostage_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc7));']);
    eval(['ace.blFastDelta_wake_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc0));']);
    
    eval(['ace.blAlpha_sws_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc3));']);
    eval(['ace.blAlpha_stg2_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc2));']);
    eval(['ace.blAlpha_stg1_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc1));']);
    eval(['ace.blAlpha_rem_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc5));']);
    eval(['ace.blAlpha_nostage_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc7));']);
    eval(['ace.blAlpha_wake_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc0));']);
    
    eval(['ace.blSigma_sws_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc3));']);
    eval(['ace.blSigma_stg2_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc2));']);
    eval(['ace.blSigma_stg1_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc1));']);
    eval(['ace.blSigma_rem_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc5));']);
    eval(['ace.blSigma_nostage_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc7));']);
    eval(['ace.blSigma_wake_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc0));']);
    
    eval(['ace.blTheta_sws_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc3));']);
    eval(['ace.blTheta_stg2_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc2));']);
    eval(['ace.blTheta_stg1_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc1));']);
    eval(['ace.blTheta_rem_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc5));']);
    eval(['ace.blTheta_nostage_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc7));']);
    eval(['ace.blTheta_wake_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc0));']);
    
    eval(['ace.blSpindleDensity_stg2_' Channels{1, Selection(i)}  '= bl_spindleDensity_2(1,i);']);
    eval(['ace.blSpindleDensity_sws_' Channels{1, Selection(i)}  '= bl_spindleDensity_3(1,i);']);
    
    eval(['ace.bl_Q1_SpindleDensity_stg2_' Channels{1, Selection(i)}  '= bl_spindleDensity_2_q1(1,i);']);
    eval(['ace.bl_Q2_SpindleDensity_stg2_' Channels{1, Selection(i)}  '= bl_spindleDensity_2_q2(1,i);']);
    eval(['ace.bl_Q3_SpindleDensity_stg2_' Channels{1, Selection(i)}  '= bl_spindleDensity_2_q3(1,i);']);
    eval(['ace.bl_Q4_SpindleDensity_stg2_' Channels{1, Selection(i)}  '= bl_spindleDensity_2_q4(1,i);']);
    
    eval(['ace.bl_Q1_SpindleDensity_stg3_' Channels{1, Selection(i)}  '= bl_spindleDensity_3_q1(1,i);']);
    eval(['ace.bl_Q2_SpindleDensity_stg3_' Channels{1, Selection(i)}  '= bl_spindleDensity_3_q2(1,i);']);
    eval(['ace.bl_Q3_SpindleDensity_stg3_' Channels{1, Selection(i)}  '= bl_spindleDensity_3_q3(1,i);']);
    eval(['ace.bl_Q4_SpindleDensity_stg3_' Channels{1, Selection(i)}  '= bl_spindleDensity_3_q4(1,i);']);
    
    eval(['ace.binSpindleDensity_stg2_Q1_' Channels{1, Selection(i)} '= sbj_hrb_spindles_densityPerMinute2_q1{1,i};']);
    eval(['ace.binSpindleDensity_stg2_Q2_' Channels{1, Selection(i)} '= sbj_hrb_spindles_densityPerMinute2_q2{1,i};']);
    eval(['ace.binSpindleDensity_stg2_Q3_' Channels{1, Selection(i)} '= sbj_hrb_spindles_densityPerMinute2_q3{1,i};']);
    eval(['ace.binSpindleDensity_stg2_Q4_' Channels{1, Selection(i)} '= sbj_hrb_spindles_densityPerMinute2_q4{1,i};']);
    
    eval(['ace.binSpindleDensity_stg3_Q1_' Channels{1, Selection(i)} '= sbj_hrb_spindles_densityPerMinute3_q1{1,i};']);
    eval(['ace.binSpindleDensity_stg3_Q2_' Channels{1, Selection(i)} '= sbj_hrb_spindles_densityPerMinute3_q2{1,i};']);
    eval(['ace.binSpindleDensity_stg3_Q3_' Channels{1, Selection(i)} '= sbj_hrb_spindles_densityPerMinute3_q3{1,i};']);
    eval(['ace.binSpindleDensity_stg3_Q4_' Channels{1, Selection(i)} '= sbj_hrb_spindles_densityPerMinute3_q4{1,i};']);
end
 
    eval(['ace.blHF_sws_' Channels{1, Selection(i)}  '= mean(hfpw(inc3));']);
    eval(['ace.blHF_stg2_' Channels{1, Selection(i)}  '= mean(hfpw(inc2));']);
    eval(['ace.blHF_stg1_' Channels{1, Selection(i)}  '= mean(hfpw(inc1));']);
    eval(['ace.blHF_rem_' Channels{1, Selection(i)}  '= mean(hfpw(inc5));']);
    eval(['ace.blHF_nostage_' Channels{1, Selection(i)}  '= mean(hfpw(inc7));']);
    eval(['ace.blHF_wake_' Channels{1, Selection(i)}  '= mean(hfpw(inc0));']);
    
    eval(['ace.blRR_sws_' Channels{1, Selection(i)}  '= mean(RRts(inc3));']);
    eval(['ace.blRR_stg2_' Channels{1, Selection(i)}  '= mean(RRts(inc2));']);
    eval(['ace.blRR_stg1_' Channels{1, Selection(i)}  '= mean(RRts(inc1));']);
    eval(['ace.blRR_rem_' Channels{1, Selection(i)}  '= mean(RRts(inc5));']);
    eval(['ace.blRR_nostage_' Channels{1, Selection(i)}  '= mean(RRts(inc7));']);
    eval(['ace.blRR_wake_' Channels{1, Selection(i)}  '= mean(RRts(inc0));']);
    
    
    eval(['ace.bL_Q1_HF_sws_' Channels{1, Selection(i)}  '= mean(hfpw(inc3_q1));']);
    eval(['ace.bl_Q1_HF_stg2_' Channels{1, Selection(i)}  '= mean(hfpw(inc2_q1));']);
    eval(['ace.bl_Q1_HF_stg1_' Channels{1, Selection(i)}  '= mean(hfpw(inc1_q1));']);
    eval(['ace.bl_Q1_HF_rem_' Channels{1, Selection(i)}  '= mean(hfpw(inc5_q1));']);
    eval(['ace.bl_Q1_HF_nostage_' Channels{1, Selection(i)}  '= mean(hfpw(inc7_q1));']);
    eval(['ace.bl_Q1_HF_wake_' Channels{1, Selection(i)}  '= mean(hfpw(inc0_q1));']);

    eval(['ace.bl_Q1_RR_sws_' Channels{1, Selection(i)}  '= mean(RRts(inc3_q1));']);
    eval(['ace.bl_Q1_RR_stg2_' Channels{1, Selection(i)}  '= mean(RRts(inc2_q1));']);
    eval(['ace.bl_Q1_RR_stg1_' Channels{1, Selection(i)}  '= mean(RRts(inc1_q1));']);
    eval(['ace.bl_Q1_RR_rem_' Channels{1, Selection(i)}  '= mean(RRts(inc5_q1));']);
    eval(['ace.bl_Q1_RR_nostage_' Channels{1, Selection(i)}  '= mean(RRts(inc7_q1));']);
    eval(['ace.bl_Q1_RR_wake_' Channels{1, Selection(i)}  '= mean(RRts(inc0_q1));']);
    
    eval(['ace.bL_Q2_HF_sws_' Channels{1, Selection(i)}  '= mean(hfpw(inc3_q2));']);
    eval(['ace.bl_Q2_HF_stg2_' Channels{1, Selection(i)}  '= mean(hfpw(inc2_q2));']);
    eval(['ace.bl_Q2_HF_stg1_' Channels{1, Selection(i)}  '= mean(hfpw(inc1_q2));']);
    eval(['ace.bl_Q2_HF_rem_' Channels{1, Selection(i)}  '= mean(hfpw(inc5_q2));']);
    eval(['ace.bl_Q2_HF_nostage_' Channels{1, Selection(i)}  '= mean(hfpw(inc7_q2));']);
    eval(['ace.bl_Q2_HF_wake_' Channels{1, Selection(i)}  '= mean(hfpw(inc0_q2));']);

    eval(['ace.bl_Q2_RR_sws_' Channels{1, Selection(i)}  '= mean(RRts(inc3_q2));']);
    eval(['ace.bl_Q2_RR_stg2_' Channels{1, Selection(i)}  '= mean(RRts(inc2_q2));']);
    eval(['ace.bl_Q2_RR_stg1_' Channels{1, Selection(i)}  '= mean(RRts(inc1_q2));']);
    eval(['ace.bl_Q2_RR_rem_' Channels{1, Selection(i)}  '= mean(RRts(inc5_q2));']);
    eval(['ace.bl_Q2_RR_nostage_' Channels{1, Selection(i)}  '= mean(RRts(inc7_q2));']);
    eval(['ace.bl_Q2_RR_wake_' Channels{1, Selection(i)}  '= mean(RRts(inc0_q2));']);

    eval(['ace.bL_Q3_HF_sws_' Channels{1, Selection(i)}  '= mean(hfpw(inc3_q3));']);
    eval(['ace.bl_Q3_HF_stg2_' Channels{1, Selection(i)}  '= mean(hfpw(inc2_q3));']);
    eval(['ace.bl_Q3_HF_stg1_' Channels{1, Selection(i)}  '= mean(hfpw(inc1_q3));']);
    eval(['ace.bl_Q3_HF_rem_' Channels{1, Selection(i)}  '= mean(hfpw(inc5_q3));']);
    eval(['ace.bl_Q3_HF_nostage_' Channels{1, Selection(i)}  '= mean(hfpw(inc7_q3));']);
    eval(['ace.bl_Q3_HF_wake_' Channels{1, Selection(i)}  '= mean(hfpw(inc0_q3));']);

    eval(['ace.bl_Q3_RR_sws_' Channels{1, Selection(i)}  '= mean(RRts(inc3_q3));']);
    eval(['ace.bl_Q3_RR_stg2_' Channels{1, Selection(i)}  '= mean(RRts(inc2_q3));']);
    eval(['ace.bl_Q3_RR_stg1_' Channels{1, Selection(i)}  '= mean(RRts(inc1_q3));']);
    eval(['ace.bl_Q3_RR_rem_' Channels{1, Selection(i)}  '= mean(RRts(inc5_q3));']);
    eval(['ace.bl_Q3_RR_nostage_' Channels{1, Selection(i)}  '= mean(RRts(inc7_q3));']);
    eval(['ace.bl_Q3_RR_wake_' Channels{1, Selection(i)}  '= mean(RRts(inc0_q3));']);

    eval(['ace.bL_Q4_HF_sws_' Channels{1, Selection(i)}  '= mean(hfpw(inc3_q4));']);
    eval(['ace.bl_Q4_HF_stg2_' Channels{1, Selection(i)}  '= mean(hfpw(inc2_q4));']);
    eval(['ace.bl_Q4_HF_stg1_' Channels{1, Selection(i)}  '= mean(hfpw(inc1_q4));']);
    eval(['ace.bl_Q4_HF_rem_' Channels{1, Selection(i)}  '= mean(hfpw(inc5_q4));']);
    eval(['ace.bl_Q4_HF_nostage_' Channels{1, Selection(i)}  '= mean(hfpw(inc7_q4));']);
    eval(['ace.bl_Q4_HF_wake_' Channels{1, Selection(i)}  '= mean(hfpw(inc0_q4));']);

    eval(['ace.bl_Q4_RR_sws_' Channels{1, Selection(i)}  '= mean(RRts(inc3_q4));']);
    eval(['ace.bl_Q4_RR_stg2_' Channels{1, Selection(i)}  '= mean(RRts(inc2_q4));']);
    eval(['ace.bl_Q4_RR_stg1_' Channels{1, Selection(i)}  '= mean(RRts(inc1_q4));']);
    eval(['ace.bl_Q4_RR_rem_' Channels{1, Selection(i)}  '= mean(RRts(inc5_q4));']);
    eval(['ace.bl_Q4_RR_nostage_' Channels{1, Selection(i)}  '= mean(RRts(inc7_q4));']);
    eval(['ace.bl_Q4_RR_wake_' Channels{1, Selection(i)}  '= mean(RRts(inc0_q4));']);

for i=1:length(Selection)
    eval(['ace.bl_Q1_Delta_sws_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc3_q1));']);
    eval(['ace.bl_Q1_Delta_stg2_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc2_q1));']);
    eval(['ace.bl_Q1_Delta_stg1_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc1_q1));']);
    eval(['ace.bl_Q1_Delta_rem_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc5_q1));']);
    eval(['ace.bl_Q1_Delta_nostage_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc7_q1));']);
    eval(['ace.bl_Q1_Delta_wake_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc0_q1));']);
    
    eval(['ace.bl_Q1_SlowDelta_sws_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc3_q1));']);
    eval(['ace.bl_Q1_SlowDelta_stg2_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc2_q1));']);
    eval(['ace.bl_Q1_SlowDelta_stg1_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc1_q1));']);
    eval(['ace.bl_Q1_SlowDelta_rem_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc5_q1));']);
    eval(['ace.bl_Q1_SlowDelta_nostage_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc7_q1));']);
    eval(['ace.bl_Q1_SlowDelta_wake_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc0_q1));']);
    
    eval(['ace.bl_Q1_FastDelta_sws_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc3_q1));']);
    eval(['ace.bl_Q1_FastDelta_stg2_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc2_q1));']);
    eval(['ace.bl_Q1_FastDelta_stg1_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc1_q1));']);
    eval(['ace.bl_Q1_FastDelta_rem_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc5_q1));']);
    eval(['ace.bl_Q1_FastDelta_nostage_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc7_q1));']);
    eval(['ace.bl_Q1_FastDelta_wake_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc0_q1));']);
    
    eval(['ace.bl_Q1_Alpha_sws_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc3_q1));']);
    eval(['ace.bl_Q1_Alpha_stg2_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc2_q1));']);
    eval(['ace.bl_Q1_Alpha_stg1_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc1_q1));']);
    eval(['ace.bl_Q1_Alpha_rem_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc5_q1));']);
    eval(['ace.bl_Q1_Alpha_nostage_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc7_q1));']);
    eval(['ace.bl_Q1_Alpha_wake_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc0_q1));']);
    
    eval(['ace.bl_Q1_Sigma_sws_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc3_q1));']);
    eval(['ace.bl_Q1_Sigma_stg2_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc2_q1));']);
    eval(['ace.bl_Q1_Sigma_stg1_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc1_q1));']);
    eval(['ace.bl_Q1_Sigma_rem_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc5_q1));']);
    eval(['ace.bl_Q1_Sigma_nostage_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc7_q1));']);
    eval(['ace.bl_Q1_Sigma_wake_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc0_q1));']);
    
    eval(['ace.bl_Q1_Theta_sws_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc3_q1));']);
    eval(['ace.bl_Q1_Theta_stg2_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc2_q1));']);
    eval(['ace.bl_Q1_Theta_stg1_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc1_q1));']);
    eval(['ace.bl_Q1_Theta_rem_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc5_q1));']);
    eval(['ace.bl_Q1_Theta_nostage_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc7_q1));']);
    eval(['ace.bl_Q1_Theta_wake_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc0_q1));']);
   
end

for i=1:length(Selection)
    eval(['ace.bl_Q2_Delta_sws_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc3_q2));']);
    eval(['ace.bl_Q2_Delta_stg2_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc2_q2));']);
    eval(['ace.bl_Q2_Delta_stg1_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc1_q2));']);
    eval(['ace.bl_Q2_Delta_rem_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc5_q2));']);
    eval(['ace.bl_Q2_Delta_nostage_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc7_q2));']);
    eval(['ace.bl_Q2_Delta_wake_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc0_q2));']);
    
    eval(['ace.bl_Q2_SlowDelta_sws_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc3_q2));']);
    eval(['ace.bl_Q2_SlowDelta_stg2_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc2_q2));']);
    eval(['ace.bl_Q2_SlowDelta_stg1_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc1_q2));']);
    eval(['ace.bl_Q2_SlowDelta_rem_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc5_q2));']);
    eval(['ace.bl_Q2_SlowDelta_nostage_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc7_q2));']);
    eval(['ace.bl_Q2_SlowDelta_wake_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc0_q2));']);
    
    eval(['ace.bl_Q2_FastDelta_sws_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc3_q2));']);
    eval(['ace.bl_Q2_FastDelta_stg2_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc2_q2));']);
    eval(['ace.bl_Q2_FastDelta_stg1_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc1_q2));']);
    eval(['ace.bl_Q2_FastDelta_rem_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc5_q2));']);
    eval(['ace.bl_Q2_FastDelta_nostage_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc7_q2));']);
    eval(['ace.bl_Q2_FastDelta_wake_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc0_q2));']);
    
    eval(['ace.bl_Q2_Alpha_sws_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc3_q2));']);
    eval(['ace.bl_Q2_Alpha_stg2_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc2_q2));']);
    eval(['ace.bl_Q2_Alpha_stg1_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc1_q2));']);
    eval(['ace.bl_Q2_Alpha_rem_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc5_q2));']);
    eval(['ace.bl_Q2_Alpha_nostage_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc7_q2));']);
    eval(['ace.bl_Q2_Alpha_wake_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc0_q2));']);
    
    eval(['ace.bl_Q2_Sigma_sws_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc3_q2));']);
    eval(['ace.bl_Q2_Sigma_stg2_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc2_q2));']);
    eval(['ace.bl_Q2_Sigma_stg1_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc1_q2));']);
    eval(['ace.bl_Q2_Sigma_rem_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc5_q2));']);
    eval(['ace.bl_Q2_Sigma_nostage_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc7_q2));']);
    eval(['ace.bl_Q2_Sigma_wake_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc0_q2));']);
    
    eval(['ace.bl_Q2_Theta_sws_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc3_q2));']);
    eval(['ace.bl_Q2_Theta_stg2_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc2_q2));']);
    eval(['ace.bl_Q2_Theta_stg1_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc1_q2));']);
    eval(['ace.bl_Q2_Theta_rem_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc5_q2));']);
    eval(['ace.bl_Q2_Theta_nostage_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc7_q2));']);
    eval(['ace.bl_Q2_Theta_wake_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc0_q2));']);
   
end

for i=1:length(Selection)
    eval(['ace.bl_Q3_Delta_sws_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc3_q3));']);
    eval(['ace.bl_Q3_Delta_stg2_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc2_q3));']);
    eval(['ace.bl_Q3_Delta_stg1_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc1_q3));']);
    eval(['ace.bl_Q3_Delta_rem_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc5_q3));']);
    eval(['ace.bl_Q3_Delta_nostage_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc7_q3));']);
    eval(['ace.bl_Q3_Delta_wake_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc0_q3));']);
    
    eval(['ace.bl_Q3_SlowDelta_sws_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc3_q3));']);
    eval(['ace.bl_Q3_SlowDelta_stg2_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc2_q3));']);
    eval(['ace.bl_Q3_SlowDelta_stg1_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc1_q3));']);
    eval(['ace.bl_Q3_SlowDelta_rem_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc5_q3));']);
    eval(['ace.bl_Q3_SlowDelta_nostage_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc7_q3));']);
    eval(['ace.bl_Q3_SlowDelta_wake_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc0_q3));']);
    
    eval(['ace.bl_Q3_FastDelta_sws_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc3_q3));']);
    eval(['ace.bl_Q3_FastDelta_stg2_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc2_q3));']);
    eval(['ace.bl_Q3_FastDelta_stg1_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc1_q3));']);
    eval(['ace.bl_Q3_FastDelta_rem_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc5_q3));']);
    eval(['ace.bl_Q3_FastDelta_nostage_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc7_q3));']);
    eval(['ace.bl_Q3_FastDelta_wake_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc0_q3));']);
    
    eval(['ace.bl_Q3_Alpha_sws_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc3_q3));']);
    eval(['ace.bl_Q3_Alpha_stg2_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc2_q3));']);
    eval(['ace.bl_Q3_Alpha_stg1_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc1_q3));']);
    eval(['ace.bl_Q3_Alpha_rem_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc5_q3));']);
    eval(['ace.bl_Q3_Alpha_nostage_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc7_q3));']);
    eval(['ace.bl_Q3_Alpha_wake_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc0_q3));']);
    
    eval(['ace.bl_Q3_Sigma_sws_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc3_q3));']);
    eval(['ace.bl_Q3_Sigma_stg2_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc2_q3));']);
    eval(['ace.bl_Q3_Sigma_stg1_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc1_q3));']);
    eval(['ace.bl_Q3_Sigma_rem_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc5_q3));']);
    eval(['ace.bl_Q3_Sigma_nostage_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc7_q3));']);
    eval(['ace.bl_Q3_Sigma_wake_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc0_q3));']);
    
    eval(['ace.bl_Q3_Theta_sws_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc3_q3));']);
    eval(['ace.bl_Q3_Theta_stg2_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc2_q3));']);
    eval(['ace.bl_Q3_Theta_stg1_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc1_q3));']);
    eval(['ace.bl_Q3_Theta_rem_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc5_q3));']);
    eval(['ace.bl_Q3_Theta_nostage_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc7_q3));']);
    eval(['ace.bl_Q3_Theta_wake_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc0_q3));']);
   
end

for i=1:length(Selection)
    eval(['ace.bl_Q4_Delta_sws_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc3_q4));']);
    eval(['ace.bl_Q4_Delta_stg2_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc2_q4));']);
    eval(['ace.bl_Q4_Delta_stg1_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc1_q4));']);
    eval(['ace.bl_Q4_Delta_rem_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc5_q4));']);
    eval(['ace.bl_Q4_Delta_nostage_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc7_q4));']);
    eval(['ace.bl_Q4_Delta_wake_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc0_q4));']);
    
    eval(['ace.bl_Q4_SlowDelta_sws_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc3_q4));']);
    eval(['ace.bl_Q4_SlowDelta_stg2_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc2_q4));']);
    eval(['ace.bl_Q4_SlowDelta_stg1_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc1_q4));']);
    eval(['ace.bl_Q4_SlowDelta_rem_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc5_q4));']);
    eval(['ace.bl_Q4_SlowDelta_nostage_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc7_q4));']);
    eval(['ace.bl_Q4_SlowDelta_wake_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc0_q4));']);
    
    eval(['ace.bl_Q4_FastDelta_sws_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc3_q4));']);
    eval(['ace.bl_Q4_FastDelta_stg2_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc2_q4));']);
    eval(['ace.bl_Q4_FastDelta_stg1_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc1_q4));']);
    eval(['ace.bl_Q4_FastDelta_rem_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc5_q4));']);
    eval(['ace.bl_Q4_FastDelta_nostage_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc7_q4));']);
    eval(['ace.bl_Q4_FastDelta_wake_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc0_q4));']);
    
    eval(['ace.bl_Q4_Alpha_sws_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc3_q4));']);
    eval(['ace.bl_Q4_Alpha_stg2_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc2_q4));']);
    eval(['ace.bl_Q4_Alpha_stg1_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc1_q4));']);
    eval(['ace.bl_Q4_Alpha_rem_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc5_q4));']);
    eval(['ace.bl_Q4_Alpha_nostage_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc7_q4));']);
    eval(['ace.bl_Q4_Alpha_wake_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc0_q4));']);
    
    eval(['ace.bl_Q4_Sigma_sws_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc3_q4));']);
    eval(['ace.bl_Q4_Sigma_stg2_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc2_q4));']);
    eval(['ace.bl_Q4_Sigma_stg1_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc1_q4));']);
    eval(['ace.bl_Q4_Sigma_rem_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc5_q4));']);
    eval(['ace.bl_Q4_Sigma_nostage_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc7_q4));']);
    eval(['ace.bl_Q4_Sigma_wake_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc0_q4));']);
    
    eval(['ace.bl_Q4_Theta_sws_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc3_q4));']);
    eval(['ace.bl_Q4_Theta_stg2_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc2_q4));']);
    eval(['ace.bl_Q4_Theta_stg1_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc1_q4));']);
    eval(['ace.bl_Q4_Theta_rem_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc5_q4));']);
    eval(['ace.bl_Q4_Theta_nostage_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc7_q4));']);
    eval(['ace.bl_Q4_Theta_wake_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc0_q4));']);
   
end
    ace.HRB_density_wake_quartile = zeros(1,q);
    ace.HRB_density_stg1_quartile = zeros(1,q);
    ace.HRB_density_stg2_quartile = zeros(1,q);
    ace.HRB_density_sws_quartile = zeros(1,q);
    ace.HRB_density_rem_quartile = zeros(1,q);
    
    for i = 1: q
        if ace.minutes_q_wake(i) ~=0
            ace.HRB_density_wake_quartile(i) = sum(ace.HRB_wake_quartile==i) / ace.minutes_q_wake(i);
        end
        if ace.minutes_q_stg1(i) ~=0
            ace.HRB_density_stg1_quartile(i) = sum(ace.HRB_stg1_quartile==i) / ace.minutes_q_stg1(i);
        end
        if ace.minutes_q_stg2(i) ~=0
            ace.HRB_density_stg2_quartile(i) = sum(ace.HRB_stg2_quartile==i) / ace.minutes_q_stg2(i);
        end
        if ace.minutes_q_sws(i) ~=0
            ace.HRB_density_sws_quartile(i) = sum(ace.HRB_sws_quartile==i) / ace.minutes_q_sws(i);
        end
        if ace.minutes_q_rem(i) ~=0
            ace.HRB_density_rem_quartile(i) = sum(ace.HRB_rem_quartile==i) / ace.minutes_q_rem(i);
        end
    end
    ace.fs=fs;
    ace.quartile_minute=ql/(fs*60);
    ace.filename=FileName;
    save([PathName 'ace' FileName(1:end-4)],'ace')
end


function [X,Channels,fs]=edf2mat_samefss(fileloc)
% [FileName,PathName] = uigetfile('/bazhlab/naji/home/EDFs_ACH_500Hz/*.edf','Select the edf data file');
% load([PathName FileName]);
fid=fopen(fileloc);% in format of [PathName FileName]
a=fread(fid,236,'*char');
ndr=fread(fid,8,'*char');
ndr=str2double(ndr'); %number of data records in sec
a=fread(fid,8,'*char');
drdur=str2double(a'); %duration of each data record in sec
ns=fread(fid,4,'*char'); ns=ns'; ns=str2double(ns);% number of signal channels
Channels=cell(1,ns);
for i=1:ns
    C=fread(fid,16,'*char');C=C';
    Channels{i}=C(find(isspace(C)==0));
end
fread(fid,ns*80,'*char'); % channel transducer type can be extracted
fread(fid,ns*8,'*char'); %channel physical dimension can be extracted
phmn=zeros(1,ns); phmx=phmn;dmn=phmn;dmx=dmn;
for i=1:ns
    pm=fread(fid,8,'*char');pm=pm';
    phmn(i)=str2double(pm);
end                         %phys min
for i=1:ns
    pm=fread(fid,8,'*char'); pm=pm';
    phmx(i)=str2double(pm);%phys max
end
for i=1:ns
    dm=fread(fid,8,'*char');dm=dm'; 
    dmn(i)=str2double(dm);
end                         %dig min
for i=1:ns
    dx=fread(fid,8,'*char'); dx=dx';
    dmx(i)=str2double(dx);
end                         %dig max
scalefac=(phmx-phmn)./(dmx-dmn);
dc=phmx-scalefac.*dmx;

fread(fid,ns*80,'*char'); % prefilters
nr=zeros(1,ns);
for i=1:ns
    nrc=fread(fid,8,'*char'); nrc=nrc';
    nr(i)=str2double(nrc); %number of samples in each data record
end
if sum(ismember(nr,nr(1)))==length(nr)
    fs=nr(1);
else
    sprintf('Data cant be stored in a single matrix')
end
fread(fid,ns*32,'*char');
ch_fs=nr/drdur;
if mean(nr)==nr(1) && mean(ch_fs)==ch_fs(1)
X=zeros(ns,nr(1)*ndr);
% for i=1:ns
%     X{i,1}=zeros(1,nr(i)*ndr);
% end
fs=ch_fs(1);
end
spins={'\\','|','/','-'};
    reverseStr = 'Reading EDF file ... ';
for i=1:ndr
    for j=1:ns
        s=fread(fid,nr(j),'int16').*scalefac(j)+dc(j);s=s';
        X(j,(i-1)*nr(j)+1:i*nr(j))=s;
    end
    
%     si=mod(i,4); 
%     if si==0
%         si=4;
%     end
%     msg = (spins{si});
%     fprintf([reverseStr, msg]);
%     reverseStr = repmat(sprintf('\b'), 1, 1);
    
end
fprintf('\n');
fclose(fid);
end

function [sbj_hrb_ind,sbj_hrb,sbj_hrbHFPW]=myhrbwindowedd(in,bnd,fs,hrbwin,RRts,hfpw,m,RR_lowerlim, RR_upperlim)
sbj_hrb_ind=[];
    sbj_hrb=[];
    sbj_hrbHFPW=[];
if ~isempty(in)
    
    for i=1:length(in)
        pmin=myHRB_finder(RRts(bnd(in(i),1):bnd(in(i),2)),fs,m, RR_lowerlim, RR_upperlim);
        for j=1:length(pmin)
            sbj_hrb_ind=[sbj_hrb_ind;pmin(j)+bnd(in(i),1)-1];
        end
    end
    if ~isempty(sbj_hrb_ind)
        rj=find(sbj_hrb_ind/fs<(hrbwin/2) | (length(RRts)-sbj_hrb_ind)/fs<(hrbwin/2));
        sbj_hrb_ind(rj)=[];
    end
    for i=1:length(sbj_hrb_ind)
        sbj_hrb=[sbj_hrb;RRts(sbj_hrb_ind(i)-floor(hrbwin/2*fs)+1:sbj_hrb_ind(i)+floor(hrbwin/2*fs))];
        sbj_hrbHFPW=[sbj_hrbHFPW;hfpw(sbj_hrb_ind(i)-floor(hrbwin/2*fs)+1:sbj_hrb_ind(i)+floor(hrbwin/2*fs))];
    end
end
end

function out = myHRB_spindle_intersect(hrb_i, spi_i, hrbwin,fs)
out = cell(1,length(spi_i));
for ch=1:length(spi_i)
    out{1,ch} = zeros(length(hrb_i),4);
    for ix=1:length(hrb_i)
        rr_seg = hrb_i(ix)- (fs*hrbwin/2+1):hrb_i(ix)+ fs*hrbwin/2;
        [~,ia,~] = intersect(rr_seg, floor(mean(spi_i{1, ch}'))); 
        for ss=1:length(ia)
            wi = 1+floor(ia(ss)/length(rr_seg)*4);
            out{1,ch}(ix,wi) = out{1,ch}(ix,wi) + 1;
        end
    end          
end
end

function out = myblspindle(in,bnd,hrbind,hrbwin,fs, spi_i)
out = zeros(1,length(spi_i));
inc=[];
for i=1:length(in)
    inc=[inc bnd(in(i),1):bnd(in(i),2)];
end
rj=[];
for i=1:length(hrbind)
    rj=[rj hrbind(i)-floor(hrbwin/2*fs)+1:hrbind(i)+floor(hrbwin/2*fs)];
end
inc(find(ismember(inc,rj)==1))=[];
for ch=1:length(spi_i)
    [~,ia,~] = intersect(inc, floor(mean(spi_i{1, ch}'))); 
    out(1,ch) = length(ia)/(length(inc)/fs)*60;
end
end

function inc=myblactivityy(in,bnd,hrbind,hrbwin,fs)
inc=[];
    for i=1:length(in)
        inc=[inc bnd(in(i),1):bnd(in(i),2)];
    end
    rj=[];
    for i=1:length(hrbind)
        rj=[rj hrbind(i)-floor(hrbwin/2*fs)+1:hrbind(i)+floor(hrbwin/2*fs)];
    end
    inc(find(ismember(inc,rj)==1))=[];
end

function hrb_EEG=myhrbwindowedd_EEG(X,hrbwin,fs,oc,sbj_hrb_ind)
hrb_EEG=cell(size(X,1),1);
for j=1:size(X,1)
    hrb_EEG{j,1}=zeros(length(sbj_hrb_ind),hrbwin*fs);
    
    for i=1:length(sbj_hrb_ind)
        hrb_EEG{j,1}(i,:)=X(j,sbj_hrb_ind(i)-floor(hrbwin/2*fs)+1:sbj_hrb_ind(i)+floor(hrbwin/2*fs));
    end
    if oc==1
    mm=max(hrb_EEG{j,1}');
    rj=find(mm>(mean(mm)+3*std(mm)));
    hrb_EEG{j,1}(rj,:)=[];
    end
end
    
end

function hrb_EEG_bins=myhrbbinnedd_EEG(x,avbin,fs,oc)
if iscell(x)
hrb_EEG_bins=cell(size(x,1),1);
for j=1:size(x,1)
    hrb_EEG_bins{j,1}=zeros(size(x{j,1},1),size(x{j,1},2)/(fs*avbin));
    
    for i=1:size(x{j,1},1)
        for ii=1:size(x{j,1},2)/(fs*avbin)
        hrb_EEG_bins{j,1}(i,ii)=mean(x{j,1}(i,(ii-1)*avbin*fs+1:ii*avbin*fs));
        end
    end
    if oc==1
        mm=max(hrb_EEG_bins{j,1}');
        rj=find(mm>(mean(mm)+3*std(mm)));
        hrb_EEG_bins{j,1}(rj,:)=[];
    end
end
else
    hrb_EEG_bins=[];
    for i=1:size(x,1)
        for ii=1:size(x,2)/(fs*avbin)
            hrb_EEG_bins(i,ii)=mean(x(i,(ii-1)*avbin*fs+1:ii*avbin*fs));
        end
    end
    if oc==1
        mm=max(hrb_EEG_bins');
        rj=find(mm>(mean(mm)+3*std(mm)));
        hrb_EEG_bins(rj,:)=[];
    end
end
 
end

function mn = mean_rj_olyr(x,idx)
if length(idx) == 0
    mn = NaN;
else
    if size(x(idx,:),1) == 0
        mn = NaN;

    elseif size(x(idx,:),1) == 1
        mn = x(idx,:);
    else
        x = x(idx,:);
        mm=max(x');
        rjc=find(mm>(mean(mm)+3*std(mm)));
        x(rjc,:)=[];
        mn = mean(x, 1);
    end
end
end

function pmin=myHRB_finder(samples,fs,m, RR_lowerlim, RR_upperlim)
% m=1.25;
allmin=find((samples(1:end-2)-samples(2:end-1))>=0 & (samples(2:end-1)-samples(3:end))<=0);
th=mean(samples(allmin))-m*std(samples(allmin)); %threshold
pmin=find((samples(1:end-2)-samples(2:end-1))>=0 & (samples(2:end-1)-samples(3:end))<=0 & samples(1:end-2)<=th);
jj=1; pmnsp=[]; rjcts=[]; cmpmn=[];

            for ii=1:length(pmin)-1
                if (pmin(ii+1)-pmin(ii))>10*fs
                    jj=jj+1;
                    pmnsp=[];
                end
                if (pmin(ii+1)-pmin(ii))<10*fs
                    pmnsp=[pmnsp pmin(ii) pmin(ii+1)];
                    cmpmn{jj,1}=pmnsp; 
                end
            end
            [rcmp,~]=size(cmpmn);
            for ii=1:rcmp
                if ~isempty(cmpmn{ii,1})
                    tempv=cmpmn{ii,1};
                    [~,kkn]=min(samples(tempv));
                    tempv(kkn)=[];
                    rjcts=[rjcts tempv];
                end
            end
            rjsmp=[];
            for ii=1:length(rjcts)
                rjsmp=[rjsmp find(pmin== rjcts(ii))];
            end
            pmin(rjsmp)=[];
%             falseR=find(samples(pmin)<0.55); pmin(falseR)=[];
            
            rj = [];
            for ii=1:length(pmin)
                win = samples(max(1,pmin - 10*fs): min(pmin+10*fs , length(samples)));
                if min(win) < RR_lowerlim || max(win)>RR_upperlim
                    rj=[rj;ii];
                end
            end
            pmin(rj)=[];
end

function [detection,bmn_s,bnd_s] = Mywamsley_spindle_detection(C3_N2,C3,fs)
% WAMSLEY Detect sleep spindles using the Wamsley algorithm.
% E. Wamsley et al. "Reduced Sleep Spindles and Spindle Coherence in 
% Schizophrenia: Mechanisms of Impaired Memory Consolidation", 
% Biol Psychiatry 71, 2012, pp 154-161
%
% First input is a structure containing EEG from C3-M2 during all NREM
% (S2) episodes, each cell in the structure is a continuous NREM
% segment [muV].
% Second input is EEG from C3-M2 from the entire night [muV].
% Last input is the sampling frequency [Hz].
% Syntax: detection = wamsley_spindle_detection(C3_N2,C3,fs)
%
% Adopted from Wamsley by M.Sc. Sabrina Lyngbye Wendt, July 2013

signalmean = threshold_wamsley(C3_N2,fs);
[detection,bmn_s,bnd_s] = wamsley(C3,fs,signalmean);

    function signalmean = threshold_wamsley(C3nrem2,fs)
        % THRESHOLD_WAMSLEY Calculates the amplitude criteria for spindle detection
        % using the wamsley method.
        % Input is a structure where each entry contains a continuous segment of
        % EEG data recorded at C3-M2 during S2. The sampling frequency is the
        % final input.
        
        %% Define parameters for the wavelet analysis
        fb = 13.5;
        fc = 0.5;
        scale = fs*fc/fb; %9.48 for 256Hz; %3.7 for 100 Hz
        
        Ltotal = 0;
        for k = 1:length(C3nrem2)
            signal = C3nrem2{k}; L = length(signal);
            %% Perform wavelet transformation
            EEGWave = cwt(signal,scale,['cmor' num2str(fb) '-' num2str(fc)]);
            EEGData = real(EEGWave.^2);
            
            %% Take Moving Average
            EEGData = EEGData.^2;
            window = ones(ceil(fs/10),1)/ceil(fs/10); % create 100ms window to convolve with
            EEGData2 = filter(window,1,EEGData); % take the moving average using the above window
            MA(Ltotal+1:Ltotal+L) = EEGData2;
            Ltotal = Ltotal+L;
        end
        
        %% Determine amplitude threshold
        signalmean = mean(MA);
    end

    function [detection,bmn_s,bnd_s] = wamsley(C3,fs,signalmean)
        % WAMSLEY Detect sleep spindles in EEG given the amplitude criteria.
        % Input is the EEG signal we wish to detect spindles in, the sampling
        % frequency and the amplitude criteria.
        % Output is a vector containing the detection of spindles.
        
        %% Define parameters for the wavelet analysis
        fb = 13.5;
        fc = 0.5;
        scale = fs*fc/fb; %9.48 for 256Hz; %3.7 for 100 Hz
        
        EEGWave = cwt(C3,scale,['cmor' num2str(fb) '-' num2str(fc)]);
        EEGData = real(EEGWave.^2);
        
        %% Take Moving Average
        EEGData = EEGData.^2;
        window = ones(ceil(fs/10),1)/ceil(fs/10); % create 100ms window to convolve with
        EEGData2 = filter(window,1,EEGData); % take the moving average using the above window
        
        %% Determine amplitude threshold
        threshold = signalmean.*4; % defines the threshold was 4.5
        
        %% Find Peaks in the MS Signal
        current_data=EEGData2;
        
        over=current_data>threshold; % Mark all points over threshold as '1'
        detection = zeros(length(current_data),1);
        detection(over) = 1;
        [begins,ends] = find_spindles(detection);
        [detection,begins,ends] = maximum_duration(detection,begins,ends,3,fs);
        [detection,begins_03,ends_03] = minimum_duration(detection,begins,ends,0.4,fs);
        
        locs_03=(zeros(1,length(current_data)))';  % Create a vector of zeros the length of the MS signal
        for i=1:((length(current_data))-ceil(fs*0.4));  % for the length of the signal, if the sum of 30 concurrent points = Fs*0.3, mark a spindle
            if sum(over(i:(i+(ceil(fs*0.4)-1))))==ceil(fs*0.4);
                locs_03(i,1)=1;
            end
        end
        
        spin_03=zeros((length(locs_03)),1);  % only mark a spindle in vector 'spin' at the end of a 300ms duration peak
        for i=1:length(locs_03);
            if locs_03(i,1)==1 && locs_03(i+1,1)==0;
                spin_03(i,1)=1;
            end
        end
        
%         for i=513:length(spin_03);%201-->513  % for every spindle marked in 'spin', delete the spindle if there is also a spindle within the second preceeding it
%             if spin_03(i,1)==1 && sum(spin_03((i-fs):(i-1)))>0;
%                 spin_03(i,1)=0;
%                 idx = find(i>=begins_03 & i<=ends_03);
%                 if isempty(idx) == 0
%                     detection(begins_03(idx):ends_03(idx)) = 0;
%                 else
%                     error('Did not find spindle beginning and ending around a spin point')
%                 end
%             end
%         end
        
        t=find((detection(2:end)-detection(1:end-1))~=0);
        t=[0 t' length(detection)];
        clear bnd_s
        for i=1:length(t)-1
          bnd_s(i,:)=[t(i)+1 t(i+1)];
        end
        bmn_s=(bnd_s(:,2)-bnd_s(:,1)+1)./fs;
        bmn_s(find(detection(bnd_s(:,1))==0),:)=[];
        bnd_s(find(detection(bnd_s(:,1))==0),:)=[];
    end

    function [begins, ends] = find_spindles(bv)
        % FIND_SPINDLES - find start and end index' of spindles.
        % Input is a binary vector bv containing ones where spindles are detected.
        % Output is vectors containing the index' of spindle beginnings and ends
        % (first sample of spindle and last sample of spindle, respectively).
        
        sise = size(bv);
        E = bv(2:end)-bv(1:end-1); % Find start and end of intervals with spindles
        
        begins = find(E==1)+1;
        if bv(1) == 1
            if sise(1) > 1
                begins = [1; begins];
            elseif sise(2) > 1
                begins = [1 begins];
            else
                error('The input signal is not one dimensional')
            end
        elseif numel(begins) == 0 && bv(1) == 0
            begins = NaN;
        end
        
        ends = find(E==-1);
        if bv(end) == 1
            if sise(1) > 1
                ends = [ends; length(bv)];
            elseif sise(2) > 1
                ends = [ends length(bv)];
            else
                error('The input signal is not one dimensional')
            end
        elseif numel(ends) == 0 && bv(end) == 0
            ends = NaN;
        end
    end

    function [bv,begins,ends] = maximum_duration(bv,begins,ends,max_dur,fs)
        % MAXIMUM_DURATION - checks the sample duration of the spindles.
        % Input is a vector containing ones in the interval where the spindle is
        % and indexs describing the start and end of the spindle. The last two
        % inputs are the maximum duration given in seconds and the sampling
        % frequency given in Hz.
        % Output is a vector containing ones in the interval where the spindle with
        % duration shorter than or equal to the maximum duration is and indexs
        % describing the start and end of the spindle.
        
        duration_samples = ends-begins+1;
        for k = 1:length(begins)
            if duration_samples(k) > max_dur*fs
                bv(begins(k):ends(k)) = 0;
                begins(k) = 0;
                ends(k) = 0;
            end
        end
        begins = begins(begins~=0);
        ends = ends(ends~=0);
    end

    function [bv,begins,ends] = minimum_duration(bv,begins,ends,min_dur,fs)
        % MINIMUM_DURATION - checks the sample duration of the spindles.
        % Input is a vector containing ones in the interval where the spindle is
        % and indexs describing the start and end of the spindle. The last two
        % inputs are the minimum duration given in seconds and the sampling
        % frequency given in Hz.
        % Output is a vector containing ones in the interval where the spindle with
        % duration longer than or equal to the minimum duration is and indexs
        % describing the start and end of the spindle.
        
        duration_samples = ends-begins+1;
        for k = 1:length(begins)
            if duration_samples(k) < ceil(min_dur*fs)
                bv(begins(k):ends(k)) = 0;
                begins(k) = 0;
                ends(k) = 0;
            end
        end
        begins = begins(begins~=0);
        ends = ends(ends~=0);
    end

end
