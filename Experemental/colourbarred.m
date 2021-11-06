clear all;
fname = 'Audi_A1_Driving_Away_30KPH.wav';
CPI=0.6;
OLF =0.5;

pfa = 10^-50;%desired Probabability of false alarm
RefWindow =50;%size of data cells
gcell=3;% number of guardcells
MakeX(fname,CPI,OLF,pfa,RefWindow,gcell)

function [] = MakeX(fname, CPI,OLF,pfa1,RefWindow1,gcell1)

%INPUT VARIABLES
% fname = 'file name
% CPI=Coherent processing interval
% OLF = overlap factor
% pfa1 = desired Probabability of false alarm
% RefWindow1 = size of one side of data cells 
% gcell1= number of guardcells

%constants
fc = 2590e6; % (Hz) I DON't know, a guess is the baseband thingy
c = 299e6;% speed of light
Maxspeed=100;



%load in file
[X,fs]=audioread(fname,'native');
x = -X(:,2);
N = length(x);% $ of samples
fs;
%Preparing the window
Ns=2^nextpow2(CPI*fs);
CPI=Ns/fs;
slide= floor(Ns*(1-OLF));
Nf=floor(((N-Ns)/slide)+1);
win = hamming(Ns);%not sure how i would go about changing the type outside the function 


%windowing
for k=0:(Nf-1)
    ywin=double(x(1+(k*slide):(k*slide)+Ns)).*win;
    ywin = ywin - mean(ywin);
    Y=(fft(ywin));%
   
    stfft(:,k+1)=fftshift(Y);
end
%fax=0:fs/Ns:2*fs;

fax=(-Ns/2:1:(Ns/2-1))*fs/Ns;
t =(Ns/2:slide:Ns/2+(Nf-1)*slide)/fs;    

%making vs speed
lambda = c/fc;
speed_m_per_sec = fax*lambda/2;
speed_km_per_hr = speed_m_per_sec*(60*60/1000);
speed_km_per_hr_Idx = find((speed_km_per_hr <= Maxspeed) & (speed_km_per_hr >= 0));

SpeedVectorOfInterest = speed_km_per_hr(speed_km_per_hr_Idx);
S_OfInterest = stfft(speed_km_per_hr_Idx, :);

%plotting spectogram
specm=max(abs(S_OfInterest));
spec=abs(S_OfInterest)./specm;
spec=20*log10(spec);
figure();
clims = [-40 0];
imagesc(t,SpeedVectorOfInterest,spec,clims);
colorbar;
colormap('jet');
xlabel("time (s)");
ylabel("speed (km/h)");
hold on

%CFAR
count1 =zeros(Nf,1);%number of detections per iteration
err=[];%pfa error if no signal otherwise the signal counts as false alarms so doesn't work
signal = transpose(stfft);%
detect=zeros(length(t),length(speed_km_per_hr_Idx)-1);
threshold=zeros(length(t),length(speed_km_per_hr_Idx)-1);
for e=1:Nf
    y_in=signal(e,:);
    [count1(e),detect(e,:),treshold(e,:)]=CACFAR(pfa1,RefWindow1,gcell1,y_in,speed_km_per_hr_Idx);
end

%marrking targets
total_fa=zeros(Nf,length(signal(1,:)));

for ix=1:Nf
    for xi=1:length(SpeedVectorOfInterest)-1
        temp=detect(ix,xi);
        if (~temp==0)
            total_fa(ix,xi)=SpeedVectorOfInterest(temp);
        end    
    end

end


total_fa(total_fa<=0)= nan;

plot(t,total_fa,'kx','MarkerSize',12);
hold off;


%%PLOT SPEED
%speedvtime(t, total_fa); 
end

function [count,intercept,thold]=CACFAR(pfa,RefWindow,gcell,y_in,indx) 
%CFAR PROCESSING
    %INTERNAL PARAMS
    sidecell=gcell+RefWindow;
    len=length(y_in);
    y_afterpow = abs(y_in).^2;
    thold=zeros(1,max(indx)-min(indx));
    count1=1;
    alpha=(pfa^(-1/(RefWindow*2))-1);
    intercept=zeros(1,max(indx)-min(indx));
    %PROCESS
    count=1;
    for i=indx
        lag=y_afterpow(i-sidecell:i-1-gcell);
        lead=y_afterpow(i+1+gcell:i+sidecell);
        mlag=sum(lag);
        mlead=sum(lead);
        gi=(mlead+mlag);%max(mlead,mlag);%change between CA-Cfar and OSGO-CFAR
        thold(i)=alpha*gi;
        if (y_afterpow(i) >= thold(i))
            intercept(count1)=count;
            count1=count1+1;
        end
count=count+1;
    end
end

function [] = speedvtime(time, points)
store = nanmean(transpose(points()));
figure(2);
hold on
plot(time,store,'r','MarkerSize',3);
xlabel("time (s)");
ylabel("speed (km/h)");
title("Speed vs. Time estimate");
hold off
end
