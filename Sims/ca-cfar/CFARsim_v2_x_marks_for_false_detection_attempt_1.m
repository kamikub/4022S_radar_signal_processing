clear all;
%PARAMETERS
times=256;%time axis ITERATIONs
pfa1 = 10^-3;%desired Probabability of false alarm
RefWindow1 = 32;%size od data cells
gcell1=3;% number of guardcells
len = 1000;%freq samples per time
l=linspace(0,len,len);


%generate noise
noisesig=randn([times,len])+1i*randn([times,len]);
signal=zeros(times,len);
signal=noisesig;


%RUN CFAR ON EACH ITERATION
count1 =[];%number of detections per iteration
err=[];%pfa error if no signal otherwise the signal counts as false alarms so doesn't work
detect=zeros(times,11);
threshold=zeros(times,len);
for e=1:times
    y_in=signal(e,:);
    [count1(e),err(e),detect(e,:),treshold(e,:)]=CACFAR(pfa1,RefWindow1,gcell1,y_in);
end
avg_pfa_error=sum(err)/times;%averagepfa error only useful with no signal

%%arbitrary axis for the spectogram of the simulated noise
tax=0:times;
fax=0:len;

%%plotting it
imagesc(tax,fax,20*log10(abs(signal)));%it looked prettier this way to get the signal dB instead of the square wave detector thing
colorbar;
colormap('jet');
xlabel("time (s)");
ylabel("frequency (kHz)");
hold on;
%marking the detections
total_fa=0;
for ix=1:times
    for xi=1:11
        temp=detect(ix,xi);
        if temp==0
            
        else
            plot(ix,detect(ix,xi),'kx','MarkerSize',12);
            total_fa=total_fa+1;
        end
        
    end
end

%CFAR-PROCESSING
function [count,pfaerr,intercept,thold]=CACFAR(pfa,RefWindow,gcell,y_in) 
%CFAR PROCESSING
    %INTERNAL PARAMS
    sidecell=gcell+RefWindow;
    len=length(y_in);
    y_afterpow = abs(y_in).^2;
    thold=zeros(1,len);
    count=0;
    alpha=(pfa^(-1/(RefWindow*2))-1);
    aacl=0;%because the initial data point are skipped the number of CUTs that there are will be less
    intercept=zeros(1,11);
    %PROCESS
    for i=1+sidecell:len-sidecell
        lag=y_afterpow(i-sidecell:i-1-gcell);
        lead=y_afterpow(i+1+gcell:i+sidecell);
        mlag=sum(lag);
        mlead=sum(lead);
        gi=(mlead+mlag);%max(mlead,mlag);%change between CA-Cfar and OSGO-CFAR
        thold(i)=alpha*gi;
        if (y_afterpow(:,i) >= thold(i))
            intercept(count+1)=i;
            count=count+1;
        end
        aacl=aacl+1;
    end
    apfa = count/aacl;
    expfa = pfa*len;
    pfaerr=0;
    pfaerr = abs(((pfa-apfa)/pfa)*100);
    
end
