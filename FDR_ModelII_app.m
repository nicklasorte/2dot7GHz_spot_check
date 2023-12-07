function [FDR_dB,ED,VD,OTR,DeltaFreq,single_fdr_loss,trans_mask] = FDR_ModelII_app(app,freq_sep,array_tx_rf,array_rx_if,array_tx_mask,array_rx_insert_loss,TXSlope,RXSlope)
LogLin = 1; 

CF = 0;
Delta_Freq_Step=1; %%%%%%%%1 MHz Step Size
temp_freq_sep=max([max(array_rx_if)+1,max(array_tx_rf)+1,abs(freq_sep)]);
end_freq=abs(temp_freq_sep);
num_freq_steps=(1/Delta_Freq_Step)*end_freq*2+1;
start_freq = -1*end_freq; 

DeltaFreq=start_freq:Delta_Freq_Step:end_freq;
band_center_idx=ceil(length(DeltaFreq)/2);
NumStep=length(DeltaFreq);

DF{:}=DeltaFreq;

%Set Data Length
LF = length(DF{:});
FDR_dB = zeros(LF,1);

%Collect Curve Data
BW_TX = [-array_tx_rf flip(array_tx_rf(1:(find(array_tx_rf==0,1)-1)))]';
BW_RX = [-array_rx_if flip(array_rx_if(1:(find(array_rx_if==0,1)-1)))]';
Mask_TX = [array_tx_mask flip(array_tx_mask(1:(find(array_tx_mask==0,1)-1)))]';
Reject_RX = -[array_rx_insert_loss flip(array_rx_insert_loss(1:(find(array_rx_insert_loss==0,1)-1)))]';

%Calculate BandWidths
idx1 = find(BW_RX == 0,1);
RXidx = idx1-1;
RX_BW = array_rx_if(RXidx)*2;
idx2 = find(BW_TX == 0,1);
TXidx = idx2-1;
TX_BW = array_tx_rf(TXidx)*2;



%Find Smallest Frequency Offset
Fmin =  min(min(array_tx_rf(1:(TXidx))),min(array_rx_if(1:(RXidx))));

%Calculate Factor
dF = 10^(floor(log10(Fmin)))/10;

if dF > 1
    dF=1;
end

%Calculate some values to deal with NaNs
[BRx,BwRx]=Extrap0dB(app,flip(-array_rx_insert_loss(1:RXidx)),flip(array_rx_if(1:RXidx)),RXSlope);
[BTx,BwTx]=Extrap0dB(app,flip(-array_tx_mask(1:TXidx)),flip(array_tx_rf(1:TXidx)),TXSlope);

for zz = 1:LF
    zz/LF*100
    %Collect Tune Frequency to evaluate
    RXCF = DF{1}(zz); %Can use a different column number to normalize?
    TXCF = CF;
    
    [Txfl,Txfh]= flfh(app,TXCF,BTx,BwTx,TXSlope);
    [Rxfl,Rxfh]= flfh(app,RXCF,BRx,BwRx,RXSlope);
    
    Fl=min(Txfl,Rxfl); Fh=max(Txfh,Rxfh);
    
    [m,~] = size(BRx);
    sumRx = 0;
    
    %Calculate TX Curve
    F_Tx = SetRange(app,TXCF,Fl,Fh,dF);
    S_Tx = makespectrumII(app,F_Tx,TXCF,BTx,BwTx,TXSlope,LogLin);

    if TXCF<=Fl || TXCF>=Fh
        TXCF = Fo(app,F_Tx,S_Tx);
    end
    
    F_Rx = SetRange(app,RXCF,Fl,Fh,dF);
    
    if LogLin %use linear or log extrap
        tmp = interp1(F_Tx,S_Tx,F_Rx,'linear');
    else
        tmp = interp1(F_Tx,S_Tx,F_Rx,'linear','extrap');
    end
    
    S_Tx = tmp;
    
    % Calculate Rx curve
    F_Rx = SetRange(app,RXCF,Fl,Fh,dF);
    S_Rx = zeros(m,length(F_Rx));
    
    for i=1:m
        S_Rx(i,:)= makespectrumII(app,F_Rx,RXCF,BRx(i,:),BwRx(i,:),RXSlope,LogLin);
        sumRx = sumRx + S_Rx(i,:);
    end
    
    CFmin = TXCF-(RX_BW/2); CFmax = TXCF+(RX_BW/2);
    
    if abs(BRx(1)) == 3
        RX_BW = BwRx(1)*2;
    end
    
    
    if RX_BW >= TX_BW
        if RXCF >= CFmin && RXCF <= CFmax
            FDR_dB(zz) = 0;
            continue
        end
    end
    
    %Calculate Tx signal that gets through Rx
    output = S_Tx + sumRx;
    
    %Calculate FDR
    FDR_dB(zz) = 10*log10(sum(10.^(S_Tx/10))/sum(10.^(output/10)));  %Calculate curve by stepping RXCF thru DF (Add OTR Factor)
    
    
end

OTR = min(FDR_dB);
BW_TX
Mask_TX
ED = [ BW_TX, Mask_TX];
VD = [ BW_RX, Reject_RX];

ind_fdr_loss=nearestpoint_app(app,freq_sep,DeltaFreq);
single_fdr_loss=FDR_dB(ind_fdr_loss);

trans_mask=horzcat(F_Tx',S_Tx');

end

function [b,bw]=Extrap0dB(app,b,bw,Slope)
% b is vector of X dB levels, bw is associated full bandwidths
% If b does not start with 0 this function produces 0 dB bw by extrapolation
% any rows or columns that are all nan should be deleted
% for rows where 1st col is nan, extrapolate if at least 2 non-nan

[m,n]=size(b);

for i=1:m
    nans=isnan(bw(i,2:n))|isnan(b(i,2:n));
    keep=logical([1 ~nans]);
    b(i,1:n)=[b(i,keep) nan(1,sum(nans))];
    bw(i,1:n)=[bw(i,keep) nan(1,sum(nans))];
end

for i=m:-1:1
    if isnan(bw(i,1))
        if isempty(Slope)
            if sum(~isnan(bw(i,2:n)))<2 % delete row if < 2 non-NaN
                bw=bw(1:i-1,:);
                b=b(1:i-1,:);
            end
        else
            if sum(~isnan(bw(i,2:n)))<1 % delete row if < 1 non-NaN
                bw=bw(1:i-1,:);
                b=b(1:i-1,:);
            end
        end
    else
        if isempty(Slope)
            if sum(~isnan(bw(i,2:n)))<1
                bw=bw(1:i-1,:);
                b=b(1:i-1,:);
            end
        end
    end
end

[m,~]=size(b);

for i=1:m
    if isnan(bw(i,1))
        if isnumeric(Slope) && ~isempty(Slope)
            b1 = b(i,1:find(isnan(b(i,2:end)), 1 ));
            bw1= bw(i,1:find(isnan(bw(i,2:end)), 1 ));
            b1 = [b1 b1(end)+Slope]; %#ok
            bw1= [bw1 10*bw1(end)];  %#ok
            if sum(~isnan(bw1(2:3)))==2
                b(i,1)=0;
                bw(i,1)=10^interp1(b1(2:3),log10(bw1(2:3)),b(i,1),'linear','extrap');
            end
        else
            if sum(~isnan(bw(i,2:3)))==2
                b(i,1)=0;
                bw(i,1)=10^interp1(b(i,2:3),log10(bw(i,2:3)),b(i,1),'linear','extrap');
            end
        end
    end
end
end

function [fl,fh] = flfh(app,Fc,B,Bw,Slope)
% this returns high and low frequency to ensure at least minB dB bw
% assumes monotonically decreasing spectrum envelope

minB= -80;

[r,~]=ind2sub(size(Bw),find(Bw==max(max(Bw))));
B=B(r,:);
Bw=Bw(r,:);
keep=~isnan(B)&~isnan(Bw);
B=B(keep);
Bw=Bw(keep);
if min(B)>minB
    B=[B minB];
    Bw=[Bw 0];
    if isempty(Slope)
        Bw(end)=10^interp1(B(end-2:end-1),log10(Bw(end-2:end-1)),B(end),'linear','extrap');
    else
        Bw(end)=10^((B(end)-B(end-1))/Slope+log10(Bw(end-1)));
    end
end
tmp1=Fc-max(max(Bw/2));%fl
tmp2=Fc+max(max(Bw/2));%fh
for n=0:5
    if round(tmp1*10^n)/10^n<round(tmp2*10^n)/10^n
        break
    end
end
fl=round(tmp1*10^n)/10^n;
fh=round(tmp2*10^n)/10^n;
end

function f = SetRange(app,fc,fl,fh,df)
ff = fc:df:fh;
fff = fc:-df:fl;
f = [fff(length(fff):-1:2) ff];
end

function S = makespectrumII(app,f,fc,b,bw,dB_dec,LogLin)
% fc is center freq, fr is freq range, fs is freq step
% b are X dB levels, bw are associated full bandwidths
% b assumed to be monotonically decreasing
% S is power with 0 dB as peak, f is freq vector
% dB_dec is roll-off after last b
% dB_dec < 0 for decreasing slope.  dB_dec > 0 for increasing
% Use dB_dec = [] to continue same slope between b(end-1) and b(end)

if b(1)~=0 % insert Zero if needed
    b=[0 b];
    bw=[nan bw];
end

tmp=[0 isnan(b(2:end))|isnan(bw(2:end))]; %check for NaNs
if sum(tmp)
    b=b(1:find(tmp==1, 1 )-1);
    bw=bw(1:find(tmp==1, 1 )-1);
end

if isnumeric(dB_dec) && ~isempty(dB_dec) %Add Extrapolation point with a bw point at + one decade
    b=[b b(end)+dB_dec];
    bw=[bw 10*bw(end)];
end

if isnan(bw(1)) % If 1 is NaN perform interpolation between ponints 2 and 3
    bw(1)=10^interp1(b(2:3),log10(bw(2:3)),0,'linear','extrap');
elseif bw(1)==0 % If 1 is a zero make it equal to the step size
    bw(1)=(f(2)-f(1));
end

if LogLin==1
    S = zeros(size(f));
    %S(f~=fc)=interp1(log10(bw/2),b,log10(abs(f(f~=fc)-fc)),'linear','extrap');
    S(f~=fc)=interp1(log10(bw),b,log10(abs(f(f~=fc)-fc)),'linear','extrap');
else
    S=interp1((bw/2),b,(abs(f-fc)),'linear','extrap');
end

S(S>0)=0;

end

function Y = Fo(app,F,S)
% This function guesses the center frequency of measured spectrum
Y = round(F(round((find(S==max(S), 1, 'last' )+find(S==max(S), 1 ))/2)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
