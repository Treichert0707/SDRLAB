%% 2 Section 1 
%cc = 1; % not needed
clc; clear;
b = 10;
x=5;
n= 2;% provide an option for choosing the no. of users, and then make it two. 
if n==0
    error('No User for service');
end
if x==0
    error('ERROR:Require atleast a PN spread of 1');
end
if b==0
    error('ERROR: Cannot have no data bits');
end
%% 
Data=randi([0 1],b,n); %creates random data
PN=randi([0 1],x,n); %each column specifies a user
%RPN=randi([0 1],x,1);
PNMatrix = PN;
%RandomPN=RPN;
%% create the spreading code for the length of the data bits
for ii=2:b
    PNMatrix = [PNMatrix ; PN]; %replicate itself "b" times
end
%%
for ii=b*x:-1:1 % further extend length by a 100 for visualization
    for jj=1:n
        PNMatrix((ii-1)*100+1:ii*100,jj)=PNMatrix(ii,jj);
    end
end
%%
for ii=b:-1:1 % extending data by the same specification
    for jj=1:n
        DataMatrix((ii-1)*x*100+1:ii*x*100,jj)=Data(ii,jj);
    end
end
%% 2 section c) XOR The Data!
signal = xor(DataMatrix,PNMatrix); % XOR the data matrix with the PN matrix
Signal1=ones(x*b*100,n); % if it was a 0, we want to send a -1
for jj=1:n
    for ii=1:x*b*100
        if signal(ii,jj)==0
           Signal1(ii,jj)=-1; % now we are sending a -1
        end
    end
end
%% 2 Section d) If we wanted to send the signals together 
TransSignal=2*sum(signal,2)'/n; % sending the bits of data of all streams at the same time 
TransSignal=TransSignal-1;
X=length(PNMatrix(:,1));
figure(1)
    subplot(4,1,1)
    h5 = plot(DataMatrix(:,1),'linewidth',3);
    legend([h5], {'Data 1'},'Location','northeast')
    title('Data 1')
    set(gca,'XTick',0:x*100:X); %creates markers every length of a data pulse
    set(gca,'XTickLabel',0:b); % number of data periods
    ylabel('Data')
    subplot(4,1,2)
    h6 = plot(PNMatrix(:,1),'linewidth',3);
    legend([h6], {'PN 1'},'Location','northeast')
    title('PN 1')
    set(gca,'XTick',0:x*100:X); %creates markers every length of a data pulse
    set(gca,'XTickLabel',0:b); % number of data periods
    ylabel('Data')
    subplot(4,1,3)
    h7 = plot(Signal1(:,1),'linewidth',3);
    ylim([-1.5 1.5]);
    legend([h7], {'Signal 1'},'Location','northeast')
    title('Signal 1')
    set(gca,'XTick',0:x*100:X); %creates markers every length of a data pulse
    set(gca,'XTickLabel',0:b); % number of data periods
    ylabel('Data')
    subplot(4,1,4)
    h8 = plot(TransSignal,'linewidth',3);
    ylim([-1.5 1.5]);
    legend([h8], {'Transmitted Signal'},'Location','northeast')
    title('Transmitted Signal')
    set(gca,'XTick',0:x*100:X); %creates markers every length of a data pulse
    set(gca,'XTickLabel',0:b); % number of data periods
    ylabel('Data')
    sgtitle('User 1 Data, PN, Signal, and Concatentated Signal')
    hold off;
figure(2)
    subplot(4,1,1);
    h13 = plot(DataMatrix(:,2),'linewidth',3);
    legend([h13], {'Data 2'},'Location','northeast')
    title('Data 2')
    set(gca,'XTick',0:x*100:X); %creates markers every length of a data pulse
    set(gca,'XTickLabel',0:b); % number of data periods
    ylabel('Data')
    subplot(4,1,2);
    h14 = plot(PNMatrix(:,2),'linewidth',3);
    legend([h14], {'PN 2'},'Location','northeast')
    title('PN 2')
    set(gca,'XTick',0:x*100:X); %creates markers every length of a data pulse
    set(gca,'XTickLabel',0:b); % number of data periods
    ylabel('Data')
    subplot(4,1,3);
    h15 = plot(Signal1(:,2),'linewidth',3);
    legend([h15], {'Signal 2'},'Location','northeast')
    title('Signal 2')
    ylim([-1.5 1.5]);
    set(gca,'XTick',0:x*100:X); %creates markers every length of a data pulse
    set(gca,'XTickLabel',0:b); % number of data periods
    ylabel('Data')
    subplot(4,1,4);
    h16 = plot(TransSignal,'linewidth',3);
    ylim([-1.5 1.5]);
    legend([h16], {'Transmitted Signal'},'Location','northeast')
    title('Transmitted Signal')
    set(gca,'XTick',0:x*100:X); %creates markers every length of a data pulse
    set(gca,'XTickLabel',0:b); % number of data periods
    ylabel('Data')
    sgtitle('User 2 Data, PN, Signal, and Concatentated Signal')
    hold off;

    
% Now we move to the receiver
%% 3 secion 1 
pn_2_use = PNMatrix(:,1); %the PN code that should be utilized depending on the scenarios...
noise = wgn(length(TransSignal),1,0.1,'linear');
rec_signal = TransSignal + noise';% the received signal, change this based on the scenarios
variance = var(noise)
X=length(pn_2_use);
figure(3)
plot(rec_signal);
set(gca,'XTick',0:x*100:X); %creates markers every length of a data pulse
set(gca,'XTickLabel',0:b); % number of data periods
axis tight;
ylabel('Data')
title('Noisy Signal')

%%
%Regenerating the Received Signal if we had any noise...
for ii=1:X
    if rec_signal(ii)>=.5
        rec_signal(ii)=1;
    elseif rec_signal(ii)<=-.5
        rec_signal(ii)=-1;
    else
        rec_signal(ii)=0;
    end
end     
%%
%Recovering the signal from the PN code
decode_signal=rec_signal;
for ii=1:X
    if pn_2_use(ii)==1
        decode_signal(ii)=-rec_signal(ii);
    end
end
%%
%Integration to evaluate the correct data  
for ii=1:b
    var=(x*100*(ii-1)+1):x*100*ii;
    integrand_output(var)=cumtrapz(decode_signal(var));    
end
%%
%The threshold may need to be changed based on the noise floor, 
%since we are choosing the different cases in which 
%this is the best case
ii=0;
threshold=10;
for jj=1:b
    if integrand_output((ii+1):(100*x*jj))<= threshold
        data_eval=0;
    elseif integrand_output((ii+1):(100*x*jj))>=(-threshold)
        data_eval=1;
    else
        data_eval=2;
    end
    switch (data_eval)
        case 0
            data((ii+1):(100*x*jj))=0;
        case 1
            data((ii+1):(100*x*jj))=1;
        otherwise
            data((ii+1):(100*x*jj))=NaN;            
    end
    ii=100*x*jj;
end 
%%
%An example plot to help you in creating plots...
figure(4)
subplot(2,1,1)
plot(DataMatrix(:,1),'linewidth',3);
set(gca,'XTick',0:x*100:X); %creates markers every length of a data pulse
set(gca,'XTickLabel',0:b); % number of data periods
axis tight;
legend('Original Data');
ylabel('Data')
subplot(2,1,2)
plot(data,'linewidth',3);
set(gca,'XTick',0:x*100:X); %creates markers every length of a data pulse
set(gca,'XTickLabel',0:b); % number of data periods
axis tight;
legend('Decoded Data');
ylabel('Data')
sgtitle('Decoding User 1s Data Using The Concatenated Signal with noise')
BER = (length(DataMatrix)-sum(DataMatrix(:,1)==data'))/(length(data));
if  BER >= 0.5
    BER = 1 - BER
else
    BER = BER
end
figure(5)
plot(integrand_output, 'linewidth',3)
set(gca,'XTick',0:x*100:X); %creates markers every length of a data pulse
set(gca,'XTickLabel',0:b); % number of data periods
axis tight;
legend('Integrand Output');
figure(6)
plot(noise)
set(gca,'XTick',0:x*100:X); %creates markers every length of a data pulse
set(gca,'XTickLabel',0:b); % number of data periods
axis tight;
legend('Noise');
title('Noise')