function PSDanalysis(TR,signal,N)
% N = number of time points

Fs = 1/TR; % Sampling frequency (1/TR)
% Lower bound is DC (0 Hz) and the Upper bound is Nyquist Frequency,
% that is, heoretically minimum number of points that one needs to sample in
% order to see the real frequency of the sample. --> 
hz = linspace(0,Fs/2,N/2+1); %The formula to convert indices to Hz 

fhat = fft(signal);
P2 = abs(fhat);
P1 = P2(1:floor(N/2)+1,:);
P1(2:end-1,:) = P1(2:end-1,:).^2;   

figure;
plot(hz,P1);


% figure;
% c = copper(rois);
% 
% rois = size(signal,2);
% for seed=1:rois
%         % This gaussian filter ithink is to smooth out the actual fequencies, will have to check
%         Power_Areas(:,seed)=gaussfilt(hz,P1(:,seed)',0.01);
%         plot(hz,Power_Areas(:,seed),'color',c(seed,:));
%         hold on
% end
% grid on
% mP = mean(Power_Areas,2);
% 
% plot(hz,mP,'color','black','LineWidth',2.5);
% hold off
% set(gca,'XTickLabel',[],'YTickLabel',[]);
