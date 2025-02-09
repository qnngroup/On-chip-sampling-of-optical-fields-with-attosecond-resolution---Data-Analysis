% Data Analysis for the 200nm antenna
% This script does the data analysis and plotting to reproduce the results
% shown in our publication "On-chip sampling of optical fields with
% attosecond resolution" - https://arxiv.org/abs/2009.06045.
% 
% 
%
% Developed by: Felix Ritzkowsky, September 2020 

clear all;
close all;
spectrumFile = "W0039.TXT"; 
%load('Data200nm2');
data2 =  load('Data200nmV3.mat');
load('retrieved_group_delay_v2.mat')
dir='figureOutput';
mkdir(dir)
c=2.998.*1e8;%speed of light in m/s





%%

dt=0.01e-15;%in s
L = 80000;
epsilon0=625000/(22468879468420441*pi); %electric constant in F/m
w=3e-6;
duration=10e-15;
wvlgth=1.179e-6;
energy=50e-12;
% duration1=10e-15;
% wvlgth1=1.179e-6;
energy1=50e-12;
cep=0;%pi/2+0.4*pi;
tau=0;
%200nm Antenna
pulse1=makePulseV3(L,1e-17,'W0039.TXT',w,energy,tau,cep+0.36.*pi); 
pulse4=makePulseV3(L,1e-17,'W0039.TXT',w,energy,tau,cep+0.36.*pi); 



% Including Marco's Phase

frqzGD=c./(wavelengths.*1e-9);
            [xData, yData] = prepareCurveData( frqzGD, group_delay );

% Set up fittype and options.
ft = fittype( 'spline' );
excludedPoints = (xData < 150000000000000) | (xData > 350000000000000);
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Exclude = excludedPoints;

% Fit model to data.
[polyFit, gof] = fit( xData, yData, ft, opts );
 figure;
plot(polyFit,frqzGD,group_delay)
xlim([150e12, 350e12])
ylim([-50, 50])
print(strcat(dir,'/','GDfit'),'-dpng','-r400')
pulse1.gdFit=polyFit;   
pulse4.gdFit=polyFit;
%pulse1.addGD();
figure(1)
plot(pulse1.t0.*1e15,fftshift(real(pulse1.E)))
hold on
pulse1.addGD();
pulse4.addGD();
plot(pulse1.t0.*1e15,fftshift(real(pulse1.E)))
hold off
xlim([-200,200])

%%
figure;
plot(pulse4.frqz,angle(pulse4.spectrum));
xlim([150e12, 350e12])
xlabel('Frequency (Hz)')
ylabel('Phase (rad)')
print(strcat(dir,'/','spectralPhasePlot'),'-dpng','-r1200')
%%
%240nm Antenna
pulse2=makePulseV3(L,1e-17,"SC111919.TXT",w,energy,tau,cep+0.36.*pi); 



%% 
pulse3=makePulseV3(L,1e-17,spectrumFile,w,energy,tau,cep); 



%% Generating 200nm Antenna Pulse
 [E200nm Atenhance f200nm Spectrum200nm]=pulsefilter(pulse1.t0,pulse1.E,'Field_Enhancement_200nm_roc_5nm_deformed.csv');
[E240nm Atenhance f240nm Spectrum240nm]=pulsefilterOrig(pulse1.t0,pulse2.E,'ConnectedTriangle_PHZ051G02_X148G55H247p6_normalE.csv');

  
%% Zeroing Phase at 250THz
FS=1e17;
padData=[zeros(1e7,1);data2.dataAvg.*tukeywin(length(data2.dataAvg),0.2);zeros(1e7,1)]; % Put in windowing tukey window 

L=length(padData);
frqz=(FS*(0:(L/2))/L);
spectrumNew=fft(fftshift(padData));
[~,indexFrqz0]=min(abs(frqz-250e12));
[~,indexFrqz0_2]=min(abs(pulse1.frqz-250e12));
frqzSim=pulse1.frqz;
phaseCC=unwrap(angle(spectrumNew(1:length(spectrumNew)/2+1)));
phaseCC=phaseCC-phaseCC(indexFrqz0);
%% Linear Fit of Phase 
fitPhaseCC=linearPhaseFit(frqz,phaseCC,200e12,300e12);
phaseCC=phaseCC-fitPhaseCC(frqz);
%% Sim

phase200nm=unwrap(angle(pulse1.spectrum));
phase200nm=phase200nm-phase200nm(indexFrqz0_2);
phase240nm=unwrap(angle(pulse2.spectrum));
phase240nm=phase240nm-phase240nm(indexFrqz0_2);

fitPhase200nm=linearPhaseFit(frqzSim,phase200nm,200e12,300e12);
phase200nm=phase200nm-fitPhase200nm(frqzSim)';

fitPhase240nm=linearPhaseFit(frqzSim,phase240nm,200e12,300e12);
phase240nm=phase240nm-fitPhase240nm(frqzSim)';

%% Time Domain Paper
% Find max
Esim=fftshift(E200nm);
tsim=pulse1.t0;
Eosa=fftshift(pulse3.E);
E2dsi=-fftshift(pulse4.E);
tosa = pulse3.t0;
[~,index1] = max(abs(Esim).^2);
[~,index3] = max(abs(Eosa).^2);
[~,indexCC] = max(abs(data2.dataAvg).^2);%max(real(backTrafo)./max((real(backTrafo))));
[~,index2dsi] = max(abs(E2dsi).^2);

figure
set(gcf,'PaperUnits','centimeters','PaperSize',[8,5],'PaperPosition',[0 0 8 5])

plot(data2.timeDelayResampled.*1e15-data2.timeDelayResampled(indexCC).*1e15,-real(data2.dataAvg)./max((real(data2.dataAvg))),tsim.*1e15-tsim(index1).*1e15,(real(Esim))./max((real(Esim))),tosa.*1e15-tosa(index2dsi).*1e15,(real((E2dsi)))./max((real(E2dsi))))
ylabel('Electric Field (V/m)')
xlabel('Time (fs)')
ylim([-1.1 1.1])
xlim([-60 100])
legend('Measurement','Simulation CC (200nm)','OSA')
 ax=gca;
ax.XMinorTick='on';
ax.YMinorTick='on';
set(gca,'FontSize',7);
box on
grid off
 print(strcat(dir,'/','PaperFieldPlot'),'-dpng','-r400')
  print(strcat(dir,'/','PaperFieldPlot'),'-depsc')
%%
figure
set(gcf,'PaperUnits','centimeters','PaperSize',[3.5,1.5],'PaperPosition',[0 0 3.5 1.5])

plot(data2.timeDelayResampled.*1e15-data2.timeDelayResampled(indexCC).*1e15,-real(data2.dataAvg)./max((real(data2.dataAvg))),tsim.*1e15-tsim(index1).*1e15,(real(Esim))./max((real(Esim))),tosa.*1e15-tosa(index2dsi).*1e15,(real((E2dsi)))./max((real(E2dsi))))
ylabel('Electric Field (V/m)')
xlabel('Time (fs)')
ylim([-1.1 1.1])
xlim([0 35])
legend('Measurement','Simulation CC (200nm)','OSA')
 ax=gca;
ax.XMinorTick='on';
ax.YMinorTick='on';
set(gca,'FontSize',6);
box on
grid off
 print(strcat(dir,'/','PaperFieldPlotInset'),'-dpng','-r400')
  print(strcat(dir,'/','PaperFieldPlotInset'),'-depsc')
%% Frequency Downsampling of Simulation
fPl=f200nm;
[~,indexF1] = min(abs(tsim+50e-15));
[~,indexF2] = min(abs(tsim-50e-15));
Lwindow=indexF2-indexF1+1;
EsimFilt= Esim;
EsimFilt(indexF1:indexF2)=tukeywin(Lwindow,0.2)'.*EsimFilt(indexF1:indexF2);
EsimFilt(1:indexF1-1)=0;
EsimFilt(indexF2+1:end)=0;
spectrumFilt=fft((EsimFilt));

plot(fPl,abs(fftshift(spectrumFilt)))
yyaxis right
plot(fPl,unwrap(angle(fftshift(spectrumFilt))))
%% Frequency Downsampling of 2DSI
f2dsi=pulse4.frqz;
[~,index3] = min(abs(pulse4.t0+50e-15));
[~,index4] = min(abs(pulse4.t0-50e-15));
Lwindow2=index4-index3+1;
E2dsiFilt= fftshift(pulse4.E);
E2dsiFilt(index3:index4)=tukeywin(Lwindow2,0.2)'.*E2dsiFilt(index3:index4);
E2dsiFilt(1:index3-1)=0;
E2dsiFilt(index4+1:end)=0;
spectrum2dsiFilt=fft(fftshift(E2dsiFilt));

figure;
yyaxis left
plot(f2dsi,abs(fftshift(spectrum2dsiFilt)))
yyaxis right
plot(f2dsi,unwrap(angle(fftshift(spectrum2dsiFilt))))

%%
spectrumPlasmon= fftshift(spectrumFilt);%Spectrum200nm;
phasePlasmon = unwrap(angle((spectrumPlasmon)));
fitPhase200nm=linearPhaseFit(fPl,phasePlasmon,200e12,300e12);
phasePlasmon = phasePlasmon - fitPhase200nm(fPl)';
phase2dsi = angle(spectrum2dsiFilt);
fit2dsi = linearPhaseFit(pulse4.frqz,phase2dsi,200e12,300e12);
phase2dsi = phase2dsi - fit2dsi(pulse4.frqz)';
[~,indexPhaseCC]=min(abs(frqz-248e12));
[~,indexPhasePlasmon]=min(abs(fPl-248e12));
[~,indexPhase2dsi] = min(abs(pulse4.frqz-248e12));

figure
set(gcf,'PaperUnits','centimeters','PaperSize',[8,5],'PaperPosition',[0 0 8 5])
yyaxis left
plot(frqz,abs(spectrumNew(1:length(spectrumNew)/2+1))./max(abs(spectrumNew(1:length(spectrumNew)/2+1))),fPl,abs(spectrumPlasmon)./max(abs(spectrumPlasmon)),pulse3.frqzOSA,smooth(normalize(pulse3.spectrumOSA,'range')))%pulse1.frqz,abs(spectrumBackup)./max(abs(spectrumBackup)),
ylabel('Spectral Amplitude (arb. Units')
ylim([0 1.1])
yyaxis right
plot(frqz,phaseCC-phaseCC(indexPhaseCC),fPl,phasePlasmon-phasePlasmon(indexPhasePlasmon),pulse4.frqz,phase2dsi-phase2dsi(indexPhase2dsi))
ylabel('Phase (rad)')
ylim([-3 3])
%legend('Plasmonic Enhancement')
xlim([150e12 350e12])
%lgd=legend('Measured CC','Pl. Enhanced Spectrum 200nm', 'OSA','Measured CC','Phase 200nm')
%lgd.Location='northeast'
xlabel('Frequency (Hz)')
 ax=gca;
ax.XMinorTick='on';
ax.YMinorTick='on';
set(gca,'FontSize',6);
box on
grid off
 print(strcat(dir,'/','SpectrumPhaseCC'),'-dpng','-r400')
  print(strcat(dir,'/','SpectrumPhaseCC'),'-depsc')

  
%% Pulse Duration
% Measured Data
spectrumNewCMPLX = spectrumNew;
spectrumNewCMPLX(end/2:end) = 0;
Ecmplx = fftshift(fft(spectrumNewCMPLX));
Ecmplx= Ecmplx(1e7+1:end-1e7);
[~,locCC] = max(Ecmplx);
% Simulated Data
spectrumPlasmonCMPLX = spectrumPlasmon;
spectrumPlasmonCMPLX(end/2:end) = 0;
EplasmonCMPLX = (fft(spectrumPlasmonCMPLX));
[~,locPL] = max(abs(EplasmonCMPLX).^2);
%OSA based Data
spectrumOSACMPLX = pulse3.spectrum;
spectrumOSACMPLX(end/2:end) = 0;
EosaCMPLX = fftshift(fft(spectrumOSACMPLX));
[~,locOSA] = max(abs(Eosa).^2);
figure
set(gcf,'PaperUnits','centimeters','PaperSize',[8,5],'PaperPosition',[0 0 8 5])

plot(data2.timeDelayResampled.*1e15-data2.timeDelayResampled(locCC).*1e15,abs(Ecmplx).^2./max(abs(Ecmplx).^2),tsim.*1e15-tsim(locPL).*1e15,abs(EplasmonCMPLX).^2./max(abs(EplasmonCMPLX).^2),tosa.*1e15-tosa(locOSA).*1e15,abs(EosaCMPLX).^2./max(abs(EosaCMPLX).^2))

ylabel('Intensity (arb. u.)')
xlabel('Time (fs)')
%ylim([0 1.1])
xlim([-20 20])
%legend('Measurement','Simulation CC (200nm)','OSA')
 ax=gca;
ax.XMinorTick='on';
ax.YMinorTick='on';
set(gca,'FontSize',6);
box on
grid off
print(strcat(dir,'/','PaperFieldPlotFWHM'),'-dpng','-r400')
print(strcat(dir,'/','PaperFieldPlotFWHM'),'-depsc')


  
   %%
 [~,indexConf] = max(data2.dataAvg.^2);
  
fig6=figure;
set(gcf,'PaperUnits','centimeters','PaperSize',[15,8],'PaperPosition',[0 0 15 8]);
confplot((data2.timeDelayResampled-data2.timeDelayResampled(indexConf)).*(1e15),-data2.dataAvg./max(-data2.dataAvg),data2.datastd./max(-data2.dataAvg));
hold on
plot(tsim.*1e15-tsim(index1).*1e15,(real(Esim))./max((real(Esim))),tosa.*1e15-tosa(index2dsi).*1e15,(real((E2dsi)))./max((real(E2dsi))))
hold off
ax=gca;
ax.XLabel.String='Delay (fs)';
%ax.XLimMode='manual';
ax.XTick=[-50:10:50];
%ax.YLimMode='manual';
axis tight;
ax.XLim=[-50 50];
ax.YLim=[-1.2 1.2 ];
ax.YLabel.String=' Electric Field (arb. u.)';


set(ax,'FontSize',11);
box on
ax.YMinorTick='on';
ax.XMinorTick='on';
print(fig6,strcat(dir,'/','ConfPlot'),'-dpng','-r1200')
print(fig6,strcat(dir,'/','ConfPlot'),'-depsc')

  


function fitresult = linearPhaseFit(frqz,phase,lowerBound,upperBound)
%% Linear Fit of Phase 
[xData, yData] = prepareCurveData( frqz, phase );

% Set up fittype and options.
ft = fittype( 'poly1' );
excludedPoints = (xData < lowerBound) | (xData > upperBound);
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft, opts );

end

