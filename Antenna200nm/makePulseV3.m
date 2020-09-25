classdef makePulseV3 < handle
    properties
        L % vector size
        t0% Time Vector
        dt
        w % Beam width
        specFilename % Spectrum File name ".txt"
        spectrumFit % Spline Fit of Spectrum
        wvlgth % Imported wvlgth Axis
        counts % imported counts
        spectrumOSA %weighted Spectrum measured with the OSA
        frqzOSA % frqz Axis for spectrumOSA
        % Phase Import
        gdFit % Fit function for the phase
        
        energy
        tau
        cep %in radian
        E
        pow2
        PeakField
        spectrum %Complex Spectrum -> is never modified
        spectrumOriginal % Original Spectrum
        %         spectralAmp % spectral Amplitude
        %         phase % spectral phase
        frqz % Frequency Axis
        
        sod
        tod
        fod
        hbar=1.055e-34%joule*seconds - Js
        
    end
    properties (SetAccess=immutable)
        epsilon0=625000/(22468879468420441*pi); %electric constant in F/m
        c=2.998.*1e8;%speed of light in m/s
        
    end
    methods
        function obj=makePulseV3(L,dt,spectrumFile,w,energy,tau,cep)
            obj.L = L;
            obj.dt = dt;
            obj.t0= linspace(-(obj.L/2-1).*dt,(obj.L/2).*dt,obj.L);
            obj.frqz = [(0:L/2-1),-(L/2):-1]./(dt.*L);         
            obj.E=zeros(obj.pow2,1);
            obj.w=w;
            obj.energy=energy;
            obj.tau=tau;
            obj.cep=cep;
            obj.specFilename = spectrumFile;
            %% Import Spectrum
            [obj.wvlgth,obj.counts]=obj.loadSpec();
            obj.rebin();
            [obj.spectrumFit,~]=obj.createSpectrumFit();
            %% Define Pulse
            
            obj.spectrumOriginal=(obj.spectrumFit(abs(obj.frqz)))';
            index=find(obj.spectrumOriginal<0);
            obj.spectrumOriginal(index)=0;
            obj.spectrumOriginal=obj.spectrumOriginal./max(obj.spectrumOriginal)
            obj.spectrum=(obj.spectrumOriginal).*exp(1i.*cep.*sign(obj.frqz));
            obj.E=ifft(obj.spectrum);
            obj.E=(obj.E);
            obj.calcE();
            obj.fft();
            obj.inverseFT();
            
            
        end
        function n=sellmeierEqt(obj,material)
            
            % Wavelength in um!
            switch material
                case 'BBO'
                    n=@(x) sqrt(2.7405+0.0184./(x.^2-0.0179)-0.0155.*x.^2);
                case 'YAG'
                    n=@(x) sqrt(1+2.28200./(1-0.01185./x.^2)+3.27644./(1-282.734./x.^2));
                case 'SF10'
                    n=@(x) sqrt(1+1.62153902./(1-0.0122241457./x.^2)+0.256287842./(1-0.0595736775./x.^2)+1.64447552./(1-147.468793./x.^2))
                case 'Te02'
                    n =@(x) sqrt(1 + 3.71789.*x.^2./(x.^2 - 0.19619.^2) + 0.07544.*x.^2/(x.^2-4.61196.^2)); %TeO2 ordinary
                case 'SiO2'
                    n=@(x) sqrt(1+0.6961663./(1-(0.0684043./x).^2)+0.4079426./(1-(0.1162414./x).^2)+0.8974794./(1-(9.896161./x).^2));
                    
            end
        end
        function addMaterial(obj,inString,thickness)
            n=obj.sellmeierEqt(inString);
            x=(obj.c./obj.frqz).*1e6;
            phaseMat=(real(n(x)).*(2*pi)./(obj.c./obj.frqz)).*thickness;
            obj.spectrum=obj.spectrum.*exp(1i.*(phaseMat))';
            obj.inverseFT()
        end
        
        function fft(obj)
            obj.spectrum=fft((obj.E));
            obj.spectrum=obj.spectrum;
        end
        
        function inverseFT(obj)
            obj.E=ifft(obj.spectrum);
            obj.E=(obj.E);
        end
        %         function inverseFT2(obj)
        %             obj.E2=(ifft([obj.spectrum2;flip(conj(obj.spectrum2))]));
        %             obj.E2=fftshift(obj.E2);
        %         end
        function calcE(obj)
            A=trapz(obj.t0,abs(obj.E(1:end)).^2);
            obj.PeakField=sqrt(obj.energy.*2./(2*pi*((obj.w^2)./8).*A.*obj.epsilon0.*obj.c));%sqrt((2./(c.*epsilon0)).*(log(2)./pi).^(3/2) .*energy./(width.^2 *tau)); %V/m Formular from Alexs Masterthesis page 118
            obj.E=obj.E.*obj.PeakField;
            
        end
        function shiftZero(obj)
            [~,index]=max(abs(obj.E).^2);
            obj.shiftPos(-obj.t0(index));
        end
        function shiftPos(obj,t)
            %[~,index]=max(abs(obj.E2));
            %obj.spectrum=obj.spectrum.*exp(1i.*(2*pi/(1/t)).*obj.frqz');
            obj.inverseFT();
        end
        
        function out=freqzp(obj)
            out=(0:obj.L/2-1)./(obj.L.*obj.dt);
        end
        
        function out=freqzn(obj)
            out=(-obj.L/2:-1)./(obj.L.*obj.dt);
        end
        
        
        function [wvlgths, counts] = loadSpec(obj)
            % Script for importing data from the following text file:
            %
            %    filename: /Users/felix/Documents/MATLAB/AutoCorrDataAF/DataAnalysis/SC111919.TXT
            %
            % Auto-generated by MATLAB on 09-Apr-2020 09:42:52
            
            %% Setup the Import Options
            opts = delimitedTextImportOptions("NumVariables", 2);
            
            % Specify range and delimiter
            opts.DataLines = [4, 1004];
            opts.Delimiter = ",";
            
            % Specify column names and types
            opts.VariableNames = ["wvlgths", "counts"];
            opts.VariableTypes = ["double", "double"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";
            
            % Import the data
            tbl = readtable(obj.specFilename, opts);
            
            %% Convert to output type
            wvlgths = tbl.wvlgths;
            counts = tbl.counts;
            
            %% Clear temporary variables
            clear opts tbl
        end
        function [fitresult, gof] = createSpectrumFit(obj)
            %CREATEFIT(FRQZNEAR,NEARFSPECTRUM)
            %  Create a fit.
            %
            %  Data for 'untitled fit 1' fit:
            %      X Input : frqzNear
            %      Y Output: nearFSpectrum
            %  Output:
            %      fitresult : a fit object representing the fit.
            %      gof : structure with goodness-of fit info.
            %
            %  See also FIT, CFIT, SFIT.
            
            %  Auto-generated by MATLAB on 20-May-2020 10:00:22
            
            
            %% Fit: 'untitled fit 1'.
            [xData, yData] = prepareCurveData(obj.frqzOSA, obj.spectrumOSA );
            
            % Set up fittype and options.
            ft = fittype( 'smoothingspline' );
            opts = fitoptions( 'Method', 'SmoothingSpline' );
            opts.SmoothingParam = 1.48180537806491e-38;
            
            % Fit model to data.
            [fitresult, gof] = fit( xData, yData, ft, opts );
            
            % Plot fit with data.
%             figure( 'Name', 'untitled fit 1' );
%             h = plot( fitresult, xData, yData );
%             legend( h, 'nearFSpectrum vs. frqzNear', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
%             % Label axes
%             xlabel( 'frqzNear', 'Interpreter', 'none' );
%             ylabel( 'nearFSpectrum', 'Interpreter', 'none' );
%             grid on
            
        end
        function rebin(obj)
            binSize=[];
            f=(obj.c./(obj.wvlgth.*1e-9));
            
            
            for i=2:(length(obj.wvlgth)-1)
                binSize(i)=-((f(i)-f(i-1))./2)+(f(i+1)-f(i))/2;
            end
            
            
            % Use this if you do not want to weigh the spectrum in the
            % focus
            
            obj.spectrumOSA=sqrt((((10.^(obj.counts(2:end)./10))./binSize')./5e-13));
            obj.spectrumOSA(isinf(obj.spectrumOSA))=NaN;
            % Use this for weighing in the focus - does not work at the moment ???
            
            %             nearFSpectrum=sqrt((((10.^(obj.counts(2:end)./10))./binSize')./5e-13).*(pi*f(2:end)/(8*0.25.*obj.c.^2))./(obj.w.*1e3));
            %             obj.spectrumOSA(isinf(nearFSpectrum))=NaN;
            obj.frqzOSA=f(2:end);
            
        end
        function loadPhase(obj,filename)
            %% Setup the Import Options
            opts = delimitedTextImportOptions("NumVariables", 2);
            
            % Specify range and delimiter
            opts.DataLines = [1, Inf];
            opts.Delimiter = " ";
            
            % Specify column names and types
            opts.VariableNames = ["VarName1", "VarName2"];
            opts.VariableTypes = ["double", "double"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";
            opts.ConsecutiveDelimitersRule = "join";
            opts.LeadingDelimitersRule = "ignore";
            
            % Import the data
            tbl = readtable(filename, opts);
            
            %% Convert to output type
            wvlgth = tbl.VarName1;
            gd= tbl.VarName2;
            
            %% Clear temporary variables
            clear opts tbl
            frqz1=obj.c./(wvlgth.*1e-9);
            [xData, yData] = prepareCurveData( frqz1, gd );
            
            % Set up fittype and options.
            ft = fittype( 'smoothingspline' );
            opts = fitoptions( 'Method', 'SmoothingSpline' );
            opts.SmoothingParam = 2.3908479310290142E-39;%4.80214643964101e-38;
            
            % Fit model to data.
            [splineFit, gof] = fit( xData, yData, ft, opts );
            
            % Plot fit with data.
%             figure( 'Name', 'untitled fit 1' );
%             h = plot( splineFit, xData, yData );
%             legend( h, 'VarName2 vs. frqz', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
%             % Label axes
%             xlabel( 'frqz', 'Interpreter', 'none' );
%             ylabel( 'VarName2', 'Interpreter', 'none' );
            grid on
            
%             %Polynomial Fit
%             %% Fit: 'untitled fit 1'.
%             [xData, yData] = prepareCurveData( frqz1, gd );
%             
%             % Set up fittype and options.
%             ft = fittype( 'poly6' );
%             index = find( frqz1 < 210e12);
%             index = [index ; find(frqz1 > 320e12)];
%             excludedPoints = excludedata( xData, yData, 'Indices', index );
%             opts.Exclude = excludedPoints;
%             % Fit model to data.
%             [polyFit, gof] = fit( xData, yData, ft );
%             
%             % Plot fit with data.
%             figure( 'Name', 'untitled fit 1' );
%             h = plot( polyFit, xData, yData, excludedPoints );
%             legend( h, 'VarName2 vs. frqz', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
%             % Label axes
%             xlabel( 'frqz', 'Interpreter', 'none' );
%             ylabel( 'VarName2', 'Interpreter', 'none' );
%             grid on
            obj.gdFit = splineFit;
            
            
            
            
        end
        
        function addGD(obj)
            L = length(obj.spectrum);
            [~,index] = min(abs(obj.frqz(1:L/2+1)-(obj.c./(1179e-9))));
            gdF = obj.gdFit(abs(obj.frqz(1:L/2+1))).*1e-15;
            phase = cumtrapz((1/(L.*obj.dt)),gdF).*2.*pi;% units of s
            phase = phase - phase(index);
            disp(size(phase))
            obj.spectrum = [obj.spectrum(1:L/2).*exp(1i.*phase(1:L/2))', obj.spectrum(L/2+1:end).*flip(conj(exp(1i.*phase(2:L/2+1))))'];

            obj.inverseFT();
        end
        
    end
end
