function [Output_spect] = PspectWtapers(varargin)
%{
What mtspectrumc() from Chronux is doing
changes the data to be accross column with each row a dataset
calculates the nfft
creates a grid of frequencies
dpass is run to generate the Slepian sequences
binning the data using tapers
run fft
filters the target frequencies
squeezes the data
determines the error if requested
%}
switch size(varargin,2)
    case  1
        disp('Not enough imput arguments')
        Output_spect = [];
        return
    case  2
        data = varargin{1,1};
        fpass = varargin{1,2}.freqRange;
        tapers = varargin{1,2}.tapers;
        Fs = varargin{1,2}.Fs;
        overlapSize = varargin{1,2}.overlapSize;
        binSize = varargin{1,2}.binSize;
    otherwise
        disp('Too many input arguments')
        Output_spect = [];
        return
end

STL = 1;    %starting location for the data
% pad = 1;
Stepping = floor(Fs*binSize);

stepSize = floor(overlapSize*Stepping);
AtMost = floor((size(data,1)-Stepping)/stepSize);
StepsbyStep = floor(linspace(stepSize,(size(data(:,1),1)-Stepping),AtMost));


%creating tapers
[tapers,~]=dpss(Stepping,tapers(1),tapers(2));
tapers = tapers*sqrt(Fs);
%
% %running the fft
[~,C]=size(data); % size of data
K=size(tapers,2); % size of tapers
tapers=tapers(:,:,ones(1,C)); % add channel indices to tapers


% tic
for hj = 1:AtMost-1
    amtUp = StepsbyStep(hj);

    if hj > 1
        dataRun = data(STL(1)+amtUp:STL(1)+amtUp+Stepping-1,1);
    else
        dataRun = data(STL(1):STL(1)+Stepping-1,1);
    end
        
    Nsp=sum(dataRun,1); % number of spikes in each channel
    Msp=Nsp'./Stepping; % mean rate for each channel
    dataRun=dataRun(:,:,ones(1,K));% add taper indices to the data
    dataRun=permute(dataRun,[1 3 2]); % permute data to be of the same dimensions as H
    dataRun_proj=dataRun.*tapers; % multiply data by the tapers
    %determining power using pspectrum with conditions
    leak = 1;
    %finding the upper and lower frequency limits
    
    
    [J_taper, ~] = pspectrum(tapers, ...
                Fs, ...
                'power', ...
                'FrequencyLimits',[fpass(1) fpass(2)], ...
                'Leakage',leak);
                
    meansp_pspec =Msp(:,ones(1,K),ones(1,size(J_taper,1)));  % add taper and frequency indices to meansp
    meansp_pspec =permute(meansp_pspec,[3,2,1]); % permute to get meansp with the same dimensions as H
    
    
    [J_ps, freq_ps] = pspectrum(dataRun_proj, ...
                Fs, ...
                'power', ...
                'FrequencyLimits',[fpass(1) fpass(2)], ...
                'Leakage',leak);

    
    dc_remove = J_ps-J_taper.*meansp_pspec;
    secondMethod = dc_remove(:,:,:);
    smS = squeeze(mean(conj(secondMethod),2));
    %determining power using pspectrum without conditions
    if hj == 1
        pspecWdindows = nan([size(J_ps,1),AtMost-1]);
    end
    pspecWdindows(:,hj) = smS;
    if hj == 2
        windowSize = size(dataRun,1);
        priorHz = windowSize/Fs;
        newHz_pspect = size(pspecWdindows,1)/priorHz;
    end

end

Output_spect = struct;

Output_spect.WindowSettings.windowsize = Stepping;
Output_spect.WindowSettings.bincount = AtMost-1;
Output_spect.Pspect.spect_data = pspecWdindows;
Output_spect.Pspect.freq = freq_ps;
Output_spect.Pspect.hz = newHz_pspect;


% toc


