%{
This script is used to show an example for using PspectWtapers
Inputs:
filtered data
setting input
    frequency rang = [3 150]
    tapers = [5 9];
    Fs = 1375;
    overlapSize = 0.01;
    binSize = 0.5;
%}

%asking the user for the loaction of the data
data = DataSet.PD02D1E1.Fullfiltered(1:8000);

setting.freqRange = [13 30];
setting.tapers = [5 9];
setting.Fs = 1375;
setting.overlapSize = 0.01;
setting.binSize = 0.5;

OutputSpectrum = PspectWtapers(data,setting);