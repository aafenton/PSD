%computes Power Spectral Density from files in directory
%optionally processes only the good data
%exports values as uV^2 / Hz
%doesn't use artifact rejection
clear; close all; clc;

%input directory
drMAT = 'MAT/';
%output directory
drOUT = 'PSD/';

%list files
files = dir([drMAT '*.mat']);

%number of channels to process
Nch = 16;

eegFS = 2000; %sample rate
winLen = 4*eegFS; %analysis window

Twindow = 1; %subwindow (pwelch algorithm)
Nwindow = Twindow * eegFS; %length of subwindow in samples
Noverlap = round ( Nwindow /2); % 50% overlap
NFFT = 2^11; %FFT resolution
[Pxx,f] = pwelch(rand(1,winLen), Nwindow, Noverlap, NFFT, eegFS); %find frequencies
kf = f <= 180; %limit frequencies to 180Hz
f = f(kf);

%through files
for fI = 1:length(files)
    
    %load MAT file
    load([drMAT files(fI).name]);
    
    %display file loaded
    disp(files(fI).name);
    
    resultsPSD = cell(1,Nch);
  
    %through channels
    for chI = 1:Nch
        
        %select channel data - values in volts
        eeg = eegData(chI,:);
        
        %convert to microvolts
        eeg = eeg * 1e6;
        
        PxxCH = zeros(1,length(f)); %to hold all results
            
        %how many windows fits in
        Nwin = floor(length(eeg) / winLen);
            
        %through subwindows
        for i = 1:Nwin
            st = (i-1)*winLen + 1;
            ed = i*winLen;
                
            %pwelch algorithm
            [Pxx,~] = pwelch(eeg(st:ed), Nwindow, Noverlap, NFFT, eegFS);
            Pxx = Pxx(kf); %select the frequency of interest
                
            Pxx = Pxx';                           
            PxxCH = PxxCH + Pxx;
        end
        
        %take average
        PxxCH = PxxCH / Nwin;
        
        %store results
        resultsPSD{chI} = PxxCH;
        
    end %ch
    
    save([drOUT files(fI).name], 'resultsPSD', 'f', 'NFFT','Nwindow','Noverlap');
    
end %files

