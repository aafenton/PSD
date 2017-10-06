%computes Power Spectral Density from files in directory
%optionally processes only the good data
clear all; close all; clc;

%input directory
drMAT = 'MAT/';
%artifact directory
drART = 'ART/';
%output directory
drOUT = 'PSD/';

%list files
files = dir([drMAT '*.mat']);

%number of channels to process
Nch = 14;

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
    
    %load artifact file
    load([drART files(fI).name]);
    
    %display file loaded
    disp(files(fI).name);
    
    resultsPSD = cell(1,Nch);
  
    %through channels
    for chI = 1:Nch
        
        %select channel data
        eeg = eegData(chI,:);
        
        %select channel artifacts
        winOK = windowsOK{chI};
              
        %get length of available windows
        d = winOK(:,2) - winOK(:,1);
        
        %select only long enough good windows
        k = d >= winLen;
        winOK = winOK(k,:);
        
        PxxCH = zeros(1,length(f)); %to hold all results
        winCount = 0; %to hold number of windows averaged
        
        %through good data windows
        for wI = 1:size(winOK,1)
         
            stW = winOK(wI,1); %start and end of the window
            edW = winOK(wI,2);
            d = edW-stW+1; %length of the window
            
            %how many windows fits in
            Nwin = floor(d / winLen);
            
            %through subwindows
            for i = 1:Nwin
                st = (i-1)*winLen + 1 + stW;
                ed = i*winLen + stW;
                
                %pwelch algorithm
                [Pxx,~] = pwelch(eeg(st:ed), Nwindow, Noverlap, NFFT, eegFS);
                Pxx = 20*log10(Pxx); %convert to decibels
                Pxx = Pxx(kf); %select the frequency of interest
                
                Pxx = Pxx';
                                                           
                PxxCH = PxxCH + Pxx;
                winCount = winCount + 1;
            end
        end
        
        %take average
        PxxCH = PxxCH / winCount;
        
        %store results
        resultsPSD{chI} = PxxCH;
        
    end %ch
    
    save([drOUT files(fI).name], 'resultsPSD', 'f', 'NFFT','Nwindow','Noverlap');
    
end %files

