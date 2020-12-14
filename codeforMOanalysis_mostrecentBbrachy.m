%TDMS PeakFinder With Duration  Reads TDMS Formatted Long EOD Files and Extracts Peaks.
%   Re-written for doing morpholino data analysis from new recording
%   software (dtd_eod.exe) written in LabView by JRG.
%
%   Written by J. Gallant
%   $Revision 1.00$  $Date: 2016/06/09 14:22$

final_data=[];
prompt = '     > Please Enter a Threshold Value for this Recording: ';
prompt2= '     > Please Enter the Gain for this Recording: ';
prompt3= '     > Please Enter the Average EOD Duration For this Species (seconds) ';
prompt4= '     > Please Enter the End Time for this Recording: ';
%set the directory to open data from- this needs to be 1 level up from the
%BBRACH folder you want to analyze
directory=uigetdir('C:/Users/Savvas Constantinou/Dropbox (MSU Efish Lab)/Projects/Morpholinos/','Select directory for all individual fish');

%Look for species directories (start with B*)"  %with this code as written,
%it is looking for a directory within the directory you provide that starts
%with BBRACH - 
subjectdirectory=fullfile(directory,'BBRACH*');
listing = dir(subjectdirectory);
%Extract the Names
specieslist={listing.name};
%Count the subjects
[~,specieslistsize]=size(specieslist);
eodpreview=figure();
set(eodpreview, 'WindowStyle','docked')

%for each species directory...
for a=1:specieslistsize
    treatdirectory=fullfile(directory,specieslist(a),'*.tdms'); %create the path
    treatlisting=dir(treatdirectory{1}); % list the files with the tdms extension
   % treatlisting=treatlisting(3:end); %remove .. and . directories
    treatlist={treatlisting.name}; %get the filenames
    [~,treatlistsize]=size(treatlist); %construct a list of the files
    
    fprintf('Starting step %s: Analysis of %s!\n',num2str(a),specieslist{a}); %Starting step 1: Analysis of BBRACH!
    
    %this was off in jasons code
    %raw_data(a).name=subjectlist{a};
    
    for b=1:treatlistsize  %for each file....
        
        %filedirectory=fullfile(directory,specieslist(a),treatlist(b),'*.tdms'); %create the path
        %filelisting=dir(filedirectory{1}); % list the files with the tdms extension
        %filelist={filelisting.name}; %get the filenames
        %[~,filelistsize]=size(filelist); %construct a list of the files
        
        %fprintf('  Step %s.%s: Analysis of %s %s started...\n',num2str(a),num2str(b),specieslist{a},treatlist{b});  %Step 1.1 Analysis: Analysis of BBRACH CRISPR
        
        %for c=1:filelistsize
            fprintf('   Step %s.%s: Analysis of %s started...\n',num2str(a),num2str(b), treatlist{b});  %Step 1.1.1 %Analyis of BBRACH CRISPR file1.tdms
            eodwave=[]; %initialize variables
            wavesize=[];
            neg_peakvoltages=[];
            neg_peak_idx=[];
            pos_peakvoltages=[];
            pos_peak_idx=[];
            pos_peak_idx_size=[];
            pos_peak_idxsize=[];
            min_index=[];
            chunk=[];
            eodstart=[];
            eodend=[];
            eod_start_idxs=[];
            eod_end_idxs=[];


            filename= fullfile(directory,specieslist(a),treatlist(b));   %construct the filename
            eval(sprintf('%s_%s = TDMS_getStruct(filename{1});',specieslist{a},treatlist{b}));  %create a container based on the filename to save data to, and open the data
            eval(sprintf('eodwave=%s_%s.Untitled.Dev1_ai0.data;',specieslist{a}, treatlist{b})); %put the data in a variable called eodwave
            eval(sprintf('eodwaveinfo=%s_%s.Props;',specieslist{a}, treatlist{b})); %put the data in a variable called eodwave
            
            j1=[1,-0.98];j2=[1,-1];
            eodwave=filtfilt(j2,j1,eodwave);
            
            %baseline=mean(eodwave(1:50));
            %eodwave=eodwave-baseline;
            
            recordlength=min(length(eodwave),1000000);
            if recordlength < 1000000
                fprintf('   Step %s.%s: Analysis of %s aborted, less than 1000000 points in recording \n',num2str(a),num2str(b), treatlist{b});  %Step 1.1.1 %Analyis of BBRACH CRISPR file1.tdms
                continue
                
            end
            plot(eodwave(1:recordlength));
            gcf;
            figure(eodpreview);

            commandwindow;
            threshold=input(prompt);
            %This was on in jasons code gain=input(prompt2);
            eodwindowsize=2000; %this was the code for Jason instead of 2000
            % which is 0.02*100,000: input(prompt3)*100000;
            %t2=input(prompt4);
            

            

            %valus of 3,000,000 is 30 seconds (30 * 100,000 samples per
            %second)
            [pos_peaks,pos_peak_idx]=findpeaks(eodwave(1:3000000),'MinPeakHeight',threshold,'MinPeakDistance',100);
            %[neg_peaks,neg_peak_idx]=findpeaks(-eodwave(1:3000000),'MinPeakHeight',threshold,'MinPeakDistance',100);
            
            
            [~,pos_peak_idxsize]=size(pos_peak_idx);
            pos_peak_idx=pos_peak_idx(2:(pos_peak_idxsize-1));
            pos_peaks=pos_peaks(2:(pos_peak_idxsize-1));
            [~,pos_peak_idx_size]=size(pos_peak_idx);
            
            %[~,neg_peak_idxsize]=size(neg_peak_idx);
            %neg_peak_idx=neg_peak_idx(2:(neg_peak_idxsize-1));
            %neg_peaks=neg_peaks(2:(neg_peak_idxsize-1));
            %[~,neg_peak_idx_size]=size(neg_peak_idx);
            %neg_peaks=abs(neg_peaks);
            
            %new code finds the negative peak nearby each positive
            %peak.  use 0.02 seconds for eodwindowsize for b. brachyistius
            neg_peaks=[];
            neg_peak_idx=[];
            for e=1:length(pos_peak_idx)
                [neg_peaks(e),local_idx]=min(eodwave((pos_peak_idx(e)-eodwindowsize):(pos_peak_idx(e)+eodwindowsize)));
                neg_peak_idx(e)=local_idx+(pos_peak_idx(e)-eodwindowsize);
            end
                
                
            
            all_peaks=pos_peaks+abs(neg_peaks);
            
            all_peaks_mean=mean(all_peaks);
            all_peaks_std=std(all_peaks);

            treatment=treatlist(b);
            fishname=strcat('fish_',num2str(b));
           
            raw_data.species=specieslist(a);
            raw_data.treat=treatlist(b);
            raw_data.filename=treatlist(b);
            raw_data.fish_name=fishname;
            raw_data.amplitude=all_peaks_mean;
            raw_data.peaks=all_peaks;
            raw_data.amplitidue_std=all_peaks_std;
            raw_data.threshold=threshold;
            %This was on in jasons code: raw_data.gain=gain;
            %raw_data.t1=t1;
            %raw_data.t2=t2;
            
            
            final_data=[final_data; raw_data];
        %end
        
    end
end
fprintf('\nAll Analysis Complete!  Enjoy your data!\n');
clearvars -except final_data

bp_data=horzcat(final_data.peaks);
[numvars,~]=size(final_data);
grp=[];
for i=1:numvars
    [~,npeaks]=size(final_data(i).peaks);
    grp=[grp,ones(1,npeaks)*i];
end
boxplot(bp_data,grp)

 for i=1:length(final_data)
     writematrix(final_data(i).peaks',char(fullfile('C:/Users/Savvas Constantinou/Desktop/Rawvoltages',strcat(final_data(i).filename,'_raw_voltages.csv'))))
 end