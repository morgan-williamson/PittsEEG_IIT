%% Source Waveform Data Concatenation & Preparation 

% It's important to note that the SPM extraction changes from location (32,
% -40, -18) to (32,-40,-20) after the first 50 trials

% Inputs:   all source waveform trials for faces (200) & random (200) at each
%           channel location (see 'coords'). 'all_waveforms' folder must be in PATH.
%           currently the extracted source activity waveforms are assumed
%           to be at 500Hz for a 1s period [-100ms to +900ms around
%           stimulus onset]

% Outputs:  concatenated matrices of raw source waveforms such that first
%           half in dimension 1 are 'face' trials and second half are
%           'random' trials from phase 3. 
%           - source_raw_phase3 
%           - source_bin_p3_df/ms = binarised above according to either
%           differentiation (df) or median split (ms)
%           - source_bin_2/5/10ds_p3_df/ms = downsampled by averaging 
%           - source_bin_2/5/10ds_g3_p3_df/ms = ds & grouped over 3 trials.
%
% Next part of pipeline: test the KL divergence of the probability
% distributions over the frequencies of states per timestep between
% face/random for each of these state series, and/or;
%
% Feed each of the state series into source_TPM_generator.m and generate 
% TPM estimates to then feed into PyPhi for IIT calculus. 
%
% Author: Morgan Williamson
% Last edited: 24/02/2020
% To update: create one last ds with 25+ms sample intervals for
% differention as we have used a 40Hz high pass filter.
%% Generation of 500Hz state series (median split/differentiation)

% First create the matrices to fill; the columns 1-200 are from face trials (31),
%                                    201-400 are from random trials (32).

cond = {'faces','random'};
coords = {'30 -43 -18','-32 -44 -22','5  35 -27','-6  20 -27'};
source_raw_ch1 = zeros(400,501);
source_binarised_ch1_df = zeros(400,500); source_binarised_ch1_ms = zeros(400,500);
source_raw_ch2 = zeros(400,501);
source_binarised_ch2_df = zeros(400,500); source_binarised_ch2_ms = zeros(400,500);
source_raw_ch3 = zeros(400,501);
source_binarised_ch3_df = zeros(400,500); source_binarised_ch3_ms = zeros(400,500);
source_raw_ch4 = zeros(400,501);
source_binarised_ch4_df = zeros(400,500); source_binarised_ch4_ms = zeros(400,500);

for ch = 1:4        % tracks each of the 4 'channels' provided by Elise from SPM waveform extraction.
    for c = 1:2     % tracks condition face / random (when concatenated 1st half is face and 2nd half is random)
        for i = 1:200
            strcoords = coords{ch};
            filename = strcat('Source_Waveform_at_',strcoords,'_for_bfmSe_',cond{c},'_fspmeeg_7MJN_trial_',num2str(i));
            load(filename,'waveform') %loads the source extracted waveform provided by Elise 
            if c == 1
                add_200 = 0;
            else
                add_200 = 200; % add random trials to colomns 201-400 after face trials 1-200.
            end
            if ch == 1
                source_raw_ch1(i+add_200,:) = waveform; 
            elseif ch ==2
                source_raw_ch2(i+add_200,:) = waveform; 
            elseif ch == 3
                source_raw_ch3(i+add_200,:) = waveform; 
            else 
                source_raw_ch4(i+add_200,:) = waveform; 
            end
            
            for j = 1:500 % Going through individual timesteps to create binarised state series. 
                
                if ch == 1
                    if source_raw_ch1(i+add_200,j+1) > source_raw_ch1(i+add_200,j) % condition for state = 1 under differentiation
                        source_binarised_ch1_df(i+add_200,j) = 1;
                    end
                    if source_raw_ch1(i+add_200,j) > mean(source_raw_ch1(i,:))     % condition for state = 2 under median split
                        source_binarised_ch1_ms(i+add_200,j) = 1;
                    end 
                elseif ch ==2                    
                    if source_raw_ch2(i+add_200,j+1) > source_raw_ch2(i+add_200,j)
                        source_binarised_ch2_df(i+add_200,j) = 1;
                    end 
                    if source_raw_ch2(i+add_200,j) > mean(source_raw_ch2(i,:))
                        source_binarised_ch2_ms(i+add_200,j) = 1;
                    end
                elseif ch == 3                    
                    if source_raw_ch3(i+add_200,j+1) > source_raw_ch3(i+add_200,j)
                        source_binarised_ch3_df(i+add_200,j) = 1;
                    end 
                    if source_raw_ch3(i+add_200,j) > mean(source_raw_ch3(i,:))
                        source_binarised_ch3_ms(i+add_200,j) = 1;
                    end
                else 
                    if source_raw_ch4(i+add_200,j+1) > source_raw_ch4(i+add_200,j)
                        source_binarised_ch4_df(i+add_200,j) = 1;
                    end 
                    if source_raw_ch4(i+add_200,j) > mean(source_raw_ch4(i,:))
                        source_binarised_ch4_ms(i+add_200,j) = 1;
                    end
                end
                 
            end
            clear('waveform')
        end 
    end
   
end

source_bin_p3_df = zeros(400,500,4); %creation of concatenated binarised files.
source_bin_p3_ms = zeros(400,500,4);

source_bin_p3_df(:,:,1) = source_binarised_ch1_df;
source_bin_phase3_ms(:,:,1) = source_binarised_ch1_ms;
source_bin_p3_df(:,:,2) = source_binarised_ch2_df;
source_bin_p3_ms(:,:,2) = source_binarised_ch2_ms;
source_bin_p3_df(:,:,3) = source_binarised_ch3_df;
source_bin_p3_ms(:,:,3) = source_binarised_ch3_ms;
source_bin_p3_df(:,:,4) = source_binarised_ch4_df;
source_bin_p3_ms(:,:,4) = source_binarised_ch4_ms;

source_raw_p3 = zeros(400,501,4);
source_raw_p3(:,:,1) = source_raw_ch1;
source_raw_p3(:,:,2) = source_raw_ch2;
source_raw_p3(:,:,3) = source_raw_ch3;
source_raw_p3(:,:,4) = source_raw_ch4;


%% Downsampled data 
% I downsampled by factors of 2, 5, 10 to 4ms, 10ms and 20ms intervals
% respectively. 

source_raw_2ds_p3 = zeros(400,251,4);
source_bin_2ds_p3_df = zeros(400,250,4);
source_bin_2ds_p3_ms = zeros(400,250,4);

for i = 1:250 % 4ms intervals / 250Hz raw data
    for trial = 1:400
        for channel = 1:4
            source_raw_2ds_p3(trial,i,channel) = mean(source_raw_p3(trial,(2*(i-1)+1):(2*i),channel));
        end 
    end
end

source_raw_2ds_p3(:,251,:) = source_raw_p3(:,501,:); % I need the last point to do the differentiation and have 250 points to properly compare ms to df

for i = 1:250 % 4ms intervals / 250Hz binarised state series
    for j = 1:400
        for k = 1:4
            if source_raw_2ds_p3(j,i+1,k) > source_raw_2ds_p3(j,i,k)
                source_bin_2ds_p3_df(j,i,k) = 1;
            end
            if source_raw_2ds_p3(j,i,k) > mean(source_raw_2ds_p3(j,:,k))
                source_bin_2ds_p3_ms(j,i,k) = 1;
            end
        end
    end
end

source_raw_5ds_p3 = zeros(400,101,4);

for i = 1:100  % 10ms intervals / 100Hz raw data
    for trial = 1:400
        for channel = 1:4
            source_raw_5ds_p3(trial,i,channel) = mean(source_raw_p3(trial,(5*(i-1)+1):(5*i),channel));
        end 
    end
end
source_raw_5ds_p3(:,101,:) = source_raw_p3(:,501,:);


source_bin_5ds_p3_df = zeros(400,100,4);
source_bin_5ds_p3_ms = zeros(400,100,4);

for i = 1:100 % 10ms intervals / 100Hz binarised data
    for j = 1:400
        for k = 1:4
            if source_raw_5ds_p3(j,i+1,k) > source_raw_5ds_p3(j,i,k)
                source_bin_5ds_p3_df(j,i,k) = 1;
            end
            if source_raw_5ds_p3(j,i,k) > mean(source_raw_5ds_p3(j,:,k))
                source_bin_5ds_p3_ms(j,i,k) = 1;
            end
        end
    end
end

source_raw_10ds_p3 = zeros(400,51,4);

for i = 1:50 % 20ms intervals / 50Hz raw data
    for trial = 1:400
        for channel = 1:4
            source_raw_10ds_p3(trial,i,channel) = mean(source_raw_p3(trial,(10*(i-1)+1):(10*i),channel));
        end 
    end
end

source_raw_10ds_p3(:,51,:) = source_raw_p3(:,501,:);


source_bin_10ds_p3_ms = zeros(400,50,4);
source_bin_10ds_p3_df = zeros(400,50,4);
for i = 1:50 % 20ms intervals / 50Hz binarised data
    for j = 1:400
        for k = 1:4            
            if source_raw_10ds_p3(j,i+1,k) > source_raw_10ds_p3(j,i,k)
                source_bin_10ds_p3_df(j,i,k) = 1;
            end
            if source_raw_10ds_p3(j,i,k) > mean(source_raw_10ds_p3(j,:,k))
                source_bin_10ds_p3_ms(j,i,k) = 1;
            end
        end
    end
end

%% Grouped (n=3) across trials & downsampled 
% ERP is averaged across 3 consecutive trials to smooth out intertrial
% variation/noise. No intertrial dependencies (each trial is involved in only one average).
% downsampling is then also performed

source_raw_g3_p3 = zeros(132,501,4);
source_raw_2ds_g3_p3 = zeros(132,251,4);
source_raw_5ds_g3_p3 = zeros(132,101,4);
source_raw_10ds_g3_p3 = zeros(132,51,4);

source_bin_g3_p3_df = zeros(132,500,4);
source_bin_g3_p3_ms = zeros(132,500,4);
source_bin_2ds_g3_p3_ms = zeros(132,250,4);
source_bin_2ds_g3_p3_df = zeros(132,250,4);
source_bin_5ds_g3_p3_ms = zeros(132,100,4);
source_bin_5ds_g3_p3_df = zeros(132,100,4);
source_bin_10ds_g3_p3_ms = zeros(132,50,4);
source_bin_10ds_g3_p3_df = zeros(132,50,4);




for i = 1:66 % question: does it matter in which order I do the averaging? Should I average ERPs across trials first,
             %           or should I instead downsample average across timesteps first? 
             % answer: it makes a small difference to the raw files but no difference to the binarised files            
    source_raw_g3_p3(i,:,:) = mean(source_raw_p3((3*(i-1)+1):3*i,:,:),1);
    source_raw_g3_p3(i+66,:,:) = mean(source_raw_p3((3*(i-1)+201):200+3*i,:,:),1);

    source_raw_2ds_g3_p3(i,:,:) = mean(source_raw_2ds_p3((3*(i-1)+1):3*i,:,:),1);
    source_raw_2ds_g3_p3(i+66,:,:) = mean(source_raw_2ds_p3((3*(i-1)+201):200+3*i,:,:),1);
    
    source_raw_5ds_g3_p3(i,:,:) = mean(source_raw_5ds_p3((3*(i-1)+1):3*i,:,:),1);
    source_raw_5ds_g3_p3(i+66,:,:) = mean(source_raw_5ds_p3((3*(i-1)+201):200+3*i,:,:),1);
    
    source_raw_10ds_g3_p3(i,:,:) = mean(source_raw_10ds_p3((3*(i-1)+1):3*i,:,:),1);
    source_raw_10ds_g3_p3(i+66,:,:) = mean(source_raw_10ds_p3((3*(i-1)+201):200+3*i,:,:),1);
    
end


for i = 1:66
    for k = 1:4
        for j = 1:500
            if source_raw_g3_p3(i,j+1,k) > source_raw_g3_p3(i,j,k)
                source_bin_g3_p3_df(i,j,k) = 1;
            end
            if source_raw_g3_p3(i,j,k) > mean(source_raw_g3_p3(i,:,k))
                source_bin_g3_p3_ms(i,j,k) = 1;
            end            
        end
        for j = 1:250
            if source_raw_2ds_g3_p3(i,j+1,k) > source_raw_2ds_g3_p3(i,j,k)
                source_bin_2ds_g3_p3_df(i,j,k) = 1;
            end
            if source_raw_2ds_g3_p3(i,j,k) > mean(source_raw_2ds_g3_p3(i,:,k))
                source_bin_2ds_g3_p3_ms(i,j,k) = 1;
            end            
        end
        for j = 1:100
            if source_raw_5ds_g3_p3(i,j+1,k) > source_raw_5ds_g3_p3(i,j,k)
                source_bin_5ds_g3_p3_df(i,j,k) = 1;
            end
            if source_raw_5ds_g3_p3(i,j,k) > mean(source_raw_5ds_g3_p3(i,:,k))
                source_bin_5ds_g3_p3_ms(i,j,k) = 1;
            end            
        end
        for j = 1:50
            
            if source_raw_10ds_g3_p3(i,j+1,k) > source_raw_10ds_g3_p3(i,j,k)
                source_bin_10ds_g3_p3_df(i,j,k) = 1;
            end
            if source_raw_10ds_g3_p3(i,j,k) > mean(source_raw_10ds_g3_p3(i,:,k))
                source_bin_10ds_g3_p3_ms(i,j,k) = 1;
            end        
        end
    end
end


save('source_bin_data_complete_4ch', 'source_bin_p3_df','source_bin_p3_ms', 'source_bin_g3_p3_df','source_bin_g3_p3_ms',... 
    'source_bin_2ds_p3_df', 'source_bin_2ds_p3_ms', 'source_bin_5ds_p3_df', 'source_bin_5ds_p3_ms', ...
    'source_bin_10ds_p3_df','source_bin_10ds_p3_ms', 'source_bin_2ds_g3_p3_df','source_bin_2ds_g3_p3_ms', ...
    'source_bin_5ds_g3_p3_df', 'source_bin_5ds_g3_p3_ms', 'source_bin_10ds_g3_p3_df', 'source_bin_10ds_g3_p3_ms')
%note: no signal for trial 9,130,349