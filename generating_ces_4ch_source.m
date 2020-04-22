%% Creation of CES data for 200ms windows around each 50ms (13 windows)(For SVM)
% Analysis Parameters: 2ms sampling, 2ch IIT (126), Phase 31/32, 200ms windows
% shifting 50ms (13)
% Author: Morgan Williamson 
% Last Edited: 03/03/20
%--------------------------------------------------------------------------
% Inputs: source_binarised_phase3 & source_bin_ds_phase3 
%              - contains binarised & cut timeseries data of extracted
%              waveforms from source data at 4 channels over phase 3
%         sias - array containing the sia structs of all 16 states of the
%               4ch setup, generated from a TPM in PyPhi
%          
% Ouputs: ces31/32_2ms - arrays of IIT timeseries data where states from
%                        samp31 and a random equally sized cut of samp32
%                        (rsamp32) are replace with a vector of small phi's
%                        from the comb's cause-effect-structure in that
%                        state. 
%                        [126 channel combs, 3 ces, 213 trials, 350 time]
%         ces31/32_epochs1/2/4/600_50ms - averages of ces31/32_2ms vectors in
%               1-600ms windows that slide 50ms. Sampling is 2ms or 500Hz.
%               
%         svm_trials - number of trials to be fed from 31/32 into SVM
%         total_trials - max of 31/32 (some won't be fed to SVM)
%--------------------------------------------------------------------------
% The epochs will be fed to SVMs for training and testing. This code is
% compatible with regional_SVM_generator, as well as individual channel
% comb generators. 
% in source pipeline, feeds into svm_decoding_4ch_IIS.m
%--------------------------------------------------------------------------
%% 
%%

states4ch = [];
for a = 0:1
    for b = 0:1
        for c = 0:1
            for d = 0:1
                states4ch = [states4ch; d,c,b,a];
            end
        end
    end
end
states4ch_aligned = zeros(1,16,4);
states4ch_aligned(1,:,:) = states4ch;
%% Equalising the size of 31/32
% For training and testing purposes, we've decided to keep the number
% of trials of face vs random the same.svm_trials = min(size31(1),size32(1));
glabel = {'','_g3'};
xlabel = {'_df','_ms'};
dslabel = {'','_2ds','_5ds','_10ds'};
for ds =1:4
    for x = 1:2
        for g = 1:2
            if ds == 1 
                if x == 1
                    if g == 1
                        load('II_4ch_source_p3_TPM.mat');
                        %load('source_binarised_phase3.mat');
                        state_TS = source_bin_p3_df;
                    elseif g == 2
                        load('II_4ch_source_p3_TPM_ms.mat');
                        %load('source_binarised_phase3.mat');
                        state_TS = source_bin_p3_ms;
                    end
                end 
            elseif ds == 2
                if x == 1
                    if g == 1
                        load('II_4ch_source_2ds_p3_df.mat');
                        %load('source_bin_2ds_p3_df.mat');
                        state_TS = source_bin_2ds_p3_df;
                    elseif g == 2            
                        load('II_4ch_source_2ds_g3_p3_df.mat');
                        %load('source_bin_2ds_g3_p3_df.mat');
                        state_TS = source_bin_2ds_g3_p3_df;
                    end
                else
                    if g == 1
                        load('II_4ch_source_2ds_p3_ms.mat');
                        %load('source_bin_2ds_p3_ms.mat');
                        state_TS = source_bin_2ds_p3_ms;
                    elseif g == 2            
                        load('II_4ch_source_2ds_g3_p3_ms.mat');
                        %load('source_bin_2ds_g3_p3_ms.mat');
                        state_TS = source_bin_2ds_g3_p3_ms;
                    end
                end    
            elseif ds == 3
                if x == 1
                    if g == 1
                        load('II_4ch_source_5ds_p3_df.mat');
                        %load('source_bin_5ds_p3_df.mat');
                        state_TS = source_bin_5ds_p3_df;
                    elseif g == 2            
                        load('II_4ch_source_5ds_g3_p3_df.mat');
                        %load('source_bin_5ds_g3_p3_df.mat');
                        state_TS = source_bin_5ds_g3_p3_df;
                    end
                else
                    if g == 1
                        load('II_4ch_source_5ds_p3_ms.mat');
                        %load('source_bin_5ds_p3_ms.mat');
                        state_TS = source_bin_5ds_p3_ms; 
                    elseif g == 2            
                        load('II_4ch_source_5ds_g3_p3_ms.mat');
                        %load('source_bin_5ds_g3_p3_ms.mat');
                        state_TS = source_bin_5ds_g3_p3_ms;
                    end
                end
            elseif ds == 4
                if x == 1
                    if g == 1
                        load('II_4ch_source_10ds_p3_df.mat');
                        %load('source_bin_10ds_p3_df.mat');
                        state_TS = source_bin_10ds_p3_df;
                    elseif g == 2             
                        load('II_4ch_source_10ds_g3_p3_df.mat');
                        %load('source_bin_10ds_g3_p3_df.mat');
                        state_TS = source_bin_10ds_g3_p3_df;
                    end
                else
                    if g == 1
                        load('II_4ch_source_10ds_p3_ms.mat');
                        %load('source_bin_10ds_p3_ms.mat');
                        state_TS = source_bin_10ds_p3_ms;
                    elseif g == 2            
                        load('II_4ch_source_10ds_g3_p3_ms.mat');
                        %load('source_bin_10ds_g3_p3_ms.mat');
                        state_TS = source_bin_10ds_g3_p3_ms;
                    end
                end
            end 
            
            trials = size(state_TS,1) / 2;
            samples = size(state_TS,2);  

            %% Translation of states into cause-effect-structures
            % If you want to redefine the window parameters (200ms every 50ms for now),
            % you can keep the same ces31/32_2ms and only change the code for epochs. 
            ces31_epochs100_50ms = zeros(15,trials,23); % all of these structures will have to be changed for the 15 part IIS of 4ch ( i need to look at the sia files)
            ces32_epochs100_50ms = zeros(15,trials,23); % also 126 combinations now, and 

            ces31_2ms = zeros(15,trials,samples);   ces31_epochs200_50ms = zeros(15,trials,21);
            ces32_2ms = zeros(15,trials,samples);   ces32_epochs200_50ms = zeros(15,trials,21);

            ces31_epochs400_50ms = zeros(15,trials,16);
            ces32_epochs400_50ms = zeros(15,trials,16);


            ces31_epochs600_50ms = zeros(15,trials,11);
            ces32_epochs600_50ms = zeros(15,trials,11);
            

            for j = 1:trials
                for k = 1:samples
                    for state = 1:16
                        if isequal(states4ch_aligned(1,state,:), state_TS(j,k,:))
                            for z = 1:15
                                ces31_2ms(z,j,k) = sias{1,(states4ch(state,1)+1),(states4ch(state,2)+1),(states4ch(state,3)+1),(states4ch(state,4)+1)}.ces.concepts{1,z}.phi;
                            end 
                        end
                        if isequal(states4ch_aligned(1,state,:), state_TS(j+trials,k,:))
                            for z = 1:15
                                ces32_2ms(z,j,k) = sias{1,(states4ch(state,1)+1),(states4ch(state,2)+1),(states4ch(state,3)+1),(states4ch(state,4)+1)}.ces.concepts{1,z}.phi;
                            end
                        end
                    end
                end
            end


            % Averaged across 100/200/400/600ms windows sliding 50ms 
            % how do I split this so that the epochs are the same or at
            % least comparable?
            
            for j = 1:trials
                for k = 1:23
                    for m = 1:15
                        if ds == 1
                            ces31_epochs100_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*20 + 1):((k-1)*20 + 50)),'all');
                            ces32_epochs100_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*20 + 1):((k-1)*20 + 50)),'all');
                        elseif ds == 2
                            ces31_epochs100_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*10+1):((k-1)*10 + 25)),'all');
                            ces32_epochs100_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*10+1):((k-1)*10 + 25)),'all');
                        elseif ds == 3
                            ces31_epochs100_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*4+1):((k-1)*4 + 10)),'all');
                            ces32_epochs100_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*4+1):((k-1)*4 + 10)),'all');
                        elseif ds == 4
                            ces31_epochs100_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*2+1):((k-1)*2 + 5)),'all');
                            ces32_epochs100_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*2+1):((k-1)*2 + 5)),'all');
                        end
                    end
                end
                for k = 1:21
                    for m = 1:15
                        if ds == 1
                            ces31_epochs200_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*20 + 1):((k-1)*20 + 100)),'all');
                            ces32_epochs200_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*20 + 1):((k-1)*20 + 100)),'all');
                        elseif ds == 2
                            ces31_epochs200_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*10+1):((k-1)*10 + 50)),'all');
                            ces32_epochs200_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*10+1):((k-1)*10 + 50)),'all');
                        elseif ds == 3
                            ces31_epochs200_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*4+1):((k-1)*4 + 20)),'all');
                            ces32_epochs200_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*4+1):((k-1)*4 + 20)),'all');
                        elseif ds == 4
                            ces31_epochs200_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*2+1):((k-1)*2 + 10)),'all');
                            ces32_epochs200_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*2+1):((k-1)*2 + 10)),'all');
                        end
                        
                    end
                end
                for k = 1:16
                    for m = 1:15
                        if ds == 1
                            ces31_epochs400_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*20 + 1):((k-1)*20+200)),'all');
                            ces32_epochs400_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*20 + 1):((k-1)*20+200)),'all');
                        elseif ds == 2
                            ces31_epochs400_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*10+1):((k-1)*10 + 100)),'all');
                            ces32_epochs400_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*10+1):((k-1)*10 + 100)),'all');
                        elseif ds == 3
                            ces31_epochs400_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*4+1):((k-1)*4 + 40)),'all');
                            ces32_epochs400_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*4+1):((k-1)*4 + 40)),'all');
                        elseif ds == 4
                            ces31_epochs400_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*2+1):((k-1)*2 + 20)),'all');
                            ces32_epochs400_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*2+1):((k-1)*2 + 20)),'all');
                        end
                        
                    end
                end
                 for k = 1:11
                     for m = 1:15
                        if ds == 1
                            ces31_epochs600_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*20 + 1):((k-1)*20 + 300)),'all');
                            ces32_epochs600_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*20 + 1):((k-1)*20 + 300)),'all');
                        elseif ds == 2
                            ces31_epochs600_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*10+1):((k-1)*10 + 150)),'all');
                            ces32_epochs600_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*10+1):((k-1)*10 + 150)),'all');
                        elseif ds == 3
                            ces31_epochs600_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*4+1):((k-1)*4 + 60)),'all');
                            ces32_epochs600_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*4+1):((k-1)*4 + 60)),'all');
                        elseif ds == 4
                            ces31_epochs600_50ms(m,j,k) = mean(ces31_2ms(m,j,((k-1)*2+1):((k-1)*2 + 30)),'all');
                            ces32_epochs600_50ms(m,j,k) = mean(ces32_2ms(m,j,((k-1)*2+1):((k-1)*2 + 30)),'all');
                        end
                         
                     end
                 end     
            end 
            save(['source_4ch_ces' dslabel{ds} glabel{g} '_p3' xlabel{x}],'ces31_2ms','ces32_2ms','ces31_epochs100_50ms','ces31_epochs200_50ms','ces31_epochs400_50ms','ces31_epochs600_50ms','ces32_epochs100_50ms','ces32_epochs200_50ms','ces32_epochs400_50ms','ces32_epochs600_50ms')
        end
    end
end
            