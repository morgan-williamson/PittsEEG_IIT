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
%% DESCRIPTION OF SCRIPT
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------
% ** Note: Helpful IDs -> 11 - Phase 1 Face, 12 - Phase 1 Random, 21 - Phase 2 Face, 22 - Phase 2 Random, 31 - Phase 3 Face, 32 - Phase 3 Random **
% Building State-by-node TPM (Variable Name: TpM_11_A, TpM_12, TpM_21, TpM_22, TpM_31, TpM_32)
% Counts in for each a(ij) [i: row number, j: column number] is needed for weighted-average phi computation
      %                        SBN TPM
      %   aX = total occurences of sX accross trials of a condition
      %   a12 = total transitions from s1 to s2 accross trials of a
      %   condition 
      %   a1A = total transitions from s1 to any state with A on (s2 or s4)
      %
      % | State at t | Pr(A = ON) at t+delta | Pr(B = ON) at t+delta |
      % | (0,0,0)- s1  |         a1A/a1        |        a1B/a1         |
      % | (1,0,0)- s2  |         a2A/a2        |        a2B/a2         |
      % | (0,1,0)- s3  |         a3A/a3        |        a3B/a3         |
      % | (1,1,0)- s4  |         a4A/a4        |        a4B/a4         |
      % | (0,0,1)- s5  |         a1A/a1        |        a1B/a1         |
      % | (1,0,1)- s6  |         a2A/a2        |        a2B/a2         |
      % | (0,1,1)- s7  |         a3A/a3        |        a3B/a3         |
      % | (1,1,1)- s8  |         a4A/a4        |        a4B/a4         |
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------
% Inputs: samp11 and size11
% Outputs: TpM_p11, the TPM across all trials of Phase 11 and FOcombos, the
% combinatorial pairings of channels from frontal and occipital areas.
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------
% Authors: Morgan Williamson
% Last Date Modified: 20/01/20
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Finding all combos of {PFC}U{OCC}
% frontal = [22 9 23 34 20 8 21 35 33];
% occipital = [24 12 6 13 2 1 10 5 11];
% frontoccipital = [24 12 6 13 2 1 10 5 11 22 9 23 34 20 8 21 35 33];
% FOcombos = nchoosek(occipital, 4);
%% Holds binarised data for every possible 2 channel pair
xm = {'ms','df'};
dm = {'','2ds','5ds','10ds'};
g3 = {'','g3_'};
for downsample = 2:4 %already done all the stuff invovling no downsampling
    for x = 1:2
        for g = 1:2
            if downsample == 1
                %TPM_comp = source_bin; %forget about it, it's been done.
                %timesteps = 499;
            elseif downsample == 2
                if x == 1
                    if g == 1
                        TPM_comp = source_bin_2ds_p3_ms;
                    else
                        TPM_comp = source_bin_2ds_g3_p3_ms;
                    end
                else 
                    if g == 1
                        TPM_comp = source_bin_2ds_p3_df;
                    else
                        TPM_comp = source_bin_2ds_g3_p3_df;
                    end
                end
                timesteps = 249;
            elseif downsample == 3
                if x == 1
                    if g == 1
                        TPM_comp = source_bin_5ds_p3_ms;
                    else
                        TPM_comp = source_bin_5ds_g3_p3_ms;
                    end
                else 
                    if g == 1
                        TPM_comp = source_bin_5ds_p3_df;
                    else
                        TPM_comp = source_bin_5ds_g3_p3_df;
                    end
                end
                timesteps = 99;
            elseif downsample == 4
                if x == 1
                    if g == 1
                        TPM_comp = source_bin_10ds_p3_ms;
                    else
                        TPM_comp = source_bin_10ds_g3_p3_ms;
                    end
                else 
                    if g == 1
                        TPM_comp = source_bin_10ds_p3_df;
                    else
                        TPM_comp = source_bin_10ds_g3_p3_df;
                    end
                end
                timesteps = 49;
            end 

    %% TpM_p11 is accross all trials of condition 11 (assume the TPM does not change over time, individuals, or trials)
    % a11_11 is the number of times we go from state 1 to state 1.
    % a1A_11 is the number of times we go from state 1 to A being on (s2 or s4)
    % accross all trials of condition 11. We're going immediately to the sbn TPMs 
    % so ignore a11_11 etc. a1_11 is the number of time state 1 occurs.
            counter = 0;
            TpM_4ch_overall = zeros(16,4);
            TpM_4ch_counter = zeros(16,4);


            for t = 1:timesteps % indexes through the 500 timesteps (-100 to 900ms)
                for i = 1:size(TPM_comp,1) %indexes through all of the trials
                    for j = 1:16 %indexes through all possible states of 4 binary channels
                        if isequal(states4ch_aligned(1,j,:), TPM_comp(i,t,:)) %checks which state we're currently in
                            counter = counter + 1;
                            for k = 1:4
                                if TPM_comp(i,t+1,k) == 1
                                    TpM_4ch_counter(j,k) = TpM_4ch_counter(j,k) + 1; % collects the transition count                          
                                end
                            end
                        end
                    end 
                end
            end 

            TpM_4ch_overall(:,:) = TpM_4ch_counter(:,:) ./ counter ;
            save(['source_TPM_4ch_' dm{downsample} '_' g3{g} xm{x}],'TpM_4ch_overall')
        end
    end
end
    %counter = 0;


