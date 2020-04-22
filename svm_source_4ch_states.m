%% Creation of Training and Testing Data 100Ms by 50ms per REGION
% Analysis Parameters: 2ms sampling, 2ch IIT (153), Phase 31/32, 100ms windows
% shifting 50ms (13)
% Author: Morgan Williamson 
% Last Edited: 21/03/20
% -------------------------------------------------------------------------
% Here we've decided to take 70% of the data to train, and the remaining
% 30% will be for testing purposes. 
%
% Inputs: source_bin_(2/5/10ds)(_g3)_p3_ms/df -  [153,3,213,13] epochised IIS data.
%         
%         
% Outputs: AccPer_Region - [13,3,100] stores accuracies of models fed
%                        instance data of IIS from all channel pairs in
%                        each region. 1 = F (36x3), 2 = O (36x3), 3 = FO
%                        (81x3 instance properties). Cross-validation is
%                        performed n = 100 times. 
% This code calls svmtrain and svmpredict from LIBSVM.
%
 
%% Preparation of files for SVM analysis.
% See LIBSVM matlab README documentation for more details. For m instances
% with n properties and k = 2 classifications, we need to create an [m 1]
% sized labelling file with each point containing the classification of
% that instance (1 = Face / 2 = Random). We do this for training and
% testing separately as they involve a separate number of trials. 

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

features = 16;
states4ch_aligned = zeros(1,16,4);
states4ch_aligned(1,:,:) = states4ch;

glabel = {'','_g3'};
xlabel = {'_df','_ms'};
dslabel = {'','_2ds','_5ds','_10ds'};

for ds = 1:4
    for x = 1:2
        for g = 1:2
            
            if g == 1
                svm_trials = 200;
            elseif g == 2
                svm_trials = 66;
            end

            train_num = round(svm_trials*0.7);                 
            test_num = (svm_trials - train_num);
            training_label_vector = [];
            testing_label_vector = [];

            for i = 1:(train_num*2)
                if i <= train_num                       %first half of data points are 31 (face), second half are 32 random 
                    training_label_vector = [training_label_vector 1];
                else
                    training_label_vector = [training_label_vector 2];
                end
            end
            training_label_vector = training_label_vector.' ;

            for i = 1:(test_num*2)
                if i <= test_num
                    testing_label_vector = [testing_label_vector 1];
                else
                    testing_label_vector = [testing_label_vector 2];
                end
            end
            testing_label_vector = testing_label_vector.' ;

            testing_instance_matrix = zeros(test_num*2,features); %trials x instance properties (15 x 126)
            training_instance_matrix = zeros(train_num*2,features);
            
            
            params = [100,23; 200,21; 400,16; 600,11];
            for q = 1:4
                epoch = params(q,1);
                epoch_num = params(q,2);
                
                Acc_states_source = zeros(epoch_num,100);%%
                for n = 1:100 
                    train_choice = randperm(svm_trials,train_num); % = a random set of 149 trials from the 213, used for SVM training

                    % For now the choice is that the same random permutation will be applied to
                    % 31 and 32 sets. Perhaps this matters, and the random permutation
                    % applied should be independent. 

                    train31 = zeros(features, train_num, epoch_num);
                    test31 = zeros(features, test_num, epoch_num);
                    train32 = zeros(features, train_num, epoch_num);
                    test32 = zeros(features, test_num, epoch_num);

                    test_count = 1;
                    add_to_test = 1;
                    
                    [e_1,e_2] = svm_state_freq_epoch_generator(eval(['source_bin' dslabel{ds} glabel{g} '_p3' xlabel{x}]),epoch_num,ds,states4ch_aligned);


                    for i = 1:svm_trials % = 100 through trials
                        for j = 1:train_num % for each of the training trials (140)
                            if i == train_choice(j) % is this trial one of out random training 70%? 
                                
                                    train31(:,j,:) = e_1(:,i,:);
                                    train32(:,j,:) = e_2(:,i,:);

                                add_to_test = 0;
                                break
                            end
                        end
                        if add_to_test == 1
                            
                                test31(:,test_count,:) = e_1(:,i,:);
                                test32(:,test_count,:) = e_2(:,i,:);
                            
                            test_count = test_count + 1;
                        else
                            add_to_test = 1;
                        end
                    end %assigned 32 to test and train arrays
                    clear('train_choice','add_to_test','test_count');

                    for k = 1:epoch_num
                        for i = 1:train_num*2
                            if i <= train_num
                                training_instance_matrix(i,:) = train31(:,i,k);
                            else 
                                training_instance_matrix(i,:) = train32(:,(i-train_num),k);
                            end
                        end
                        for i = 1:test_num*2
                            if i <= test_num
                                testing_instance_matrix(i,:) = test31(:,i,k);
                            else 
                                testing_instance_matrix(i,:) = test32(:,(i-test_num),k);
                            end
                        end

                            to_scale = zeros(svm_trials*2,16);
                            to_scale(1:(train_num*2),:) = training_instance_matrix(:,:);
                            to_scale((train_num*2 +1):svm_trials*2,:) = testing_instance_matrix(:,:);
                            rescaler = rescale(to_scale);
                            training_instance_matrix = rescaler(1:(train_num*2),:);
                            testing_instance_matrix = rescaler((train_num*2 +1):svm_trials*2,:);
                            clear('to_scale','rescaler');

                            model = svmtrain(training_label_vector, training_instance_matrix, '-s 0 -t 0 -c 2^(-20) -q' );
                            [predicted_label, accuracy, prob_estimates] = svmpredict(testing_label_vector, testing_instance_matrix, model);

                        Acc_states_source(k,n) = accuracy(1);
                        clear('model','predicted_label','accuracy','prob_estimates');
                        n
                    end
                end
                
                save(['acc_source_4ch_states_' num2str(epoch) 'ms_epochs' dslabel{ds} glabel{g} '_p3' xlabel{x}], 'Acc_states_source');
            end 
        end
    end
end


function [e_1, e_2] = svm_state_freq_epoch_generator(state_TS,epoch_num,ds,states_4ch_aligned)
%%% Inputs a state-series (currently 4ch), epoch_num to convert to, ds =
%%% downsampling tracker (1-4), and states_4ch_aligned (because 4 channels
%%% -> 16 states). Outputs an epochised average state series converted to
%%% state (0-16).
% state_TS is trials x times x 4 ch

trials = size(state_TS,1);
times = size(state_TS,2);
source_state_ts = zeros(trials,times,16);
states = size(states_4ch_aligned,2);

if ds == 3
    ds = 5;
elseif ds == 4
    ds = 10;
end

for i = 1:trials
    for j = 1:times
        for k = 1:states
            if isequal(states_4ch_aligned(1,k,:),state_TS(i,j,:))
                source_state_ts(i,j,k) = 1;
                break 
            end
        end
    end
end

e_1 = zeros(16, trials / 2, epoch_num); e_2 = zeros(16, trials / 2, epoch_num);
% e200_1 = zeros(16, trials / 2, 21); e200_2 = zeros(16, trials / 2, 21);
% e400_1 = zeros(16, trials / 2, 16); e400_2 = zeros(16, trials / 2, 16);
% e600_1 = zeros(16, trials / 2, 11); e600_2 = zeros(16, trials / 2, 11);

%isequal(states4ch_aligned(1,j,:), source_binarised_phase3(i,t,:))
state31_2ms = zeros(16, (trials / 2), times);
state32_2ms = zeros(16, (trials / 2), times);

for j = 1:(trials/2)
    for k = 1:times
        for i = 1:states
             
            state31_2ms(i,j,k) = source_state_ts(j,k,i);
            state32_2ms(i,j,k) = source_state_ts((j+trials/2),k,i);
        end
    end
end

if epoch_num == 23
    if ds == 1
        epoch_steps = 50;
        slide_step = 20;
    elseif ds == 2 
        epoch_steps = 25;
        slide_step = 10;
    elseif ds == 5
        epoch_steps = 10;
        slide_step = 4;
    elseif ds == 10
        epoch_steps = 5;
        slide_step = 2;
    end
elseif epoch_num == 21
    if ds == 1
        epoch_steps = 100;
        slide_step = 20;
    elseif ds == 2 
        epoch_steps = 50;
        slide_step = 10;
    elseif ds == 5
        epoch_steps = 20;
        slide_step = 4;
    elseif ds == 10
        epoch_steps = 10;
        slide_step = 2;
    end
elseif epoch_num == 16
    if ds == 1
        epoch_steps = 200;
        slide_step = 20;
    elseif ds == 2 
        epoch_steps = 100;
        slide_step = 10;
    elseif ds == 5
        epoch_steps = 40;
        slide_step = 4;
    elseif ds == 10
        epoch_steps = 20;
        slide_step = 2;
    end
elseif epoch_num == 11
    if ds == 1
        epoch_steps = 300;
        slide_step = 20;
    elseif ds == 2 
        epoch_steps = 150;
        slide_step = 10;
    elseif ds == 5
        epoch_steps = 60;
        slide_step = 4;
    elseif ds == 10
        epoch_steps = 30;
        slide_step = 2;
    end
end 

% Averaged across 100/200/400/600ms windows sliding 50ms 

for j = 1:trials/2
    for k = 1:epoch_num
        for m = 1:states
            e_1(m,j,k) = mean(state31_2ms(m,j,((k-1)*20 + 1):((k-1)*slide_step + epoch_steps)),'all');
            e_2(m,j,k) = mean(state32_2ms(m,j,((k-1)*20 + 1):((k-1)*slide_step + epoch_steps)),'all');
        end    
    end 
end
end 

