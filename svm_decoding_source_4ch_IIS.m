%% Creation of Training and Testing Data of IIS epochs for LIBSVM
% Analysis Parameters: 2ms sampling, source-localised 4ch IIT, Phase 31/32, 1/2/4/600ms windows
% shifting 50ms (19/17/13/9 epochs total)
% Author: Morgan Williamson 
% Last Edited: 03/03/20
% -------------------------------------------------------------------------
% Here I've decided to take 70% of the data to set params by crosss-validation / train, and the remaining
% 30% will be for testing purposes. 
%
% Inputs: ces31/32_epochs1/2/4/600_50ms -  epochised IIS data.
%         LIBSVM training & predicting functions (see LIBSVM documentation)   
%         svm_trials - 200 from generating_ces_data files
% Outputs: Acc_4ch_source - [epoch_num,100] stores accuracies of models fed
%                        instance data of IIS 
% This code calls svmtrain and svmpredict from LIBSVM.
% feeds to plotter and from 
 
%% Preparation of files for SVM analysis.
% See LIBSVM matlab README documentation for more details. For m instances
% with n properties and k = 2 classifications, we need to create an [m 1]
% sized labelling file with each point containing the classification of
% that instance (1 = Face / 2 = Random). We do this for training and
% testing separately as they involve a separate number of trials. 

param_cv = 0; %determines if the script will attempt parameter selection via grid-search cross validation (=1)
glabel = {'','_g3'};
xlabel = {'_df','_ms'};
dslabel = {'','_2ds','_5ds','_10ds'};
for ds = 1:4
    for x = 1:2
        for g = 1:2
            load(['source_4ch_ces' dslabel{ds} glabel{g} '_p3' xlabel{x}]);
            if g == 1
                svm_trials = 200;
            elseif g == 2
                svm_trials = 66;
            end
            
            train_num = round(svm_trials*0.7);     test_num = (svm_trials - train_num);
            training_label_vector = []; testing_label_vector = [];      param_label_vector = [];
            
            
            if param_cv == 1
             %alters whether I'm doing cross-validation to select the best cost (takes way more time)
                for i = 1:(train_num)
                    if i <= train_num/2 %first half of data points are 31 (face), second half are 32 random 
                        training_label_vector = [training_label_vector 1];
                        param_label_vector = [param_label_vector 1]; %splits train set into param selection set and train set
                    else
                        training_label_vector = [training_label_vector 2];
                        param_label_vector = [param_label_vector 2];
                    end
                end
            else % if we're not doing parameter selection and just setting cost = 2^-20
                for i = 1:(train_num*2)
                    if i <= train_num %first half of data points are 31 (face), second half are 32 random 
                        training_label_vector = [training_label_vector 1];            
                    else
                        training_label_vector = [training_label_vector 2];           
                    end
                end
            end


            for i = 1:(test_num*2)
                if i <= test_num
                    testing_label_vector = [testing_label_vector 1];
                else
                    testing_label_vector = [testing_label_vector 2];
                end
            end

            training_label_vector = training_label_vector.'; param_label_vector = param_label_vector.'; testing_label_vector = testing_label_vector.';
            
            params = [100,23; 200,21; 400,16; 600,11];
            for q = 1:4
                epoch = params(q,1);
                epoch_num = params(q,2);

                Acc_4ch_source = zeros(epoch_num,100); % The array that contains the accuracy percentages for each of the pairs(153) for each of the 200 repeats.
                 % = 140, 70% of the available SVM trials.
                 %60



                testing_instance_matrix = zeros(test_num*2,(size(ces31_2ms,1))); %trials x instance properties (15 x 126)
                if param_cv == 1   
                    training_instance_matrix = zeros(train_num,(size(ces31_2ms,1)));
                    param_instance_matrix = zeros(train_num, (size(ces31_2ms,1)));
                else 
                    training_instance_matrix = zeros(train_num*2,size(ces31_2ms,1));
                end
                if param_cv == 1
                    cv_params = zeros(100,epoch_num,2);
                end
                %%
                 for n = 1:100
                    train_choice = randperm(svm_trials,train_num); % = a random set of 149 trials from the 213, used for SVM training
                    param_select = randperm(train_num, train_num/2);
                    % For now the choice is that the same random permutation will be applied to
                    % 31 and 32 sets. Perhaps this matters, and the random permutation
                    % applied should be independent. 

                    train31 = zeros(size(ces31_2ms,1), train_num, epoch_num);
                    test31 = zeros(size(ces31_2ms,1), test_num, epoch_num);
                    train32 = zeros(size(ces31_2ms,1), train_num, epoch_num);
                    test32 = zeros(size(ces31_2ms,1), test_num, epoch_num);

                    test_count = 1;
                    add_to_test = 1;

                    for i = 1:svm_trials % = 200 through trials
                        for j = 1:train_num % for each of the training trials (140)
                            if i == train_choice(j) % is this trial one of out random training 70%? 
                                if epoch == 100
                                    train31(:,j,:) = ces31_epochs100_50ms(:,i,:);
                                    train32(:,j,:) = ces32_epochs100_50ms(:,i,:);
                                elseif epoch == 200
                                    train31(:,j,:) = ces31_epochs200_50ms(:,i,:);
                                    train32(:,j,:) = ces32_epochs200_50ms(:,i,:);
                                elseif epoch ==400
                                    train31(:,j,:) = ces31_epochs400_50ms(:,i,:);
                                    train32(:,j,:) = ces32_epochs400_50ms(:,i,:);
                                else 
                                    train31(:,j,:) = ces31_epochs600_50ms(:,i,:);
                                    train32(:,j,:) = ces32_epochs600_50ms(:,i,:);                       
                                end
                                add_to_test = 0;
                                break
                            end
                        end
                        if add_to_test == 1
                            if epoch == 100
                                test31(:,test_count,:) = ces31_epochs100_50ms(:,i,:);
                                test32(:,test_count,:) = ces32_epochs100_50ms(:,i,:);
                            elseif epoch == 200
                                test31(:,test_count,:) = ces31_epochs200_50ms(:,i,:);
                                test32(:,test_count,:) = ces32_epochs200_50ms(:,i,:);
                            elseif epoch ==400
                                test31(:,test_count,:) = ces31_epochs400_50ms(:,i,:);
                                test32(:,test_count,:) = ces32_epochs400_50ms(:,i,:);
                            else 
                                test31(:,test_count, :) = ces31_epochs600_50ms(:,i,:);
                                test32(:,test_count,:) = ces32_epochs600_50ms(:,i,:);                       
                            end

                            test_count = test_count + 1;
                        else
                            add_to_test = 1;
                        end
                    end %assigned 32 to test and train arrays
                    clear('train_choice','add_to_test','test_count');
                 %% 

                    for k = 1:epoch_num

                        if param_cv == 1 
                            add_to_train = 1;
                            train = 1;
                            for i = 1:train_num % index through the original 70% training partition
                                for z = 1:(train_num/2) %index through the 50% of the 70% to be partitioned into param selection/training sets
                                    if i == param_select(z) % is this trial one of the selected for param?
                                        param_instance_matrix(z,:) = train31(:,i,k);
                                        param_instance_matrix((z + train_num/2),:) = train32(:,i,k);
                                        add_to_train = 0;
                                        break
                                    end

                                end
                                if add_to_train == 1
                                    training_instance_matrix(train,:) = train31(:,i,k);
                                    training_instance_matrix((train + train_num/2), :) = train32(:,i,k);
                                    train = train + 1;
                                else 
                                    add_to_train = 1;
                                end            
                            end
                        else
                            for i = 1:train_num
                                training_instance_matrix(i,:) = train31(:,i,k);
                                training_instance_matrix(i+train_num,:) = train32(:,i,k);
                            end
                        end
                        for i = 1:test_num

                                testing_instance_matrix(i,:) = test31(:,i,k);
                                testing_instance_matrix(i+test_num,:) = test32(:,i,k);

                        end

                        clear('add_to_train','train');
                 %%   

                        to_scale = zeros(svm_trials*2,15);
                        if param_cv == 1
                            to_scale(1:train_num,:) = param_instance_matrix;
                            to_scale((train_num+1):(train_num*2),:) = training_instance_matrix(:,:);
                            to_scale((train_num*2 +1):svm_trials*2,:) = testing_instance_matrix(:,:);

                            rescaler = rescale(to_scale);
                            param_instance_matrix = rescaler(1:train_num,:);
                            training_instance_matrix = rescaler((train_num+1):(train_num*2),:);
                        else
                            to_scale(1:(train_num*2),:) = training_instance_matrix(:,:);
                            to_scale((train_num*2 +1):svm_trials*2,:) = testing_instance_matrix(:,:);
                            rescaler = rescale(to_scale);
                            training_instance_matrix = rescaler(1:(train_num*2),:);
                        end
                        testing_instance_matrix = rescaler((train_num*2 +1):svm_trials*2,:);
                        clear('to_scale','rescaler');

            %             # grid of parameters
            %             folds = 5;
            %             [C] = meshgrid(-20:5:20);
            %             [C,gamma] = meshgrid(-20:5:20, -20:5:20);
            % 
            %             # grid search, and cross-validation
            %             cv_acc = zeros(numel(C),1);
            %             for i=1:numel(C)
            %                 cv_acc(i) = svmtrain(param_label_vector, param_instance_matrix, sprintf('-s 0 -t 0 -c %f -v %d -q', 2^C(i), folds));
            %                 cv_acc(i) = svmtrain(param_label_vector, param_instance_matrix, sprintf('-s 2 -t 0 -c %f -g %f -v %d -q', 2^C(i), 2^gamma(i), folds));
            %             end
            % 
            %             # pair (C,gamma) with best accuracy
            %             [~,idx] = max(cv_acc);
            %             [~,idx] = max(cv_acc);
            %             best_C = 2^C(idx);

                        %%%% STORE THE PARAMETERS IF YOU CAN %%%%%

                        %best_gamma = 2^gamma(idx);

                        if param_cv == 1 %this is the cross-validation step
                            bestcv = 0;
                            for log2c = -20:10:20 %this is the parameter grid to search
                                cmd = ['-v 5 -s 0 -q  -t 0 -c ', num2str(2^log2c)];
                                cv = svmtrain(param_label_vector,param_instance_matrix, cmd);
                                if (cv >= bestcv)
                                    bestcv = cv; best_C = 2^log2c;
                                end
                            end
                            clear('cv','cmd');
                            cv_params(n,k,:) = [best_C bestcv]; % storing parameters chosen by grid-search cv             
                        else
                            best_C = 2^(-20); %otherwise we fix cost
                        end 
                        cmd1 = ['-s 0 -t 0 -c ', num2str(best_C),' -q '];
                        model = svmtrain(training_label_vector, training_instance_matrix, cmd1); % was manually -c 2^(-20) in the past
                        [predicted_label, accuracy, prob_estimates] = svmpredict(testing_label_vector, testing_instance_matrix, model);

                        Acc_4ch_source(k,n) = accuracy(1);
                        clear('model','bestcv','predicted_label','accuracy','prob_estimates');
                        
                    end
                    n
                end
            save(['acc_source_4ch_' num2str(epoch) 'ms_epochs' dslabel{ds} glabel{g} '_p3' xlabel{x}], 'Acc_4ch_source') %if cv also save the params
            end 
            
        end
    end
end