function [conditions, r2p, r2n] = GetConditions(Expt, FileType, PrefDir)

ResponseToPositive = 0;
ResponseToNegative = 0;

conditions = logical([]);
switch upper(FileType) 
    case 'DPI'
        conditions(1,:) = ones(1, length([Expt.Trials(:)]));
    
    case 'DT'
        conditions(1,:) = ones(1, length([Expt.Trials(:)]));

    case 'FLID'
        if isfield(Expt.Trials, 'sM')
            conditions(1,:) = [Expt.Trials(:).sM]==16; % flat
            conditions(2,:) = [Expt.Trials(:).sM]==17; % cylinder
        else
            conditions(1,:) = [Expt.Trials(:).flatsurf]==1;
            conditions(2,:) = [Expt.Trials(:).flatsurf]==0;
        end
    case'BDID'
        if(mean([Expt.Trials([Expt.Trials(:).Id]>0).RespDir])>0)
            ResponseToPositive = 1;
            ResponseToNegative = -1;
        else
            ResponseToPositive = -1;
            ResponseToNegative = 1;
        end
        if(PrefDir == 2)
            conditions(1,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]~=0;
            conditions(2,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]~=0;
            conditions(3,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(4,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(5,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(6,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
            % 7 and 8 are the main conditions here and BD conditions are just to see if when he was doing better the 7-8 difference was smaller or not
            conditions(7,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]~=0; 
            conditions(8,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]~=0; 
            conditions(9,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir]~=0; 
            conditions(10,:)= [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir]~=0;             
            % 11 and 12 are the the two bds near zero disparity. 
            conditions(11,:)= [Expt.Trials(:).bd] == max([Expt.Trials([Expt.Trials(:).bd]<0).bd]);
            conditions(12,:)= [Expt.Trials(:).bd] == min([Expt.Trials([Expt.Trials(:).bd]>0).bd]);
            
            conditions(15,:) = [Expt.Trials(:).RespDir]~=0; % All completed trials for normalization
        else
            conditions(1,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]~=0;
            conditions(2,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]~=0;
            conditions(3,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(4,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(5,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(6,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(7,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]~=0;
            conditions(8,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]~=0;
            conditions(9,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir]~=0; 
            conditions(10,:)= [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir]~=0; 
            conditions(11,:)= [Expt.Trials(:).bd] == min([Expt.Trials([Expt.Trials(:).bd]>0).bd]);
            conditions(12,:)= [Expt.Trials(:).bd] == max([Expt.Trials([Expt.Trials(:).bd]<0).bd]);

            conditions(15,:) = [Expt.Trials(:).RespDir]~=0; % All completed trials for normalization
        end
    
    case 'TWO'
        if(mean([Expt.Trials([Expt.Trials(:).dx]>0).RespDir])>0)
            ResponseToPositive = 1;
            ResponseToNegative = -1;
        else
            ResponseToPositive = -1;
            ResponseToNegative = 1;
        end
        if (max([Expt.Trials(:).bd])~=max([Expt.Trials(:).dx]) || min([Expt.Trials(:).bd]) ~= min([Expt.Trials(:).dx]))
            disp('max min dx bd something does not match');
        end
        if(PrefDir == 2)
            conditions(1,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]<0 & [Expt.Trials(:).RespDir]~=0;
            conditions(2,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]>0 & [Expt.Trials(:).RespDir]~=0;
            conditions(3,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(4,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative; % it was ~= until 10/13/10
            conditions(5,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive; % it was ~= until 10/13/10
            conditions(6,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(7,:) = [Expt.Trials(:).bd]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip
            conditions(8,:) = [Expt.Trials(:).bd]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip pair to compare 
            conditions(9,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]>0 & [Expt.Trials(:).RespDir]~=0;
            conditions(10,:)= [Expt.Trials(:).bd]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip
            conditions(11,:)= [Expt.Trials(:).bd]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip pair to copmare
            if (sum(conditions(7,:)) ==0 )
                conditions(7,:) = [Expt.Trials(:).bd]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip
                conditions(8,:) = [Expt.Trials(:).bd]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip pair to compare 
                conditions(10,:)= [Expt.Trials(:).bd]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip
                conditions(11,:)= [Expt.Trials(:).bd]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip pair to copmare
            end
        else
            conditions(1,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]>0 & [Expt.Trials(:).RespDir]~=0;
            conditions(2,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]<0 & [Expt.Trials(:).RespDir]~=0;
            conditions(3,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(4,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive; % it was ~= until 10/13/10
            conditions(5,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative; % it was ~= until 10/13/10
            conditions(6,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(7,:) = [Expt.Trials(:).bd]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip
            conditions(8,:) = [Expt.Trials(:).bd]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip pair to copmare
            conditions(9,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]>0 & [Expt.Trials(:).RespDir]~=0;
            conditions(10,:)= [Expt.Trials(:).bd]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip
            conditions(11,:)= [Expt.Trials(:).bd]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip pair to compare 
            if (sum(conditions(7,:)) ==0 )
                conditions(7,:) = [Expt.Trials(:).bd]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip
                conditions(8,:) = [Expt.Trials(:).bd]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip pair to compare 
                conditions(10,:)= [Expt.Trials(:).bd]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip
                conditions(11,:)= [Expt.Trials(:).bd]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0; %flip pair to copmare
            end
        end
    case {'SRID', 'DRID'}
        ORs = unique([Expt.Trials(:).or]);
        %SPc = {Expt.Trials(:).Spikes};
        %for ss = 1: length(SPc)
        %    Spcount(ss) = sum(SPc{ss}<5000);
        %end
        %if  mean(Spcount([Expt.Trials(:).or]==ORs(1)))>mean(Spcount([Expt.Trials(:).or]==ORs(end))) %mean([Expt.Trials([Expt.Trials(:).or]==ORs(1)).count])>mean([Expt.Trials([Expt.Trials(:).or]==ORs(end)).count])
        if mean([Expt.Trials([Expt.Trials(:).or]==ORs(1)).count])>mean([Expt.Trials([Expt.Trials(:).or]==ORs(end)).count])
            Por = ORs(1);
            Nor = ORs(end);
        else
            Nor = ORs(1);
            Por = ORs(end);
        end            
        if(PrefDir == 2)
            conditions(1,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).or]==Por;
            conditions(2,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).or]==Por;
            conditions(3,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).or]==Nor;
            conditions(4,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).or]==Nor;
            conditions(6,:) = 0;
        else
            conditions(1,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).or]==Por;
            conditions(2,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).or]==Por;
            conditions(3,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).or]==Nor;
            conditions(4,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).or]==Nor;
            conditions(6,:) = 0;
        end
    case {'DID', 'DIDB'}
        if(mean([Expt.Trials([Expt.Trials(:).dx]>0).RespDir])>0)
            ResponseToPositive = 1;
            ResponseToNegative = -1;
        else
            ResponseToPositive = -1;
            ResponseToNegative = 1;
        end
        TP = NotDoingVeryGood(Expt);
        if(PrefDir == 2) % DID , ...
            conditions(1,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(2,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(3,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;%7
            conditions(4,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;%8
            %conditions(5,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;%9
            conditions(5,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;%10
            %conditions(6,:)= [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;%11
            conditions(6,:)= [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;%12
            conditions(7,:) = [Expt.Trials(:).dx]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0;
            conditions(8,:)= [Expt.Trials(:).dx]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0;
            % 9, 10 are to be used with 11 ,12  for zscored CP calculation
            conditions(9,:)= [Expt.Trials(:).dx]> min([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]< 0 & TP & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(10,:)=[Expt.Trials(:).dx]> min([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]< 0 & TP & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(11,:)=[Expt.Trials(:).dx]< max([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]> 0 & TP & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(12,:)=[Expt.Trials(:).dx]< max([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]> 0 & TP & [Expt.Trials(:).RespDir]==ResponseToPositive;
            %conditions(11,:)=(conditions(9,:) | conditions(10,:)) & [Expt.Trials(:).RespDir]==ResponseToNegative;
            %conditions(12,:)=(conditions(9,:) | conditions(10,:)) & [Expt.Trials(:).RespDir]==ResponseToPositive;
            % 13 and 14 are the the two bds near zero disparity. 
            conditions(13,:)= [Expt.Trials(:).dx] == max([Expt.Trials([Expt.Trials(:).dx]<0).dx]);
            conditions(14,:)= [Expt.Trials(:).dx] == min([Expt.Trials([Expt.Trials(:).dx]>0).dx]);
            conditions(15,:) = [Expt.Trials(:).RespDir]~=0; % All completed trials for normalization

            % 16 - 19 are the the two bds near zero disparity. 
            conditions(16,:)= [Expt.Trials(:).dx] == max([Expt.Trials([Expt.Trials(:).dx]<0).dx]) & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(17,:)= [Expt.Trials(:).dx] == max([Expt.Trials([Expt.Trials(:).dx]<0).dx]) & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(18,:)= [Expt.Trials(:).dx] == min([Expt.Trials([Expt.Trials(:).dx]>0).dx]) & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(19,:)= [Expt.Trials(:).dx] == min([Expt.Trials([Expt.Trials(:).dx]>0).dx]) & [Expt.Trials(:).RespDir]==ResponseToPositive;

        else 
            conditions(1,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(2,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(3,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;%7
            conditions(4,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;%8
            %conditions(5,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;%9
            conditions(5,:)= [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;%10
            %conditions(7,:)= [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;%11
            conditions(6,:)= [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;%12
            conditions(7,:)= [Expt.Trials(:).dx]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0;
            conditions(8,:)= [Expt.Trials(:).dx]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0;
            % 9, 10 are to be used with 11 ,12  for zscored CP calculation
            conditions(9,:)= [Expt.Trials(:).dx]< max([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]> 0 & TP & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(10,:)=[Expt.Trials(:).dx]< max([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]> 0 & TP & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(11,:)=[Expt.Trials(:).dx]> min([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]< 0 & TP & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(12,:)=[Expt.Trials(:).dx]> min([Expt.Trials(:).dx]) & [Expt.Trials(:).dx]< 0 & TP & [Expt.Trials(:).RespDir]==ResponseToNegative;
            %conditions(11,:)=(conditions(9,:) | conditions(10,:)) & [Expt.Trials(:).RespDir]==ResponseToPositive;
            %conditions(12,:)=(conditions(9,:) | conditions(10,:)) & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(13,:)= [Expt.Trials(:).dx] == min([Expt.Trials([Expt.Trials(:).dx]>0).dx]);
            conditions(14,:)= [Expt.Trials(:).dx] == max([Expt.Trials([Expt.Trials(:).dx]<0).dx]);
            conditions(15,:) = [Expt.Trials(:).RespDir]~=0; % All completed trials for normalization

            % 16 - 19 are the the two bds near zero disparity. 
            conditions(16,:)= [Expt.Trials(:).dx] == min([Expt.Trials([Expt.Trials(:).dx]>0).dx]) & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(17,:)= [Expt.Trials(:).dx] == min([Expt.Trials([Expt.Trials(:).dx]>0).dx]) & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(18,:)= [Expt.Trials(:).dx] == max([Expt.Trials([Expt.Trials(:).dx]<0).dx]) & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(19,:)= [Expt.Trials(:).dx] == max([Expt.Trials([Expt.Trials(:).dx]<0).dx]) & [Expt.Trials(:).RespDir]==ResponseToNegative;
        end
    otherwise
        if(mean([Expt.Trials([Expt.Trials(:).dx]>0).RespDir])>0)
            ResponseToPositive = 1;
            ResponseToNegative = -1;
        else
            ResponseToPositive = -1;
            ResponseToNegative = 1;
        end
        if(PrefDir == 2)
            conditions(1,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(2,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(3,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(4,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(5,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(6,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(7,:) = [Expt.Trials(:).dx]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0;
            conditions(8,:) = [Expt.Trials(:).dx]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0;
        else
            conditions(1,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(2,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(3,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(4,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(5,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(6,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(7,:) = [Expt.Trials(:).dx]== max([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0;
            conditions(8,:) = [Expt.Trials(:).dx]== min([Expt.Trials(:).dx]) & [Expt.Trials(:).RespDir]~=0;
        end
end
if (nargout==3)
    r2p = ResponseToPositive;
    r2n = ResponseToNegative;
end
end