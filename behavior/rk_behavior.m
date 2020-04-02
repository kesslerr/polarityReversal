%% process self-written Psychopy logfile(s) of CHF experiment


%% analysis folder
analysis_folder='/media/cth/Samsung USB/CHF/processing_new';
cd(analysis_folder);
subcounter=1;

%% loop over subjects
for nsub=4:10
    if nsub<10,     nsubs=strcat('0',num2str(nsub));
    else,           nsubs=num2str(nsub);
    end

    cd(strcat('RK_200',nsubs,'/log'));
    
    %% import files
    load('log_act_color.txt');
    log_event=textread('log_event.txt','%s');
    log_hemifield=textread('log_hemifield.txt','%s');
    load('log_stimfreq.txt');
    load('log_time_block.txt');
    load('log_time_exp.txt');

    %% hitrates

    % get indices of colorchange 1
    %ind = find(log_act_color==1);

    % partition the colorchanges

    %% hitrates
    % First colorchange
    color=0; counter=1; hitrate=[]; RT=[];
    for i=1:length(log_event)
        if log_act_color(i) == 1
            if color==0;
                color=1;
                color_start_ind = i;
            elseif color==1
            end
        elseif log_act_color(i) == 0
            if color == 0;
            elseif color==1
                color = 0;
                color_stop_ind = i;
                % look in the time window for responses
                for j=color_start_ind:color_stop_ind
                    if strcmp(log_event{j},'keypress_1_correct')
                        RT(counter)=log_time_exp(j)-log_time_exp(color_start_ind);
                        hitrate(counter)=1;     counter=counter+1; break,
                    elseif strcmp(log_event{j},'keypress_1_incorrect')
                        RT(counter)=log_time_exp(j)-log_time_exp(color_start_ind);
                        hitrate(counter)=0;     counter=counter+1; break,
                    end
                    % if no hit at the last index, it is a miss
                    if j ==color_stop_ind
                        hitrate(counter) = 0;
                        RT(counter)=log_time_exp(j)-log_time_exp(color_start_ind);
                        counter=counter+1; break,
                    end
                end
            end
        end  
    end
    hitrate1=mean(hitrate);
    % hitrates>1 are caused by pauses or so, so delete them
    RT1=mean(RT(find(RT<1)));

    % Second colorchange
    color=0; counter=1; hitrate=[]; RT=[];
    for i=1:length(log_event)
        if log_act_color(i) == 2
            if color==0;
                color=1;
                color_start_ind = i;
            elseif color==1
            end
        elseif log_act_color(i) == 0
            if color == 0;
            elseif color==1
                color = 0;
                color_stop_ind = i;
                % look in the time window for responses
                for j=color_start_ind:color_stop_ind
                    if strcmp(log_event{j},'keypress_2_correct')
                        RT(counter)=log_time_exp(j)-log_time_exp(color_start_ind);
                        hitrate(counter)=1;     counter=counter+1; break,
                    elseif strcmp(log_event{j},'keypress_2_incorrect')
                        RT(counter)=log_time_exp(j)-log_time_exp(color_start_ind);
                        hitrate(counter)=0;     counter=counter+1; break,
                    end
                    % if no hit at the last index, it is a miss
                    if j ==color_stop_ind
                        hitrate(counter) = 0;
                        RT(counter)=log_time_exp(j)-log_time_exp(color_start_ind);
                        counter=counter+1; break,
                    end
                end
            end
        end  
    end
    
    %% verrechnen
    hitrate2=mean(hitrate);
    RT2=mean(RT(find(RT<1)));

    hitrate=(hitrate1+hitrate2)/2;
    RT=(RT1+RT2)/2;
    
    TOTALhitrate(subcounter)=hitrate;
    TOTALRT(subcounter)     =RT;
    TOTALhitrate1(subcounter)=hitrate1;
    TOTALRT1(subcounter)     =RT1;
    TOTALhitrate2(subcounter)=hitrate2;
    TOTALRT2(subcounter)     =RT2;
    
    subcounter=subcounter+1;
    
    cd(analysis_folder);
    
end

TOTALhitrate=TOTALhitrate'; TOTALRT=TOTALRT';
TOTALhitrate1=TOTALhitrate1'; TOTALRT1=TOTALRT1';
TOTALhitrate2=TOTALhitrate2'; TOTALRT2=TOTALRT2';

clearvars -except TOTALhitrate TOTALRT TOTALhitrate1 TOTALRT1 TOTALhitrate2 TOTALRT2
