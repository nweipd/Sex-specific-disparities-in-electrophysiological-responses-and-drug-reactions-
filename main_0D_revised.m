close all
clear all
load Xfinal_total.mat
%% I changed index for ko and fkatp, and make them be one index. 
constants_0D
CL=300; % ms 
beats=1000;
beats_recorded=4;
options=[];%options for ode solver
dt=0.1; % output time step 

flag_ode=1;
flag_drug_vector=1:14;
n_drug=length(flag_drug_vector); % choose 4 dose Bepridil, Diltiazem, Nifedipine 2,


case_total=zeros(4,5); % each row represent a trace; each column represent a parameter 
case_total(:,1)=0; % endo 
case_total(:,2)=[0;2;0;4]; % hormone 
case_total(:,3)=[1;1;2;2];% sex 
case_total(:,4)=0.2; %fkatp 
case_total(:,5)=9; % ko 




APD_total_drug=zeros(2,size(case_total,1),n_drug); 
Xfinal_total_drug=zeros(41,size(case_total,1),n_drug);
X_total_drug=zeros(beats_recorded*CL/dt+1,41,size(case_total,1),n_drug); % output for X to plot the traces for currents and voltages. 

% this for loop only record X for the new beat and replace. 

num_drug=1
for flag_drug=flag_drug_vector
    num=1
    for row=1:size(case_total,1)
        flag_cell=case_total(row,1);
        flag_hormone=case_total(row,2);
        flag_gender=case_total(row,3);
        fkatp=case_total(row,4);
        ko=case_total(row,5);
        X0=Xfinal_total(:,num);
        iter=0
        APDlast2=100;
        APDevenlast2=100;
        APDoddlast2=100;
        flag_alternans=-1;
        while flag_alternans==-1 && iter<=iter_max
            iter=iter+1
            APD_test=zeros(beats_recorded,1); % save last 4 beasts of each iter to compare 
            X_output=zeros(beats_recorded*CL/dt+1,41); % save X of the last 4 beats of each iter
            num_beat=1;
            for n=[1:beats]
                [time X]=ode15s(@model_0D,[0:dt:CL],X0,options,flag_gender,flag_cell,flag_ode,flag_hormone,ko,fkatp,flag_drug);
                X0=X(end,:);
                if n>=beats-3
                    v=X(:,1);
                    APD=compute_APD_revised(time,v);
                    APD_test(num_beat)=APD;
                    X_output((num_beat-1)*CL/dt+1:num_beat*CL/dt+1,:)=X;
                    num_beat=num_beat+1;
                end
            end
            APDlast2=abs(APD_test(end-1)-APD_test(end));
            APDevenlast2=abs(APD_test(end-2)-APD_test(end));
            APDoddlast2=abs(APD_test(end-3)-APD_test(end-1));
            if APDlast2<alternans_thr
                flag_alternans=0;
            elseif APDevenlast2<alternans_thr && APDoddlast2<alternans_thr
                flag_alternans=1;
            else
                flag_alternans=-1;
            end
        end
        if flag_alternans==0 || flag_alternans==1
            APD_total_drug(:,num,num_drug)=APD_test(end-1:end); % use num as index, not n. 
            Xfinal_total_drug(:,num,num_drug)=X(end,:); 
            if flag_alternans==0
                warning('1:1 behavior')
            elseif flag_alternans==1
                warning('alternans')
            end
        elseif flag_alternans==-1
            warning('max iter is reached, neither')
        end
        X_total_drug(:,:,num,num_drug)=X_output;
        num=num+1
    end
    num_drug=num_drug+1  
end
   



save Xfinal_total_drug.mat Xfinal_total_drug
save X_total_drug.mat X_total_drug 
save APD_total_drug.mat APD_total_drug




                




