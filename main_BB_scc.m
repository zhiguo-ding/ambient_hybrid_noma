clear  
figure
plalpha_ray = 3; 
plalpha_ric = 2; 
K_ric = 10;
sigma2_dbm= -40;%+10*log10(BW)+Nf; %Thermal noise in dBm
sigma = 10^((sigma2_dbm-30)/10);
R=4; %target data rate
eps0 = exp(R)-1; 
M = 4; %number of   users
 
ct=1000;
 
%eps_BB = 0.001;
max_itr=1000; % for K=4, 1000 is needed
snrdb = [0: 10: 30];
Mvec = [1:5];
 
for icD = 1: length(Mvec) 
    %Ptx = 10^((snrdb(icD)-30)/10);
    M = Mvec(icD);
    for ict = 1 : ct 

        % generate locations of M user
        rP = 2; %size of the square
        rcenter = 20; 
        coodP = [rP*rand(M,1).*sign(randn(M,1)) rP*rand(M,1).*sign(randn(M,1))]; %locations
        coodP = coodP +[ rcenter rcenter]; %the center of the squire is shifted
        rkP = max(1,sqrt(abs(coodP(:,1)).^2+abs(coodP(:,2)).^2));

        %the channel vectors of the   users
        hmfading = abs(complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1))).^2; %[3;2]/10;%
        hmfading = max(hmfading,0.01);% as in the paper of Design of Downlink Hybrid NOMA Transmission, we do clip to avoid singularity
        hm_unoder = hmfading./rkP.^plalpha_ray; %add the path loss
        [hm, ind_order] = sort(hm_unoder,'descend');        
        hm = sort(hm,'descend');
        coodP = coodP(ind_order,:);
        %hm = hm_unoder;
        %the channel vectors between the   users
        for i = 1 : M
            for j = 1 : M
                temp1 = coodP(i,:)-coodP(j,:);
                distij(i,j) = max(1, sqrt(temp1*temp1')); %avoid distance smaller than one
                fadingij = abs(sqrt(K_ric/(K_ric+1))*complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1)) ...
                    +sqrt(1/(K_ric+1))*1 ).^2;
                fadingij = max(fadingij,0.01);% just to avoid signularity 
                gmi(i,j) = fadingij/distij(i,j)^plalpha_ric;% just to avoid signularity 
            end
        end
        %gmi = abs(complex(sqrt(0.5)*randn(M,M),sqrt(0.5)*randn(M,M))).^2; %[1 2 ;3 4]/10;%
        %gmi = max(gmi,0.01);
        gammami = hm.*gmi;  
        gammami = gammami - diag(diag(gammami)); 
        gammami = gammami + diag(hm);   
        gammami = gammami/sigma; %noise nomalization   
        hm = hm/sigma;

        %initilization    
        %tic
        for m = 1 : M
            P_oma(m,1) = (exp(R)-1)/hm(m);              
        end 
       x_ini =  sum(P_oma)*ones(M,1);
        eps_BB = sum(P_oma)/100;

        Bset=[zeros(M,1) x_ini*10];
        if feasibility(Bset(:,2), gammami,M,R)
            Upk =  sum(Bset(:,2));
            Lowk = sum(Bset(:,1));
        else
            Upk = sum(x_ini);Lowk = sum(x_ini);
        end
        Bset = [Bset; Lowk Upk];%odd position is for lower bound
        % the last row of Bset is for low and upper bound values

        bb_k=1;
        low_ind = 1; %index of the set whose lower bound is highest
        while (Upk(end)-Lowk(end)>eps_BB) & (bb_k<max_itr) 
            vec_upper = Bset(end,1:2:end);% lower bound of each set, odd positions %Bset(end,2:2:end);% upper bound of each set, even positions
            [temp,upper_ind] = min(vec_upper);            
            B_temp = Bset(1:end-1,2*(upper_ind-1)+1:2*(upper_ind-1)+2);%find the corresponding set 
            Bset(:,2*(upper_ind-1)+1:2*(upper_ind-1)+2) = [];%remove the set
            length_temp = B_temp(1:end,2)-B_temp(1:end,1);%the last row is not used
            [tempbb, ind_edge] = max(length_temp);%find the longest edge

            Bset_low_half = B_temp;  
            Bset_low_half(ind_edge, 2) = sum(Bset_low_half(ind_edge, :))/2;
            if feasibility(Bset_low_half(:,2), gammami,M,R)  % use f_max to test
                Lowk1 = sum(Bset_low_half(:,1));
                Upk1 = sum(Bset_low_half(:,2));
            else
                Lowk1 = sum(x_ini)*2;
                Upk1 = sum(x_ini)*2;
            end
            Bset_up_half = B_temp;  
            Bset_up_half(ind_edge, 1) = sum(Bset_up_half(ind_edge, :))/2;
            if feasibility(Bset_up_half(:,2), gammami,M,R) 
                Lowk2 = sum(Bset_up_half(:,1));
                Upk2 = sum(Bset_up_half(:,2));
            else
                Lowk2 = sum(x_ini)*2;
                Upk2 = sum(x_ini)*2;
            end        

            Bset = [Bset [Bset_low_half; Lowk1 Upk1]]; %add this new set
            Bset = [Bset [Bset_up_half; Lowk2 Upk2]]; %add this new set
            %Bounding 
            Lowknew = min(Bset(end,1:2:end));%the new lower bound
            Lowk = [Lowk Lowknew]; %  need to remember this index
            Upnew = min(Bset(end,2:2:end));%the new upper bound
            Upk = [Upk Upnew]; % need to remember this index
            %[Lowk; Upk]
            %prunning
            vec_pru = Bset(end,1:2:end);% lower bound of each set, odd positions
            ind_pru = find(vec_pru>Upk(end)); % **** new
            if ~isempty(ind_pru)
                Bset(:,[2*(ind_pru-1)+1 2*(ind_pru-1)+2])=[];
            end
            bb_k = bb_k+1;

            %Bset
            %[Upk;Lowk]
            %sum(P_oma)
        end

        %recover the x
        [tempbb,opt_ind] = min(Bset(end,2:2:end));%find which box contributes the upper bound
        Boptimal = Bset(1:end-1,2*(opt_ind-1)+1:2*(opt_ind-1)+2);%extract the box

        x= Boptimal(:,2);
        etami_final = findeta(x, gammami,M,R); 
        powerallx(ict) = sum(x);
        pict_oma(ict) = sum(P_oma);
        pict_bb(ict) = sum(x);
        %time_ictbb(ict) = toc;
        %[x P_oma; sum(x) sum(P_oma)]
        %etami_final 

                
        %initilization of SCA     
        %tic
        % for m = 1 : M
        %     P_oma(m,1) = (exp(R)-1)/hm(m);            
        % end 
        Pini = zeros(M*(M+1)/2,1);
        for m = 1 : M
            Pini((m-1)*((m-1)+1)/2+1,1) = P_oma(m,1);            
        end 
        Plast = Pini;
        for ica = 1 : 5
            A = []; % No other constraints
            b = [];
            Aeq = [];%[eye(M,M);-eye(M,M)];%-eye(2*Dsize);
            beq = [];%[C;zeros(M,1)];%zeros(2*Dsize,1);
            lb = [];
            ub = [];    
            P0 = Pini;%zeros(length(x_taylor),1); %[x_1 ... x_K real(f_1) image(f_1)... ... f_K]
            options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
            Pnew = fmincon(@(xca) myobjsca(xca,M),P0,A,b,Aeq,beq,lb,ub,@(xca) myconssca(xca, Plast,M,gammami,R) ,options);
            %[Pini Plast Pnew]
            Plast = Pnew;
            
        end
        Puser = zeros(M,1);
        etamsca = zeros(M,M);
        for m = 1 : M
            Puser(m) = Pnew(m*(m+1)/2,1);
            for i = 1 : m
                etamsca(m,i) = Pnew((m-1)*((m-1)+1)/2+i,1)/Pnew(i*(i+1)/2,1);
            end
        end
        pict_sca(ict) = sum(Puser);        
        % %etamsca
        % time_ictsca(ict) = toc;
 
     

    end

    p_bb(icD) = mean(pict_bb)    
    p_oma(icD) = mean(pict_oma)
    p_sca(icD) = mean(pict_sca) 
    %time_sca(icD) = mean(time_ictsca);
    %time_bb(icD) = mean(time_ictbb);
end

plot(Mvec,p_oma,Mvec, p_bb, Mvec,p_sca )

%plot(Mvec,time_sca,Mvec,time_scc)
 
% %check accumulated interference
% [gammami(1,1)*x(1)+gammami(2,1)*etami_final(2,1)*x(1)+gammami(3,1)*etami_final(3,1)*x(1) ...
%     gammami(2,2)*etami_final(2,2)*x(2)+gammami(3,2)*etami_final(3,2)*x(2)...
%     gammami(3,3)*etami_final(3,3)*x(3)
%     ]
  
 
%feasibility check 
function [temp,y] = feasibility(x, gammami,M,R) 
    %x contains M power coefficients, P_1 to P_M
    Pvec = x;
    A = []; % No other constraints
    b = [];
    Aeq = [];%[eye(M,M);-eye(M,M)];%-eye(2*Dsize);
    beq = [];%[C;zeros(M,1)];%zeros(2*Dsize,1);
    lb = [];
    ub = [];    
    % 
    temp = 1;%
    etami = zeros(M,M); %for U_M, this is empty
    for m = M :-1: 2 %starting from U_M, to U_2, there is no reflection for U_1
        %x contains [f1... fK] 
        f0 = zeros(m-1,1);%zeros(length(x_taylor),1); %[x_1 ... x_K real(f_1) image(f_1)... ... f_K]
        options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
        f = fmincon(@(y) 1,f0,A,b,Aeq,beq,lb,ub,@(f) mycons(f,m,gammami,etami,Pvec,M,R) ,options);
        %f = fmincon(@(x) sum(x(1:end-1)),f0,A,b,Aeq,beq,lb,ub,@(x) mycons(x,m,gammami,etami,Pvec,M,R),options);
        %f = fmincon(@(x) sum(min(x(1:end-1))),f0,A,b,Aeq,beq,lb,ub,@(x) mycons(x,m,gammami,etami,Pvec,M,R),options);

        [c,] = mycons(f,m,gammami,etami,Pvec,M,R) ;       
        if max(c)<=0.000001
            temp = 1; %feasible
        else
            temp = 0; %infeasible
        end
        if temp == 0 %once it is not feasible, just quit back to the main BB loop
            break;
        end
        eta_new = [f;1]; % m by 1
        etami(m,:) = [eta_new; zeros(M-m,1)]';
    end

    %for, U_1, we need to check its feasibility
    if temp~=0 %if temp is already 0, no need to discuss U_1
        etami(1,1) = 1;
        temp1 = 0;
        for j = 2 : M
             temp1 = temp1 + sum(gammami(j,1)*etami(j,1)*Pvec(j));
        end
        if log(1+gammami(1,1)*Pvec(1)/(temp1+1))<R
            temp = 0; %infeasible
        end
    end
end 
 

%find reflection coefficients
function [etami] = findeta(x, gammami,M,R) 
    %x contains M power coefficients, P_1 to P_M
    Pvec = x;
    A = []; % No other constraints
    b = [];
    Aeq = [];%[eye(M,M);-eye(M,M)];%-eye(2*Dsize);
    beq = [];%[C;zeros(M,1)];%zeros(2*Dsize,1);
    lb = [];
    ub = [];    
    
    etami = zeros(M,M); %for U_M, this is empty
    for m = M :-1: 2 %starting from U_M, to U_2, there is no reflection for U_1
        %x contains [f1... fK] 
        f0 = zeros(m-1,1);%zeros(length(x_taylor),1); %[x_1 ... x_K real(f_1) image(f_1)... ... f_K]
        options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
        f = fmincon(@(y) 1,f0,A,b,Aeq,beq,lb,ub,@(f) mycons(f,m,gammami,etami,Pvec,M,R) ,options);

        eta_new = [f;1]; % m by 1
        etami(m,:) = [eta_new; zeros(M-m,1)]';
    end

    etami(1,1) = 1;
    
end
 
function [c,ceq] = myconsbb(x, Etalast,Pvec, M,gammami,R) 

%(j-1)*((j-1)+1)/2 is the start of the j-th segment of the vector
%start directly for U_2's eta, not U_1's eta
 
for m =1 : M-1 %these users suffer interference
    sum1 = 0;
    for i =1 : m
        sum2 = 0;
        for j = m : M
            sum2 = sum2 + gammami(j,i)*Pvec(i)*x((j-1)*((j-1)+1)/2+i,1);
        end
        gammamix = zeros(M*(M+1)/2,1);
        for j = m+1 : M
            gammamix((j-1)*((j-1)+1)/2+i,1) = gammami(j,i)*Pvec(i);
        end
        Imi = log(gammamix'*Etalast+1) + gammamix'*(x-Etalast)/(gammamix'*Etalast+1);
        sum1 = sum1 + log(1+sum2) - Imi;
    end

    c(m,1) = R -sum1;    
end

%for user M
    sum1 = 0;
    for i =1 : M
        sum2 = 0;
        for j = M : M
            sum2 = sum2 + gammami(j,i)*Pvec(i)*x((j-1)*((j-1)+1)/2+i,1);
        end 
        sum1 = sum1 + log(1+sum2)  ;
    end

    c(M,1) = R -sum1;    
 

c = [c; x- ones(M*(M+1)/2,1)];%to ensure eta_mi<=1
 
c = [c;-x];%to ensure eta_mi>=0
ceq = [];
end
 

function [c,ceq] = myconsbb_wotsca(x, Etalast,Pvec, M,gammami,R) 

%(j-1)*((j-1)+1)/2 is the start of the j-th segment of the vector
%start directly for U_2's eta, not U_1's eta
 
for m =1 : M-1 %these users suffer interference
    sum1 = 0;
    for i =1 : m
        sum2 = 0;
        for j = m : M
            sum2 = sum2 + gammami(j,i)*Pvec(i)*x((j-1)*((j-1)+1)/2+i,1);
        end
        sum3 = 0;
        for j = m+1 : M
            sum3 = sum3 + gammami(j,i)*Pvec(i)*x((j-1)*((j-1)+1)/2+i,1);
        end 
        sum1 = sum1 + log(1+sum2) - log(1+sum3);
    end

    c(m,1) = R -sum1;    
end

%for user M
    sum1 = 0;
    for i =1 : M
        sum2 = 0;
        for j = M : M
            sum2 = sum2 + gammami(j,i)*Pvec(i)*x((j-1)*((j-1)+1)/2+i,1);
        end 
        sum1 = sum1 + log(1+sum2)  ;
    end

    c(M,1) = R -sum1;    
 

c = [c; x- ones(M*(M+1)/2,1)];%to ensure eta_mi<=1
 
c = [c;-x];%to ensure eta_mi>=0
ceq = [];
end
   
function [objsca] = myobjsca(xca,M)

objsca = 0;
for m = 1 : M
    objsca = objsca + xca(m*(m+1)/2);
end
end

 

function [c,ceq] = myconssca(x, Plast,M,gammami,R) 
%(j-1)*((j-1)+1)/2 is the start of the j-th segment of the vector
for m =1 : M-1 %these users suffer interference
    sum1 = 0;
    for i =1 : m
        sum2 = 0;
        for j = m : M
            sum2 = sum2 + gammami(j,i)*x((j-1)*((j-1)+1)/2+i,1);
        end
        hmi = zeros(M*(M+1)/2,1);
        for j = m+1 : M
            hmi((j-1)*((j-1)+1)/2+i,1) = gammami(j,i);
        end
        Imi = log(hmi'*Plast+1) + hmi'*(x-Plast)/(hmi'*Plast+1);
        sum1 = sum1 + log(1+sum2) - Imi;
    end

    c(m,1) = R -sum1;    
end

%for user M
    sum1 = 0;
    for i =1 : M
        sum2 = 0;
        for j = M : M
            sum2 = sum2 + gammami(j,i)*x((j-1)*((j-1)+1)/2+i,1);
        end 
        sum1 = sum1 + log(1+sum2)  ;
    end

    c(M,1) = R -sum1;    
 

Pmmvec = zeros(M*(M+1)/2,1); 
for m = 1 : M
    for i = 1 : m
        Pmmvec((m-1)*((m-1)+1)/2+i,1) = x(i*(i+1)/2,1);
    end
end
c = [c; x- Pmmvec];%to ensure Pmi<=Pii
 
c = [c;-x];%to ensure Pmi>=0
ceq = [];
end
 


function [c,ceq] = myconssccM(x, M, R,gammami)  
%(j-1)*((j-1)+1)/2 is the start of the j-th segment of the vector
 
sum1 = 0;
for i =1 : M
    sum1 = sum1 + log(1+gammami(M,i)*x(i));
end
c(1,1) = R -sum1;    
c = [c;-x];
ceq = [];
end


function [c,ceq] = myconssccother(x, n, M, R,gammami,delta_all,Pmi_all)  
%(j-1)*((j-1)+1)/2 is the start of the j-th segment of the vector
deltan = x(1:n); %the first n elements are delta
Pni = x(n+1:n+n);% the rest n for the U_n's reflecting coefficients 
sum1 = 0;
for i =1 : n
    sum2 = 0;
    for j = n+1: M
        sum2 = sum2 + gammami(j,i)*Pmi_all(j,i);
    end
    sum1 = sum1 + log(1+gammami(n,i)*Pni(i)/(1+sum2));
end
c(1,1) = R -sum1;    
c = [c;-x]; %nonnegative
c = [c; Pni-Pmi_all(M,1:n)' - deltan - sum(delta_all(1:n,n:M-1),2)];%each row of delta_all is for one Pmm
c = [c; Pmi_all(n+1:M,n) - Pni(n)];
ceq = [];
end

function [c,ceq] = mycons(f,m,gammami,etami,Pvec,M,R) 
 %f is eta
 f = [f;1]; % the last eta is one
 
 sum1=0;
for i = 1 : m     
    temp1 = 0;
    for j = m+1 : M
         temp1 = temp1 + sum(gammami(j,i)*etami(j,i)*Pvec(j));
    end
    sum1 = sum1 + log(1+gammami(m,i)*f(i)*Pvec(i)/(temp1+1));     
end 
c(1,1) = R-sum1;  
 
c = [c; f-ones(m,1)];
c = [c; -f];
ceq = [];
end
 