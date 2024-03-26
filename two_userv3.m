clear 
M=2;
R=1;
K_ric = 10;
sigma2_dbm= -40;%+10*log10(BW)+Nf; %Thermal noise in dBm
sigma = 10^((sigma2_dbm-30)/10);
plalpha_ray = 3; 
plalpha_ric = 2; 
ct = 5000;
case_result = [];
Rvec = [3 : 0.5 : 5];
 
for icD =  1: length(Rvec)
    R = Rvec(icD);
    eps0 = exp(R)-1;
    for ict = 1:ct
        % generate locations of M user
        rP = 5; %size of the square
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
        

        gamma0 = gammami(2,1);% 
        gamma1 = gammami(1,1); gamma2 = gammami(2,2);

        %test pure noma
        %gamma2=gamma1/eps0/(1+eps0)-0.1; gamma0 = gamma1/(1+eps0)-0.1; %+ means second case, - means first case
        %gamma2 = gamma0/(1+eps0)-0.1; gamma1 = gamma0*(1+eps0)+0.1;
        %[gamma1 gamma2 gamma0/gamma2]

        A = []; % No other constraints
        b = [];
        Aeq = [];%
        beq = [];%zeros(M+1,1); 
        lb = [];
        ub = [];
        x0 =  ones(3,1); %initilization 
        options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
        x = fmincon(@(x) sum(x(1:2)),x0,A,b,Aeq,beq,lb,ub,@(x) mycons(x,M,R,gamma1,gamma2, gamma0,eps0),options);

        %%case 1 lambda3=0
        temp1 = sqrt(eps0*exp(R)/gamma2/gamma1);
        x1 = [temp1; temp1-1/gamma2; gamma1/eps0/gamma0*temp1-1/gamma0];
        con_x1 = [real(log(1+gamma0*x1(3))+log(1+gamma2*x1(2))-R); -(eps0*gamma0*x1(3)+eps0-gamma1*x1(1));x1(1)-x1(3)];

        %case 2 lambda3~=0 lambda2=0
        temp2 = sqrt(exp(R)/gamma0/gamma2);
        x2 = [temp2-1/gamma0; temp2-1/gamma2; temp2-1/gamma0];
        con_x2 = [real(log(1+gamma0*x2(3))+log(1+gamma2*x2(2))-R); -(eps0*gamma0*x2(3)+eps0-gamma1*x2(1));x2(1)-x2(3)];

        %case 3 lambda3~=0 lambda2~=0
        x3 = [eps0/(gamma1-eps0*gamma0); exp(R)/(eps0*gamma0*gamma2/(gamma1-eps0*gamma0)+gamma2)-1/gamma2 ; eps0/(gamma1-eps0*gamma0)];
        lam1 = exp(R)*(gamma1-eps0*gamma0)/gamma1/gamma2;
        lam2 = gamma0*exp(R)*(eps0*gamma0-gamma1)/gamma1^2/gamma2 - 1/(eps0*gamma0-gamma1);
        lam3 = 1-gamma1*lam2;
        con_x3 = [x3 ; lam2 ;lam3;real(log(1+gamma0*x3(3))+log(1+gamma2*x3(2))-R); -(eps0*gamma0*x3(3)+eps0-gamma1*x3(1));x3(1)-x3(3)];
        % OMA
        x4 = [ eps0/gamma1 ; eps0/gamma2;0];
        %pure NOMA
        pure1 = [eps0/gamma0 ;0; eps0/gamma0];
        con_pure1 = [1/(1+eps0)-gamma2/gamma0; gamma1/gamma0-(1+eps0)];
        pure2 = [eps0*(1+eps0)/gamma1;0; eps0/gamma0];
        con_pure2 = [1/(1+eps0)/eps0-gamma2/gamma1; (1+eps0)-gamma1/gamma0];

        con_all = [min(con_pure1) min(con_pure2) min(con_x1) min(con_x2) min(con_x3)];

        all = [x pure1 pure2 x1 x2 x3  ]; 
        [z1,z2] = max(con_all);
        ana_ict(ict) = sum(all(1:2,z2+1));
        sim_ict(ict) = sum(x(1:2));  
        oma_ict(ict) = sum(x4);
        
        case_ict(icD,ict) = z2;

        

    
    end
    
    for m =1 : 5
        temp1 = find(case_ict(icD,:)==m);
        case_icD(icD,m) = length(temp1)/ct;
    end
    ana(icD) = mean(ana_ict)
    sim(icD) = mean(sim_ict)
    oma(icD) = mean(oma_ict)
end

plot(Rvec, oma, Rvec, sim, Rvec, ana)

%plot(Rvec, 10*log10(oma/10)+30, Rvec,  10*log10(sim/10) +30, Rvec,10*log10(ana/10) +30 )

%%%%%%%%%%%%%%%%%%
function [c,ceq] = mycons(x,M,R,gamma1,gamma2, gamma0,eps0)

 

c(1,1) = R - log(1+gamma0*x(3)) - log(1+gamma2*x(2));
c(2,1) = eps0*gamma0*x(3)+eps0-gamma1*x(1);%R - log(1+gamma1*x(1)/(1+gamma0*x(3)));
c(3,1) = x(3)-x(1);
c(4:6,1) = -x;
ceq = [];
end