%% Load data and process
load product_conc.mat
load Mf.mat
load constant.mat

%% Fmincon fitting
x = zeros(1,1);
score = zeros(1,1);
R1 = []; %Secretion or uptake rate for donor 1
R2 = []; %Secretion or uptake rate for donor 2
R3 = []; %Secretion or uptake rate for donor 3
R4 = []; %Secretion or uptake rate for donor 4
S1 = ones(size(P1D1))*1e9; %Score for donor 1
S2 = ones(size(P1D2))*1e9; %Score for donor 2
S3 = ones(size(P1D3))*1e9; %Score for donor 3
S4 = ones(size(P1D4))*1e9; %Score for donor 4
Fitted1 = []; %Simulated nM values for donor 1
Fitted2 = []; %Simulated nM values for donor 2
Fitted3 = []; %Simulated nM values for donor 3
Fitted4 = []; %Simulated nM values for donor 4
r1 = []; %Secretion or uptake rate for donor 1 over first 24 hours
r2 = []; %Secretion or uptake rate for donor 2 over first 24 hours
r3 = []; %Secretion or uptake rate for donor 3 over first 24 hours
r4 = []; %Secretion or uptake rate for donor 4 over first 24 hours

%Donor1
for j = size(P1D1,2):-1:1 %samples
    for m = size(P1D1,1):-1:1 %all amino acids, glucose, lactate, pyruvate, hydroxyproline, myo-inositol, ammonia
        for i = 1:150             
            constantsv = constantD1(j,:);
            Mfv = Mf(m,1); %nM in fresh media RPMI
            Observed = [P2D1(m,j); P3D1(m,j)]; %experimental nmoles
            
            P1v = P1D1(m,j);
            Rus0 = -100+200*rand(1,1);
            A = Mf(m,5); %constraints: R<1 for essential amino acids and glucose, R>1 for alanine, glutamate, lactate, pyruvate, ammonia
            b = 0;
            Aeq = [];
            beq = [];

            F = @(Rus)Rsim(Rus,Observed,constantsv,Mfv,P1v);

            [x, fval] = fmincon(F,Rus0,A,b,Aeq,beq);
 
            if fval < S1(m,j)                      
                R1(m,j) = x;
                S1(m,j) = fval;
                [~, ~,Fitted1{m,j},r1(m,j)] = F(x);
            end

            [~, ~, ~, ~, P0(m,j),P1(m,j),P1s(m,j),P2(m,j),P2s(m,j),P3(m,j)] = Rsim(x, Observed, constantsv, Mfv, P1v);
        end
    end  
end

%Donor2
for j = size(P1D2,2):-1:1 %samples
    for m = size(P1D2,1):-1:1 %all amino acids, glucose, lactate, pyruvate, hydroxyproline, myo-inositol, ammonia
        for i = 1:150
            constantsv = constantD2(j,:);
            Mfv = Mf(m,2); %nM in fresh media RPMI
            Observed = P3D2(m,j); %experimental nmoles
            
            P1v = P1D2(m,j);
            Rus0 = -100+200*rand(1,1);
            A = Mf(m,5); %constraints: R<1 for essential amino acids and glucose, R>1 for alanine, glutamate, lactate, pyruvate, ammonia
            b = 0;
            Aeq = [];
            beq = [];

            F = @(Rus)Rsim(Rus,Observed,constantsv,Mfv,P1v);

            [x, fval] = fmincon(F,Rus0,A,b,Aeq,beq);

            if fval < S2(m,j)                      
                R2(m,j) = x;
                S2(m,j) = fval;
                [~, ~, Fitted2{m,j},r2(m,j)] = F(x);
            end
            
            [~, ~, ~, ~, P0(m,j),P1(m,j),P1s(m,j),P2(m,j),P2s(m,j),P3(m,j)] = Rsim(x, Observed, constantsv, Mfv, P1v);
        end
    end  
end

%Donor3
for j = size(P1D3,2):-1:1 %samples
    for m = size(P1D3,1):-1:1 %all amino acids, glucose, lactate, pyruvate, hydroxyproline, myo-inositol
        for i = 1:150
            constantsv = constantD3(j,:);
            Mfv = Mf(m,3); %mM in fresh media RPMI
            Observed = P2D3(m,j); %experimental nmoles
            P1v = P1D3(m,j);
            Rus0 = -100+200*rand(1,1);
            A = Mf(m,5); %constraints: R<1 for essential amino acids and glucose, R>1 for alanine, glutamate, lactate, pyruvate, ammonia
            b = 0;
            Aeq = [];
            beq = [];

            F = @(Rus)Rsim(Rus,Observed,constantsv,Mfv,P1v);

            [x, fval] = fmincon(F,Rus0,A,b,Aeq,beq);

            if fval < S3(m,j)                      
                R3(m,j) = x;
                S3(m,j) = fval;
                [~, ~,Fitted3{m,j},r3(m,j)] = F(x);
            end
            
            [~, ~, ~, ~, P0(m,j),P1(m,j),P1s(m,j),P2(m,j),P2s(m,j),P3(m,j)] = Rsim(x, Observed, constantsv, Mfv, P1v);
        end
    end  
end

%Donor4
for j = size(P1D4,2):-1:1 %samples
    for m = size(P1D4,1):-1:1 %all amino acids, glucose, lactate, pyruvate, hydroxyproline, myo-inositol, ammonia
        for i = 1:150
            constantsv = constantD4(j,:);
            Mfv = Mf(m,4); %nM in fresh media RPMI
            Observed = [P2D4(m,j); P3D4(m,j)]; %experimental nmoles
            P1v = P1D4(m,j);
            Rus0 = -100+200*rand(1,1);
            A = Mf(m,5); %constraints: R<1 for essential amino acids and glucose, R>1 for alanine, glutamate, lactate, pyruvate, ammonia
            b = 0;
            Aeq = [];
            beq = [];

            F = @(Rus)Rsim(Rus,Observed,constantsv,Mfv,P1v);

            [x, fval] = fmincon(F,Rus0,A,b,Aeq,beq);

            if fval < S4(m,j)                      
                R4(m,j) = x;
                S4(m,j) = fval;
                [~, ~, Fitted4{m,j},r4(m,j)] = F(x);
  
            end
                      
            [~, ~, ~, ~, P0(m,j),P1(m,j),P1s(m,j),P2(m,j),P2s(m,j),P3(m,j)] = Rsim(x, Observed, constantsv, Mfv, P1v);
        end
    end  
end
                             
%% Score function
function [score,vscore,Fitted,r,P0,P1,P1s,P2,P2s,P3] = Rsim(Rus, Observed_in, constantsv_in, Mfv_in, P1v_in)
        
    k  = 0.0416666670000000;
    C0 = constantsv_in(2);
    C1 = constantsv_in(3);
    C2 = constantsv_in(4);
    C3 = constantsv_in(5);
    V0 = constantsv_in(6);
    V1 = constantsv_in(7);
    V2 = constantsv_in(8);
    Vr_1 = constantsv_in(9);
    Vr_2 = constantsv_in(10);
    Va_1 = constantsv_in(11);
    Va_2 = constantsv_in(12);
    Mf = Mfv_in(1);
    P0 = Mf*V0;    
    R = Rus;
       
    
    
    if  size(Observed_in,1) == 1 && Vr_2 > 0 && P1v_in(1) ~= 0
        P1 = P1v_in(1);
        r = (P1-P0)*k/((C0+C1)/2);
        P1s = P1 - (Vr_1/V1)*P1 + Va_1*Mf; 
        P2 = P1s + R*((C2+C1)/2)/k;
        P2s = P2 - (Vr_2/V2)*P2 + Va_2*Mf;
        P3 = P2s + R*((C3+C2)/2)/k;
        Fitted = P3;
    elseif  size(Observed_in,1) == 1 && Vr_2 > 0 && P1v_in(1) == 0
        P1 = P0 + R*((C0+C1)/2)/k;
        r = (P1-P0)*k/((C0+C1)/2);
        P1s = P1 - (Vr_1/V1)*P1 + Va_1*Mf; 
        P2 = P1s + R*((C2+C1)/2)/k;
        P2s = P2 - (Vr_2/V2)*P2 + Va_2*Mf;
        P3 = P2s + R*((C3+C2)/2)/k;
        Fitted = P3;
    elseif size(Observed_in,1) == 1 && Vr_2 == 0
        P1 = P1v_in(1);
        r = (P1-P0)*k/((C0+C1)/2);
        P1s = P1 - (Vr_1/V1)*P1 + Va_1*Mf; 
        P2 = P1s + R*((C2+C1)/2)/k;
        P2s = 0;
        P3 = 0;
        Fitted = P2;
    elseif size(Observed_in,1) == 2
        P1 = P1v_in(1);
        r = (P1-P0)*k/((C0+C1)/2);
        P1s = P1 - (Vr_1/V1)*P1 + Va_1*Mf; 
        P2 = P1s + R*((C2+C1)/2)/k;
        P2s = P2 - (Vr_2/V2)*P2 + Va_2*Mf;
        P3 = P2s + R*((C3+C2)/2)/k;
        Fitted = [                                                                                                                                                                                                                                     
        P2                                                                                      
        P3
        ];
    end
    
    score = (Observed_in - Fitted).^2;
    score = score./Observed_in.^2;
    vscore = 0;
    score = sum(score);
end