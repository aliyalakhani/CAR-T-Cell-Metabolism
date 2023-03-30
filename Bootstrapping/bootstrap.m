load bootstrap_data.mat

%% Bootstrapping Code for citrate

output_cit = zeros(4,8);

for i = 2:5 %CD19(2), Leu16(3), Rituximab(4), RFR-LCDR(5)
    for m = 100000:-1:1

        x = floor(12*rand(1,12))+1; %Pick 12 random positions from pooled data
        y = floor(6*rand(1,6))+1; %Pick 6 random positions from pooled data
        
        EGFRt_glc = cit_glc(1,:);
        CAR_glc = cit_glc(i,:);
        EGFRt_CAR_glc = [EGFRt_glc CAR_glc];

        EGFRt_gln = cit_gln(1,:);
        CAR_gln = cit_gln(i,:);
        EGFRt_CAR_gln = [EGFRt_gln CAR_gln];
    
        %Mean of glucose measurements for EGFRt and CAR
        meanEGFRt_CAR_glc_sample1 = mean(EGFRt_CAR_glc(x(1:6)));
        meanEGFRt_CAR_glc_sample2 = mean(EGFRt_CAR_glc(x(7:12)));
        
        %Mean of glutamine measurements for EGFRt and CAR
        meanEGFRt_CAR_gln_sample1 = mean(EGFRt_CAR_gln(y(1:3)));
        meanEGFRt_CAR_gln_sample2 = mean(EGFRt_CAR_gln(y(4:6)));
        
        %Get ratio of means and iterate 100,000 times
        ratioEGFRt_CAR1(m) = meanEGFRt_CAR_gln_sample1/meanEGFRt_CAR_glc_sample1;
        ratioEGFRt_CAR2(m) = meanEGFRt_CAR_gln_sample2/meanEGFRt_CAR_glc_sample2;

        %Difference between two ratios
        differenceEGFRt_CAR(m) = ratioEGFRt_CAR1(m) - ratioEGFRt_CAR2(m);
    end

    %Determine p-values
    figure();
    histogram(differenceEGFRt_CAR);
    output_cit(i-1,1) = (sum(differenceEGFRt_CAR >= abs(citrate_obs(i-1))) + sum(differenceEGFRt_CAR <= -1*abs(citrate_obs(i-1))))/size(differenceEGFRt_CAR,2);
    output_cit(i-1,2) = (sum(differenceEGFRt_CAR <= abs(citrate_obs(i-1))) + sum(differenceEGFRt_CAR >= -1*abs(citrate_obs(i-1))))/size(differenceEGFRt_CAR,2);
    
    output_cit(i-1,3) = mean(differenceEGFRt_CAR);
    output_cit(i-1,4) = std(differenceEGFRt_CAR);
    output_cit(i-1,5) = mean(ratioEGFRt_CAR1);
    output_cit(i-1,6) = std(ratioEGFRt_CAR1);
    output_cit(i-1,7) = mean(ratioEGFRt_CAR2);
    output_cit(i-1,8) = std(ratioEGFRt_CAR2);
end

%% Bootstrapping Code for aKG

output_aKG = zeros(4,8);

for i = 2:5 %CD19(2), Leu16(3), Rituximab(4), RFR-LCDR(5)
    for m = 100000:-1:1

        x = floor(12*rand(1,12))+1; %Pick 12 random positions from pooled data
        y = floor(6*rand(1,6))+1; %Pick 6 random positions from pooled data
        
        EGFRt_glc = aKG_glc(1,:);
        CAR_glc = aKG_glc(i,:);
        EGFRt_CAR_glc = [EGFRt_glc CAR_glc];

        EGFRt_gln = aKG_gln(1,:);
        CAR_gln = aKG_gln(i,:);
        EGFRt_CAR_gln = [EGFRt_gln CAR_gln];
    
        %Mean of glucose measurements for EGFRt and CAR
        meanEGFRt_CAR_glc_sample1 = mean(EGFRt_CAR_glc(x(1:6)));
        meanEGFRt_CAR_glc_sample2 = mean(EGFRt_CAR_glc(x(7:12)));
        
        %Mean of glutamine measurements for EGFRt and CAR
        meanEGFRt_CAR_gln_sample1 = mean(EGFRt_CAR_gln(y(1:3)));
        meanEGFRt_CAR_gln_sample2 = mean(EGFRt_CAR_gln(y(4:6)));
        
        %Get ratio of means and iterate 100,000 times
        ratioEGFRt_CAR1(m) = meanEGFRt_CAR_gln_sample1/meanEGFRt_CAR_glc_sample1;
        ratioEGFRt_CAR2(m) = meanEGFRt_CAR_gln_sample2/meanEGFRt_CAR_glc_sample2;

        %Difference between two ratios
        differenceEGFRt_CAR(m) = ratioEGFRt_CAR1(m) - ratioEGFRt_CAR2(m);
    end

    %Determine p-values
    figure();
    histogram(differenceEGFRt_CAR);
    output_aKG(i-1,1) = (sum(differenceEGFRt_CAR >= abs(aKG_obs(i-1))) + sum(differenceEGFRt_CAR <= -1*abs(aKG_obs(i-1))))/size(differenceEGFRt_CAR,2);
    output_aKG(i-1,2) = (sum(differenceEGFRt_CAR <= abs(aKG_obs(i-1))) + sum(differenceEGFRt_CAR >= -1*abs(aKG_obs(i-1))))/size(differenceEGFRt_CAR,2);
    
    output_aKG(i-1,3) = mean(differenceEGFRt_CAR);
    output_aKG(i-1,4) = std(differenceEGFRt_CAR);
    output_aKG(i-1,5) = mean(ratioEGFRt_CAR1);
    output_aKG(i-1,6) = std(ratioEGFRt_CAR1);
    output_aKG(i-1,7) = mean(ratioEGFRt_CAR2);
    output_aKG(i-1,8) = std(ratioEGFRt_CAR2);
end

%% Bootstrapping Code for fumarate

output_fum = zeros(4,8);

for i = 2:5 %CD19(2), Leu16(3), Rituximab(4), RFR-LCDR(5)
    for m = 100000:-1:1

        x = floor(12*rand(1,12))+1; %Pick 12 random positions from pooled data
        y = floor(6*rand(1,6))+1; %Pick 6 random positions from pooled data
        
        EGFRt_glc = fum_glc(1,:);
        CAR_glc = fum_glc(i,:);
        EGFRt_CAR_glc = [EGFRt_glc CAR_glc];

        EGFRt_gln = fum_gln(1,:);
        CAR_gln = fum_gln(i,:);
        EGFRt_CAR_gln = [EGFRt_gln CAR_gln];
    
        %Mean of glucose measurements for EGFRt and CAR
        meanEGFRt_CAR_glc_sample1 = mean(EGFRt_CAR_glc(x(1:6)));
        meanEGFRt_CAR_glc_sample2 = mean(EGFRt_CAR_glc(x(7:12)));
        
        %Mean of glutamine measurements for EGFRt and CAR
        meanEGFRt_CAR_gln_sample1 = mean(EGFRt_CAR_gln(y(1:3)));
        meanEGFRt_CAR_gln_sample2 = mean(EGFRt_CAR_gln(y(4:6)));
        
        %Get ratio of means and iterate 100,000 times
        ratioEGFRt_CAR1(m) = meanEGFRt_CAR_gln_sample1/meanEGFRt_CAR_glc_sample1;
        ratioEGFRt_CAR2(m) = meanEGFRt_CAR_gln_sample2/meanEGFRt_CAR_glc_sample2;

        %Difference between two ratios
        differenceEGFRt_CAR(m) = ratioEGFRt_CAR1(m) - ratioEGFRt_CAR2(m);
    end

    %Determine p-values
    figure();
    histogram(differenceEGFRt_CAR);
    output_fum(i-1,1) = (sum(differenceEGFRt_CAR >= abs(fumarate_obs(i-1))) + sum(differenceEGFRt_CAR <= -1*abs(fumarate_obs(i-1))))/size(differenceEGFRt_CAR,2);
    output_fum(i-1,2) = (sum(differenceEGFRt_CAR <= abs(fumarate_obs(i-1))) + sum(differenceEGFRt_CAR >= -1*abs(fumarate_obs(i-1))))/size(differenceEGFRt_CAR,2);
    
    output_fum(i-1,3) = mean(differenceEGFRt_CAR);
    output_fum(i-1,4) = std(differenceEGFRt_CAR);
    output_fum(i-1,5) = mean(ratioEGFRt_CAR1);
    output_fum(i-1,6) = std(ratioEGFRt_CAR1);
    output_fum(i-1,7) = mean(ratioEGFRt_CAR2);
    output_fum(i-1,8) = std(ratioEGFRt_CAR2);
end

%% Bootstrapping Code for malate

output_mal = zeros(4,8);

for i = 2:5 %CD19(2), Leu16(3), Rituximab(4), RFR-LCDR(5)
    for m = 100000:-1:1

        x = floor(12*rand(1,12))+1; %Pick 12 random positions from pooled data
        y = floor(6*rand(1,6))+1; %Pick 6 random positions from pooled data
        
        EGFRt_glc = mal_glc(1,:);
        CAR_glc = mal_glc(i,:);
        EGFRt_CAR_glc = [EGFRt_glc CAR_glc];

        EGFRt_gln =mal_gln(1,:);
        CAR_gln = mal_gln(i,:);
        EGFRt_CAR_gln = [EGFRt_gln CAR_gln];
    
        %Mean of glucose measurements for EGFRt and CAR
        meanEGFRt_CAR_glc_sample1 = mean(EGFRt_CAR_glc(x(1:6)));
        meanEGFRt_CAR_glc_sample2 = mean(EGFRt_CAR_glc(x(7:12)));
        
        %Mean of glutamine measurements for EGFRt and CAR
        meanEGFRt_CAR_gln_sample1 = mean(EGFRt_CAR_gln(y(1:3)));
        meanEGFRt_CAR_gln_sample2 = mean(EGFRt_CAR_gln(y(4:6)));
        
        %Get ratio of means and iterate 100,000 times
        ratioEGFRt_CAR1(m) = meanEGFRt_CAR_gln_sample1/meanEGFRt_CAR_glc_sample1;
        ratioEGFRt_CAR2(m) = meanEGFRt_CAR_gln_sample2/meanEGFRt_CAR_glc_sample2;

        %Difference between two ratios
        differenceEGFRt_CAR(m) = ratioEGFRt_CAR1(m) - ratioEGFRt_CAR2(m);
    end

    %Determine p-values
    figure();
    histogram(differenceEGFRt_CAR);
    output_mal(i-1,1) = (sum(differenceEGFRt_CAR >= abs(malate_obs(i-1))) + sum(differenceEGFRt_CAR <= -1*abs(malate_obs(i-1))))/size(differenceEGFRt_CAR,2);
    output_mal(i-1,2) = (sum(differenceEGFRt_CAR <= abs(malate_obs(i-1))) + sum(differenceEGFRt_CAR >= -1*abs(malate_obs(i-1))))/size(differenceEGFRt_CAR,2);
    
    output_mal(i-1,3) = mean(differenceEGFRt_CAR);
    output_mal(i-1,4) = std(differenceEGFRt_CAR);
    output_mal(i-1,5) = mean(ratioEGFRt_CAR1);
    output_mal(i-1,6) = std(ratioEGFRt_CAR1);
    output_mal(i-1,7) = mean(ratioEGFRt_CAR2);
    output_mal(i-1,8) = std(ratioEGFRt_CAR2);
end