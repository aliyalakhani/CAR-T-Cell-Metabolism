%% Load data
load Compat.mat

%% Normalize flux values
x_norm = []; %CAR T cell panel
y_norm = []; %Cancer cell panel
x_norm2 = []; %CAR T cell panel includes nucleobases for overflow metabolism score
x_sem = []; %CAR T cell panel
y_sem = []; %Cancer cell panel
x_sem2 = []; %CAR T cell panel includes nucleobases for overflow metabolism score
compatibility_raji = []; 
compatibility_ramos = [];
raji_sem = [];
ramos_sem = [];

%CAR T cell panel
for j = size(Ratesx,1):-1:1
    for m = size(Ratesx,2):-1:1
        if Ratesx(j,:) > 0
            x_norm(j,m) = (abs(Ratesx(j,m)) - min(abs(Ratesx(j,:)))) / (max(abs(Ratesx(j,:))) - min(abs(Ratesx(j,:))));
        elseif Ratesx(j,:) < 0
            x_norm(j,m) = (abs(Ratesx(j,m)) - min(abs(Ratesx(j,:)))) / (max(abs(Ratesx(j,:))) - min(abs(Ratesx(j,:))));
        else 
            x_norm(j,m) = Ratesx(j,m)/max(abs(Ratesx(j,:)));
        end
    end
end

%Cancer cell panel
for j = size(Ratesy,1):-1:1
    for m = size(Ratesy,2):-1:1
        if Ratesy(j,:) > 0
            y_norm(j,m) = (abs(Ratesy(j,m)) - min(abs(Ratesy(j,:)))) / (max(abs(Ratesy(j,:))) - min(abs(Ratesy(j,:))));
        elseif Ratesy(j,:) < 0
            y_norm(j,m) = (abs(Ratesy(j,m)) - min(abs(Ratesy(j,:)))) / (max(abs(Ratesy(j,:))) - min(abs(Ratesy(j,:))));
        else 
            y_norm(j,m) = Ratesy(j,m)/max(abs(Ratesy(j,:)));
        end
    end
end

%CAR T cell panel with nucleobases for overflow metabolism score
for j = size(Ratesx2,1):-1:1
    for m = size(Ratesx2,2):-1:1
        if Ratesx2(j,:) > 0
            x_norm2(j,m) = (abs(Ratesx2(j,m)) - min(abs(Ratesx2(j,:)))) / (max(abs(Ratesx2(j,:))) - min(abs(Ratesx2(j,:))));
        elseif Ratesx2(j,:) < 0
            x_norm2(j,m) = (abs(Ratesx2(j,m)) - min(abs(Ratesx2(j,:)))) / (max(abs(Ratesx2(j,:))) - min(abs(Ratesx2(j,:))));
        else 
            x_norm2(j,m) = Ratesx2(j,m)/max(abs(Ratesx2(j,:)));
        end
    end
end

%Propagated Error for x_norm
for j = size(Ratesx,1):-1:1
    for m = size(Ratesx,2):-1:1
        if Ratesx(j,:) > 0
            minRatesx = min(abs(Ratesx(j,:)));
            maxRatesx = max(abs(Ratesx(j,:)));

            minSEMxIndex = find(abs(Ratesx) == min(abs(Ratesx(j,:))));
            maxSEMxIndex = find(abs(Ratesx) == max(abs(Ratesx(j,:))));

            term1 = (SEMx(j,m)) / (maxRatesx - minRatesx);
            term2 = (SEMx(maxSEMxIndex).*(minRatesx - abs(Ratesx(j,m)))) / (maxRatesx - minRatesx)^2;
            term3 = (SEMx(minSEMxIndex).*(abs(Ratesx(j,m))- maxRatesx)) / (maxRatesx - minRatesx)^2;

            x_sem(j,m) = sqrt(term1^2 + term2^2 + term3^2); 
        elseif Ratesx(j,:) < 0
            minRatesx = min(abs(Ratesx(j,:)));
            maxRatesx = max(abs(Ratesx(j,:)));

            minSEMxIndex = find(abs(Ratesx) == min(abs(Ratesx(j,:))));
            maxSEMxIndex = find(abs(Ratesx) == max(abs(Ratesx(j,:))));

            term1 = (SEMx(j,m)) / (maxRatesx - minRatesx);
            term2 = (SEMx(maxSEMxIndex).*(minRatesx - abs(Ratesx(j,m)))) / (maxRatesx - minRatesx)^2;
            term3 = (SEMx(minSEMxIndex).*(abs(Ratesx(j,m))- maxRatesx)) / (maxRatesx - minRatesx)^2;

            x_sem(j,m) = sqrt(term1^2 + term2^2 + term3^2);       
        else 
            x_sem(j,m) = sqrt((SEMx(j,m)./Ratesx(j,m))^2 + (SEMx(maxSEMxIndex)./maxRatesx)^2);
        end
    end
end

%Propagated Error for y_norm
for j = size(Ratesy,1):-1:1
    for m = size(Ratesy,2):-1:1
        if Ratesy(j,:) > 0
            minRatesy = min(abs(Ratesy(j,:)));
            maxRatesy = max(abs(Ratesy(j,:)));

            minSEMyIndex = find(abs(Ratesy) == min(abs(Ratesy(j,:))));
            maxSEMyIndex = find(abs(Ratesy) == max(abs(Ratesy(j,:))));

            term1 = (SEMy(j,m)) / (maxRatesy - minRatesy);
            term2 = (SEMy(maxSEMyIndex).*(minRatesy - abs(Ratesy(j,m)))) / (maxRatesy - minRatesy)^2;
            term3 = (SEMy(minSEMyIndex).*(abs(Ratesy(j,m))- maxRatesy)) / (maxRatesy - minRatesy)^2;

            y_sem(j,m) = sqrt(term1^2 + term2^2 + term3^2); 
        elseif Ratesy(j,:) < 0
            minRatesy = min(abs(Ratesy(j,:)));
            maxRatesy = max(abs(Ratesy(j,:)));

            minSEMyIndex = find(abs(Ratesy) == min(abs(Ratesy(j,:))));
            maxSEMyIndex = find(abs(Ratesy) == max(abs(Ratesy(j,:))));

            term1 = (SEMy(j,m)) / (maxRatesy - minRatesy);
            term2 = (SEMy(maxSEMyIndex).*(minRatesy - abs(Ratesy(j,m)))) / (maxRatesy - minRatesy)^2;
            term3 = (SEMy(minSEMyIndex).*(abs(Ratesy(j,m))- maxRatesy)) / (maxRatesy - minRatesy)^2;

            y_sem(j,m) = sqrt(term1^2 + term2^2 + term3^2);        
        else 
            y_sem(j,m) = sqrt((SEMy(j,m)./Ratesy(j,m))^2 + (SEMy(maxSEMyIndex)./maxRatesy)^2);
        end
    end
end

%Propagated Error for x_norm2 for overflow metabolism score
for j = size(Ratesx2,1):-1:1
    for m = size(Ratesx2,2):-1:1
        if Ratesx2(j,:) > 0
            minRatesx = min(abs(Ratesx2(j,:)));
            maxRatesx = max(abs(Ratesx2(j,:)));

            minSEMxIndex = find(abs(Ratesx2) == min(abs(Ratesx2(j,:))));
            maxSEMxIndex = find(abs(Ratesx2) == max(abs(Ratesx2(j,:))));

            term1 = (SEMx2(j,m)) / (maxRatesx - minRatesx);
            term2 = (SEMx2(maxSEMxIndex).*(minRatesx - abs(Ratesx2(j,m)))) / (maxRatesx - minRatesx)^2;
            term3 = (SEMx2(minSEMxIndex).*(abs(Ratesx2(j,m))- maxRatesx)) / (maxRatesx - minRatesx)^2;

            x_sem2(j,m) = sqrt(term1^2 + term2^2 + term3^2); 
        elseif Ratesx2(j,:) < 0
            minRatesx = min(abs(Ratesx2(j,:)));
            maxRatesx = max(abs(Ratesx2(j,:)));

            minSEMxIndex = find(abs(Ratesx2) == min(abs(Ratesx2(j,:))));
            maxSEMxIndex = find(abs(Ratesx2) == max(abs(Ratesx2(j,:))));

            term1 = (SEMx2(j,m)) / (maxRatesx - minRatesx);
            term2 = (SEMx2(maxSEMxIndex).*(minRatesx - abs(Ratesx2(j,m)))) / (maxRatesx - minRatesx)^2;
            term3 = (SEMx2(minSEMxIndex).*(abs(Ratesx2(j,m))- maxRatesx)) / (maxRatesx - minRatesx)^2;

            x_sem2(j,m) = sqrt(term1^2 + term2^2 + term3^2);       
        else 
            x_sem2(j,m) = sqrt((SEMx2(j,m)./Ratesx2(j,m))^2 + (SEMx2(maxSEMxIndex)./maxRatesx)^2);
        end
    end
end

%% Compatibility

%Raji Compatibility
for j = size(Ratesx,1):-1:1
    for m = size(Ratesx,2):-1:1
        if x_norm(j,m) >= 0 && y_norm(j,1) >= 0
            compatibility_raji(j,m) = abs(x_norm(j,m)+y_norm(j,1))-x_norm(j,m).*y_norm(j,1).*(1+abs(x_norm(j,m)+y_norm(j,1))-x_norm(j,m).*y_norm(j,1));
        elseif x_norm(j,m) < 0 && y_norm(j,1) < 0
            compatibility_raji(j,m) = abs(x_norm(j,m)+y_norm(j,1))-x_norm(j,m).*y_norm(j,1).*(1+abs(x_norm(j,m)+y_norm(j,1))-x_norm(j,m).*y_norm(j,1));
        else 
            compatibility_raji(j,m) = min(2, -x_norm(j,m).*y_norm(j,1)/abs(x_norm(j,m)+y_norm(j,1))+abs(x_norm(j,m)-y_norm(j,1)));
        end
    end
end

%Ramos Compatibility
for j = size(Ratesx,1):-1:1
    for m = size(Ratesx,2):-1:1
        if x_norm(j,m) >= 0 && y_norm(j,2) >= 0
            compatibility_ramos(j,m) = abs(x_norm(j,m)+y_norm(j,2))-x_norm(j,m).*y_norm(j,2).*(1+abs(x_norm(j,m)+y_norm(j,2))-x_norm(j,m).*y_norm(j,2));
        elseif x_norm(j,m) < 0 && y_norm(j,2) < 0
            compatibility_ramos(j,m) = abs(x_norm(j,m)+y_norm(j,2))-x_norm(j,m).*y_norm(j,2).*(1+abs(x_norm(j,m)+y_norm(j,2))-x_norm(j,m).*y_norm(j,2));
        else 
            compatibility_ramos(j,m) = min(2, -x_norm(j,m).*y_norm(j,2)/abs(x_norm(j,m)+y_norm(j,2))+abs(x_norm(j,m)-y_norm(j,2)));
        end
    end
end


%Propagated Error for Raji Compatibility
for j = size(Ratesx,1):-1:1
    for m = size(Ratesx,2):-1:1
        if x_norm(j,m) >= 0 && y_norm(j,1) >= 0
            
            a = ((-y_norm(j,1)^2).*(x_norm(j,m)+y_norm(j,1)))/(abs(x_norm(j,m)+y_norm(j,1)))^3;
            b = ((-x_norm(j,m)^2).*(x_norm(j,m)+y_norm(j,1)))/(abs(x_norm(j,m)+y_norm(j,1)))^3;
            c = (x_norm(j,m) - y_norm(j,1))/abs(x_norm(j,m) - y_norm(j,1));
            d = (y_norm(j,1) - x_norm(j,m))/abs(x_norm(j,m) - y_norm(j,1));

            raji_sem(j,m) = sqrt(a^2.*x_sem(j,m)^2 + b^2.*y_sem(j,1)^2 + c^2.*x_sem(j,m)^2 + d^2.*y_sem(j,1)^2);  
        elseif x_norm(j,m) < 0 && y_norm(j,1) < 0


            a = ((-y_norm(j,1)^2).*(x_norm(j,m)+y_norm(j,1)))/(abs(x_norm(j,m)+y_norm(j,1)))^3;
            b = ((-x_norm(j,m)^2).*(x_norm(j,m)+y_norm(j,1)))/(abs(x_norm(j,m)+y_norm(j,1)))^3;
            c = (x_norm(j,m) - y_norm(j,1))/abs(x_norm(j,m) - y_norm(j,1));
            d = (y_norm(j,1) - x_norm(j,m))/abs(x_norm(j,m) - y_norm(j,1));

            raji_sem(j,m) = sqrt(a^2.*x_sem(j,m)^2 + b^2.*y_sem(j,1)^2 + c^2.*x_sem(j,m)^2 + d^2.*y_sem(j,1)^2); 
        else
            a = (x_norm(j,m) + y_norm(j,1))/abs(x_norm(j,m) + y_norm(j,1));
            b = (x_norm(j,m) + y_norm(j,1))/abs(x_norm(j,m) + y_norm(j,1));
            c = y_norm(j,1)*(-abs(x_norm(j,m) + y_norm(j,1)) - ((x_norm(j,m).*(x_norm(j,m) + y_norm(j,1)))/abs((x_norm(j,m) + y_norm(j,1)))) + 2.*x_norm(j,m).*y_norm(j,1) -1);
            d = x_norm(j,m)*(-abs(x_norm(j,m) + y_norm(j,1)) - ((y_norm(j,1).*(x_norm(j,m) + y_norm(j,1)))/abs((x_norm(j,m) + y_norm(j,1)))) + 2.*x_norm(j,m).*y_norm(j,1) -1);
            
            raji_sem(j,m) = sqrt(a^2.*x_sem(j,m)^2 + b^2.*y_sem(j,1)^2 + c^2.*x_sem(j,m)^2 + d^2.*y_sem(j,1)^2);
        end
    end
end

%Propagated Error for Ramos Compatibility
for j = size(Ratesx,1):-1:1
    for m = size(Ratesx,2):-1:1
        if x_norm(j,m) >= 0 && y_norm(j,2) >= 0
            
            a = ((-y_norm(j,2)^2).*(x_norm(j,m)+y_norm(j,2)))/(abs(x_norm(j,m)+y_norm(j,2)))^3;
            b = ((-x_norm(j,m)^2).*(x_norm(j,m)+y_norm(j,2)))/(abs(x_norm(j,m)+y_norm(j,2)))^3;
            c = (x_norm(j,m) - y_norm(j,2))/abs(x_norm(j,m) - y_norm(j,2));
            d = (y_norm(j,2) - x_norm(j,m))/abs(x_norm(j,m) - y_norm(j,2));

            ramos_sem(j,m) = sqrt(a^2.*x_sem(j,m)^2 + b^2.*y_sem(j,2)^2 + c^2.*x_sem(j,m)^2 + d^2.*y_sem(j,2)^2);  
        elseif x_norm(j,m) < 0 && y_norm(j,2) < 0


            a = ((-y_norm(j,2)^2).*(x_norm(j,m)+y_norm(j,2)))/(abs(x_norm(j,m)+y_norm(j,2)))^3;
            b = ((-x_norm(j,m)^2).*(x_norm(j,m)+y_norm(j,2)))/(abs(x_norm(j,m)+y_norm(j,2)))^3;
            c = (x_norm(j,m) - y_norm(j,2))/abs(x_norm(j,m) - y_norm(j,2));
            d = (y_norm(j,2) - x_norm(j,m))/abs(x_norm(j,m) - y_norm(j,2));

            ramos_sem(j,m) = sqrt(a^2.*x_sem(j,m)^2 + b^2.*y_sem(j,2)^2 + c^2.*x_sem(j,m)^2 + d^2.*y_sem(j,2)^2); 
        else
            a = (x_norm(j,m) + y_norm(j,2))/abs(x_norm(j,m) + y_norm(j,2));
            b = (x_norm(j,m) + y_norm(j,2))/abs(x_norm(j,m) + y_norm(j,2));
            c = y_norm(j,2)*(-abs(x_norm(j,m) + y_norm(j,2)) - ((x_norm(j,m).*(x_norm(j,m) + y_norm(j,2)))/abs((x_norm(j,m) + y_norm(j,2)))) + 2.*x_norm(j,m).*y_norm(j,2) -1);
            d = x_norm(j,m)*(-abs(x_norm(j,m) + y_norm(j,2)) - ((y_norm(j,2).*(x_norm(j,m) + y_norm(j,2)))/abs((x_norm(j,m) + y_norm(j,2)))) + 2.*x_norm(j,m).*y_norm(j,2) -1);
            
            ramos_sem(j,m) = sqrt(a^2.*x_sem(j,m)^2 + b^2.*y_sem(j,2)^2 + c^2.*x_sem(j,m)^2 + d^2.*y_sem(j,2)^2);
        end
    end
end

            
            