% function [tau_ini,tau_conf,z_ini,z_conf,prob_ini,prob_conf]=calculate_kendall(x_ini,y_ini,x_con,y_con)

% //      
% //      denominator = sqrt(concordant_pairs + discordant_pairs+ extra-y) * sqrt(concordant_pairs + discordant_pairs+ extra-x)
% //              denominator = sqrt(LeftNum) * sqrt(RightNum)
% //      where concordant pairs are those differences (xi - xj) that have the same sign as the differences (yi - yj)
% //      extra-y are those non-zero differences (yi - yj) where (xi - xj) = 0
% //      extra-x are those non-zero differences (xi - xj) where (yi - yj) = 0
% // 
load spikesanddistancesforTAUS_40msINTLIB.mat;
spike_alldistINI_INT=spike_alldistINI;
spike_alldistCONF_INT=spike_alldistCONF;
spike_smalldistINI_INT=spike_smalldistINI;
spike_smalldistCONF_INT=spike_smalldistCONF;
animal =1
i_trial =3
                x_ini=spike_smalldistINI_PlC{animal,i_trial,x_cell}(:,2);
                y_ini=spike_smalldistINI_INT{animal,i_trial,y_cell}(:,2);
                x_con=spike_smalldistCONF_PlC{animal,i_trial,x_cell}(:,2);
                y_con=spike_smalldistCONF_INT{animal,i_trial,y_cell}(:,2);

x_cell=1;
y_cell=1
           xys = {x_ini(:,1);y_ini(:,1);x_con(:,1);y_con(:,1)};
            for trtype=1:2
                x=xys{1+(2*(trtype-1)),:};
                y=xys{2+(2*(trtype-1)),:};
                
            concordant_pairs=0;
            discordant_pairs=0;
            extray = 0;
            extrax = 0;
            N=size(x,1);
            Var = (4*N + 10) / (9*N*(N - 1));
            for i= 1:N-1
                for j = i+1:N

                    a=x(i) - x(j);
                    b=y(i) - y(j);
                    if a==0 && b==0
                       continue
                    end
                    
                   c = a*b;
                     if  c>0
                         concordant_pairs = concordant_pairs + 1;
                     elseif c<0
                         discordant_pairs = discordant_pairs + 1;
                     else % c==0
                         
                         if a==0
                             extray=extray+1;
                         else 
                             extrax=extrax+1;
                         end
                     end
                        
                         
                    end
            end
%             concordant_pairs
%             discordant_pairs
%             extrax
%             extray
            numerator = (concordant_pairs - discordant_pairs);
            denominator = sqrt(concordant_pairs + discordant_pairs+ extray) * sqrt(concordant_pairs + discordant_pairs+ extrax);
            if trtype==1
            tau_ini = numerator / denominator;
            z_ini = tau_ini/ sqrt(Var);
            tempz =abs(z_ini)/sqrt(2);
            t=1.0/(1.0+0.5*tempz);
            prob_ini=t*exp(-tempz*tempz-1.26551223+t*(1.00002368+t*(0.37409196+t*...
                (0.09678418+ t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+ t*(-0.82215223+t*0.17087277)))))))));
            else % trtype==2
            tau_conf = numerator / denominator;   
            z_conf = tau_conf / sqrt(Var);      
            tempz =abs(z_conf)/sqrt(2);
            t=1.0/(1.0+0.5*tempz);
            prob_conf=t*exp(-tempz*tempz-1.26551223+t*(1.00002368+t*(0.37409196+t*...
                (0.09678418+ t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+ t*(-0.82215223+t*0.17087277)))))))));            
            end
    end
          
% end
            
        
    