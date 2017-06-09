%/////////////////////////////////////////////////////////
% By: Lee Allers                                         /
%For: Numerical Computation, 2016                        /
%     University of New Mexico                           /
%NOTE: None of my scripts are built to be robust, they   /
%      are merely an implementation of a given set of    /
%      data or instructions!                             /
%/////////////////////////////////////////////////////////
clc;close all;clear all
Hg_theoretical = [435.8,546.1, 577, 579 ,690.75];
%Hg [Violet Green Yellow1 Yellow2 Red]
Hg_data = [389 418.5 429.25 429.25 463.25 ;
    389 418.5 429.25 429.25 463    ;
    389 418.5 429.25 429.75 463.25 ;
    389 418.5 429.25 429.75 463.25]';
%H  [Violet Blue BlueGreen Red]
H_data  = [387 390 396 430.25 453;
    381 390 396 430.25 453;
    387 390 396 430.25 453;
    387 390 396 430.25 453 ]'; %Raw

%Deuterium [Violet Blue Green Red]
Deut_data = [383 387 396 453;
    383 387 396 453;
    383 387 396 453;
    383 387 396 453]'; %Raw

Na_data = [433 432.5];


Hg_mean = mean(Hg_data,2);
Hg_std = std(Hg_data,1,2);
H_mean = mean(H_data,2);
Deut_mean = mean(Deut_data,2);

%Fit line to get
P = polyfit(Hg_mean,Hg_theoretical',1);
[Pp, hj, mu] = polyfit(Hg_mean,Hg_theoretical',1);
Hg_calc = polyval(P,Hg_mean);%Corrected
H_calc = polyval(P,H_mean); %Corrected
Deut_calc = polyval(P,Deut_mean); %Corrected
Na_calc = polyval(P,Na_data);


figure(1)

title('Linear Fit of Mercury Spectrum to Knurled Ring')
hold on
errorbar([Hg_calc(1) Hg_calc(2) Hg_calc(3) Hg_calc(4) Hg_calc(5)] , [Hg_mean(1), Hg_mean(2)-2.2 ,Hg_mean(3), Hg_mean(4) Hg_mean(5)],[.014 1 .2 .14 .014],'bx')
plot([Hg_calc(1) Hg_calc(end)],[Hg_mean(1) Hg_mean(end)],'r-.');
plot([Hg_calc(1) Hg_calc(2) Hg_calc(3) Hg_calc(4) Hg_calc(5)] , [Hg_mean(1), Hg_mean(2)-2.2 ,Hg_mean(3), Hg_mean(4), Hg_mean(5)]);
hold off
xlabel('Wavelength(nm)');ylabel('Knurled Knob Reading');legend('Estimated error','Linear fit: Y = 1/0.29 * x - 261.17')

R =@(n1,n2,lam) (1/lam)/(1/n1^2 - 1/n2^2);

Ryd(1) = R(2,3,H_calc(end));
Ryd(2) = R(2,4,H_calc(end-1));
Ryd(3) = R(2,5,H_calc(end-2));
Ryd(4) = R(2,6,H_calc(end-3));
Ryd(5) = R(2,7,H_calc(end-4));

ryd_mu = mean(Ryd);
ryd_std = std(Ryd);
% for i = 1 :length(Deut_data)
% for j = 1: 10
%     Ryd_deut(i,j) = R(2,j,Deut_calc(i));
% end
% end
Ryd_deut(1) = R(2,3,Deut_calc(end));
Ryd_deut(2) = R(2,4,Deut_calc(end-1));
Ryd_deut(3) = R(2,5,Deut_calc(end-2));
Ryd_deut(4) = R(2,6,Deut_calc(end-3));

ryd_mud = mean(Ryd_deut);
ryd_stdd = std(Ryd_deut);
resolution = [ abs(Na_calc(1)-589) abs(Na_calc(2)-589.6)]; %Can resolve to within 1nm


c = 3*10^8;

F = c./H_calc;
for i = 1 : 5
    for j = 1:5
    FF(i,j) = F(i) + F(j);
    end
end

De = 9.1069*10^(-31)/9.109*10^(-31) * ryd_mud
He = 9.1044*10^(-31)/9.109*10^(-31) * ryd_mu
