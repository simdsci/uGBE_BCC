%%
% Copyright (c) 2018, King Mongkut's University of Technology Thonburi,
% KMUTT. Producted at KMUTT. Written by Ooraphan Chirayutthanasak,
% ooraphan.chira@kmutt.ac.th.
% All rights reserved. This file is port of uGBE. For details, see
% Supplementary data section of the article by Ooraphan Chirayutthanasak, 
% Taira Okita, Somsak Dangtip, Gregory S. Rohrer, Sutatch Ratanaphan,
% and Rajchawit Sarochawikasit 
% % "Universal function for grain boundary energies in bcc metals" 

% uGBE program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% WGBE program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% Copyright (c) 2013, Lawrence Livermore National Security, 
% LLC. Produced at the Lawrence Livermore National Laboratory. 
% Written by Vasily Bulatov, bulatov1@llnl.gov; Bryan Reed, 
% reed12@llnl.gov. CODE-LLNL-CODE-646028. 
% All rights reserved. This file is part of GB5DOF. For details, 
% see Supplementary data section of the original article by V. V. Bulatov, 
% B. W. Reed and M. Kumar "Grain Boundary Energy Function for FCC Metals".

%%
function en = uGBE(fnData, GBE_coh_twin, type)
% fnData is input data file that contain 18 columns of the orientations 
    % in the two grains (P and Q).
% GBE_twin is coherent twin energy in unit of J/m2.
% type of alkali or transition-BCC metals ('transition', 'alkali').

% The grain boundary energy function that approximate the energy from 5D
% space bounaries in BCC Fe metal.  
%
% There are 4 sections in this code. First section, main section contains 
% the code of full energy approximation that based on the study of Bulatov.
% Second section, there are the functions for calculate energy in 1D, 2D,
% and 3D spaces for <100>, <110>, and <111> sets. 
% Third section contains the code for calculating rsw 
% and rsw piecewise functions that based on the study of Dette. 
% The last section is the function of scaffolding calculation 
% that duplicate from the study of Bulatov. 
% V. V. Bulatov, B. W. Reed, M. Kumar, Acta Mater. 65 (2014) 161-175.
% H. Dette, J. Goesmann, C. Greiff, R. Janisch, Acta Mater. 125 (2017) 145-153

    % The only one input of this function is the excel file, which define 
    % the orientations of the two grains (P and Q) by a laboratory frame.
    % There are 18 columns for 3x3 matrix of P and Q cube frame. 
    %
    % Example of the excel input for coherent twin boundary in Fe is as follows:
    %
    % |                  Matrix P                  |                  Matrix Q                  |
    % | x1 |    |    | y1 |    |    | z1 |    |    | x2 |    |    | y2 |    |    | z2 |    |    |
    % |  4 |  2 |  2 |  1 | -1 | -1 |  0 |  2 | -2 |  4 |  2 |  2 | -1 |  1 |  1 |  0 | -2 |  2 |
    % 
    % coherent twin energy = 0.262260282 J/m2, and type = 'transition'
    % This particular excel file pass into UGBE(). Function returns
    %% 0.2204 %% as a energy value of this coherent twin boundary in unit of J/m2.
    %%
    labFrame = readmatrix(fnData);
    [m,~]=size(labFrame);
    
    % There are 50 parameters for BCC fitting that compose of 6 weighting
    % parameters for approximate the energy in 5D space, 32 1D-parameters,
    % 6 2D-parameters, and 5 3D-parameters.
    

    params = [0.775, 0.6, 0.3875, 0.79, 0.95, 0.41, 0.78711, 0.71504, 0.8459, 0.78418, 0.81895, ...
        0.55547, 0.79707, 1.0557, 0.96797, 0.78984, 0.91816, 0.4166, 0.84414, 0.19648, 0.89844, ...
        0.70918, 0.86309, 0.59258, 1.9086, 2.7182, 0.49309, 0.44688, 0.475, 0.23867, 0.47617, ...
        0.37133, 0.97564, 0.45762, 0.8252, 0.20234, 0.5334, 0.82012, 0.24023, 0.53379, 0.91406, ...
        1.0012, 3.1816, -0.092578, -0.066992, 0.69434, -0.21406, 0.41602, 0.74668, 0.45312];
    
    for i=1:m
        
        % Normalize the cube frame of input file
        [P,Q]=convert2PQ(labFrame(i,:));
        
        % Generate geometry parameters
        geom100 = distances_to_set(P,Q,'100'); 
        geom110 = distances_to_set(P,Q,'110');
        geom111 = distances_to_set(P,Q,'111');
    
        % Calculated the grain boundary energy for all P,Q frames, the 5D
        % space boundaries are approximated by the summation of weight
        % normalized grain boundary energy of possible scaffolding in 3 sets.
        en(i) = weightallset(geom100, geom110, geom111, params);
 
    end
    
% find Ergb that is nearly the highest energy at sigma 23 {310} boundary 
% Ergb is estimated by coherent twin energy
    if type == "transition"
        ergb = (4.519*GBE_coh_twin) + 0.2232;
    elseif type == "alkali"
        ergb = (5.8306*GBE_coh_twin) + 0.0092;
    end

% normalized grain boundary energy is 
       approx_en = approx_from_rgb(en, GBE_coh_twin, ergb);
       
%     % The input file needs the column 19, which is the simulation energy
%     % list. If the 19th column is not available, the graph will not present. 
%     if n==19 
%         printGraph(en', labFrame(:,19));
%     end
%     
    fnFitEn = 'AppoxEn.csv'; %',fnData,'
    csvwrite(fnFitEn, approx_en');

end
%%
function approx_en = approx_from_rgb(en, coh_twin, ergb)
a = 0.2;
b= 1;
    approx_en(en>=0.2) = ((en(en>=0.2)-a)*(ergb - coh_twin)/(b-a)) + coh_twin;
    approx_en(en<0.2) = en(en<0.2)*coh_twin/a;
end

function en = weightallset(geom100, geom110, geom111, params)
%

    %eRGB = params(1); % The only dimensioned parameter.  The energy scale, set by the energy of a random boundary.
                        % However, this parameter are significate use in the
                        % universal function for BCC metels in future work.
    d0100 = params(1); % Maximum distance for the 100 set.  Also the distance scale for the rsw weighting function.
    d0110 = params(2); % Same for the 110 set
    d0111 = params(3); % Same for the 111 set
    weight100 = params(4); % Weight for the 100 set, relative to the random boundary
    weight110 = params(5); % Same for 110
    weight111 = params(6); % Same for 111

    % The following three energy lists are in units of eRGB. 
    ksiData = geom100(2,:);
    etaData = geom100(3,:);
    phiData = geom100(4,:);
    e100    = enMix100(ksiData, etaData, phiData, params);
    
    ksiData = geom110(2,:);
    etaData = geom110(3,:);
    phiData = geom110(4,:);
    e110    = enMix110(ksiData, etaData, phiData, params);
    
    ksiData = geom111(2,:);
    etaData = geom111(3,:);
    phiData = geom111(4,:);
    e111    = enMix111(ksiData, etaData, phiData, params);

    d100 = geom100(1,:);
    d110 = geom110(1,:);
    d111 = geom111(1,:);
    
    offset = 0.00001;  % Cutoff of weighting function at small d, purely for numerical purposes
    
    % Now calculate the weights, in a manner designed to give an rsw-like
    % function of d.  Note it calculates a weight for every representation of
    % the boundary within each set.
    s100    = sin(pi/2*d100/d0100);
    s100(d100>d0100) = 1 ;  % Weight saturates at zero above d0
    s100(d100<offset*d0100) = offset*pi/2;  % Avoid calculation of NaN's, replace with something small but finite
    w100    = (1./(s100.*(1-0.5*log(s100)))-1)*weight100;

    s110    = sin(pi/2*d110/d0110);
    s110(d110>d0110) = 1;
    s110(d110<offset*d0110) = offset*pi/2;
    w110    = (1./(s110.*(1-0.5*log(s110)))-1)*weight110;

    s111    = sin(pi/2*d111/d0111);
    s111(d111>d0111) = 1;
    s111(d111<offset*d0111) = offset*pi/2;
    w111    = (1./(s111.*(1-0.5*log(s111)))-1)*weight111;

    en = (sum(e100.*w100)+sum(e110.*w110)+sum(e111.*w111)+1)/(sum(w100)+sum(w110)+sum(w111)+1);
    
end

%
% calculate boudnary energy for <100> set
%

function en = enMix100(ksi, eta, phi, params)
%

    entwist = twist100(ksi, params); 
    entilt = atgb100(ksi, eta, params);   
    
    en = Fmix100(phi, entwist, entilt, params);

end

function en = Fmix100(phi, en1, en2, params)
%

    a = params(46); % 100 mix shape factor

    % select only twist energy larger than or equal to  tilt energy
    select = en1>=en2;
    
    en = zeros(size(phi));
    
    en(select) = en2(select) + (en1(select)-en2(select)).*rsw(phi(select),pi/2,0,a);
    en(~select) = en1(~select) + (en2(~select)-en1(~select)).*rsw(phi(~select),0,pi/2,a);

end

function en = atgb100(ksi, eta, params)
%
    
    period = pi/2 ;

    en1 = stgb100(ksi, params);         % energy of symmetirc tilt energy
    en2 = stgb100(period-ksi, params);  % energy of invert phase of symmetirc tilt energy
    
    eta(eta > period) = 2*period - eta(eta>period);
    
    a = params(43); % 100 ATGB shape factor
    
    en = Fatgb(eta, en1, en2, period, a);

end

function en = stgb100(ksi, params)
%

    % Configuring the number of change point as k.
    % Following the number of chenge point, the Ksi and Energy at change 
    % point are defined by vtheta and gamma, respectively. 
    %
    % This technique will be used in all function of 1D and 2D that call
    % Ftwgb_stgb() function
    k=6;
    vtheta=[0 params(12) 36.85*pi/180 params(13) 53*pi/180 params(14) pi/2];
    gamma=[0 params(7) params(8) params(9) params(10) params(11) 0];
%36.8 
    
    a = 0.5;
    
    en = Ftwgb_stgb(ksi, vtheta, gamma, k, a);

end

function en = twist100(ksi, params)
%

    ksi(ksi > pi/4) = pi/2 - ksi(ksi > pi/4);


    % The same with stgb100() function
    k = 3;
    vtheta = [0 params(18) 36.8*pi/180 pi/4];
    gamma = [0 params(15) params(16) params(17)];
    
    a = 0.5; 
    
    en = Ftwgb_stgb(ksi, vtheta, gamma, k, a);
    
end

%
% calculate boudnary energy for <110> set
%

function en = enMix110(ksi, eta, phi, params)
%

    entwist = twist110(ksi, params);
    entilt = atgb110(ksi, eta, params);
    
    en = Fmix110(ksi, eta, phi, entwist, entilt, params);

end

function en = Fmix110(ksi, eta, phi, entwist, entilt, params)
%

    en = zeros(size(phi));
    
    en1 = entwist;
    en2 = entilt;
    
    select = en1>=en2;
    
    a = params(47); % 110 mix shape factor
    
    en(select) = en2(select) + (en1(select)-en2(select)).*rsw(phi(select),pi/2,0,a);
    en(~select) = en1(~select) + (en2(~select)-en1(~select)).*rsw(phi(~select),0,pi/2,a);
    
    % The Ksi degree is 129.4 and the Eta degree is 0 the one section of
    % rsw is not work. The Phi degree in this case use 3 section (k=3) to
    % interpolate the energy. 
    selectC = (eta == 0 & ksi > 129*pi/180 & ksi < 130*pi/180);
    
    k=3;
    vtheta=[0 0.292 1.277 pi/2];
    gamma=[entwist(1) params(48) params(49) entilt(1)];    
    
    aC = 0.5;
    
    en(selectC) = Ftwgb_stgb(phi(selectC), vtheta, gamma, k, aC);
 
end

function en = atgb110(ksi, eta, params)
%
    
    period = pi;
    
    en1 = stgb110(ksi, params);
    en2 = stgb110(period-ksi, params);
    
    eta(eta > period) = 2*period - eta(eta>period);
    
    a = params(44); % 110 ATGB shape factor
    
    en = Fatgb(eta, en1, en2, period, a);
    
end

function en = stgb110(ksi, params)
%

    % The same with stgb100() function
    k = 6;
    vtheta = [0 params(24) 70.65*pi/180 params(25) 129.5*pi/180 params(26) pi];
    gamma = [0 params(19) params(20) params(21) params(22) params(23) 0];
    
    a = 0.7;
    
    en = Ftwgb_stgb(ksi, vtheta, gamma, k, a);

end

function en = twist110(ksi, params)
%

    ksi(ksi>pi/2) = pi - ksi(ksi>pi/2);

    % The same with stgb100() function
    k=5;
    vtheta=[0 params(32) 50.5*pi/180 params(33) 70.5*pi/180 pi/2];
    gamma=[0 params(27) params(28) params(29) params(30) params(31)];    

    a = params(34);
    
    en = Ftwgb_stgb(ksi, vtheta, gamma, k, a);

end

%
% calculate boudnary energy for <111> set
%

function en = enMix111(ksi, eta, phi, params)
%

    entwist = twist111(ksi, params); 
    entilt = atgb111(ksi, eta, params);   
    
    en = Fmix111(phi, entwist, entilt, params);

end

function en = Fmix111(phi, entwist, entilt, params)
%

    a = params(50); % parabola coefficient
    b = a - 1;  % Ensures correct value at x = 1.

    x = phi/(pi/2); 

    % This one fit well enough with a simple one-parameter parabola that the
    % more complicated power laws in the other sets weren't needed.
    en = entwist + (entilt - entwist).*(a*x - b*x.^2) ;

end

function en = atgb111(ksi, eta, params)
%

    period = pi/3;
    
    en1 = stgb111(ksi, params);
    en2 = stgb111_eta60(ksi, params);
    
    eta(eta > period) = 2*period - eta(eta>period);

    a = params(45); %111 ATGB shape factor
    
    en = Fatgb(eta, en1, en2, period, a);

end

function en = stgb111(ksi, params)
%

    ksi(ksi > pi/3) = 2*pi/3 - ksi(ksi>pi/3);

    % The same with stgb100() function
    k=2;
    vtheta=[0 params(37) pi/3];
    gamma=[0 params(35) params(36)];
    
    a = 0.5;
    
    en = Ftwgb_stgb(ksi, vtheta, gamma, k, a);

end

function en = stgb111_eta60(ksi, params)
%

    ksi(ksi > pi/3) = 2*pi/3 - ksi(ksi>pi/3);

    % The same with stgb100() function
    k=2;
    vtheta=[0 params(40) pi/3];
    gamma=[0 params(38) params(39)];
    
    a = 0.5;
    
    en = Ftwgb_stgb(ksi, vtheta, gamma, k, a);

end

function en = twist111(ksi, params)
%

    ksi(ksi>pi/3) = 2*pi/3 - ksi(ksi>pi/3);
    
    % The same with stgb100() function
    k = 1;
    vtheta = [0 pi/3];
    gamma = [0 params(41)];
    
    a = params(42); %111 twist shape factor
    
    en = Ftwgb_stgb(ksi, vtheta, gamma, k, a);

end

%
%
% support function
%
%

function [P,Q]=convert2PQ(data)
%

    P(1,:)=[data(1) data(2) data(3)];
    P(2,:)=[data(4) data(5) data(6)];
    P(3,:)=[data(7) data(8) data(9)];
    Q(1,:)=[data(10) data(11) data(12)];
    Q(2,:)=[data(13) data(14) data(15)];
    Q(3,:)=[data(16) data(17) data(18)];
    
    P=P./repmat(sum(P.^2,2).^0.5,1,3);
    Q=Q./repmat(sum(Q.^2,2).^0.5,1,3);
    
end

function printGraph(zFit, z)
%

    h=figure(1);
    plot(zFit,z,'*');
    title('Energy Interpolation versus Energy Simulation');
    xlabel('Grain Boundary Energy Interpoation (J/m^{2})');
    ylabel('Grain Boundary Energy Simulation (J/m^{2})');
    drawnow;
    
end

% function params = readParams(E)
% %
%            
% %     switch E
% %         case 'Fe'
% %         params = [1.4023 0.2375	0.525 1.975	63.95 0.175 0.15 0.77402 0.69922 0.83711 0.77988 0.80156	0.56016	0.81621	1.0838	0.95527	0.78984	0.91367	0.4209	0.84355	0.17969	0.91641	0.72793	0.85625	0.62637	1.9059	2.717	0.48047	0.43418	0.46074	0.21758	0.44863	0.46465	0.9459	0.80898	0.18828	0.54063	0.80195	0.21934	0.55781	0.93281	0.86836	2.6916	-0.069141	0.086914	0.98437	-0.2	0.40234	0.75879	0.42207]';
% %         
% %         case 'W'
% % 
% %         params = [2.901	0.50506	0.49353	2.7026	4.5924	0.67813	0.097952	0.74844	0.68984	0.80859	0.72676	0.79551	0.5457	0.76328	1.0918	0.9791	0.79453	0.91465	0.37148	0.82383	0.20371	0.85977	0.70195	0.87812	0.55664	1.8666	2.6859	0.52187	0.43652	0.48066	0.25215	0.49336	0.425	0.97564	0.86406	0.20293	0.52285	0.85605	0.25469	0.5375	0.89395	1.2201	4.2859	-0.10781	-0.039844	1.0451	-0.31738	0.40449	0.72324	0.60078]';
% % 
% %         otherwise
% %         error('Undefined element');
% %     end
%     
% end

function en = Fatgb(eta, en1, en2, period, a)
%
% This is an energy function for asymmetric tilt boundaries
%

    % select only symmetric tilt energy larger than or equal to its invert
    % phase (by period)
    select = en1>=en2;
    
    en = zeros(size(eta));
    
    en(select) = en2(select) + (en1(select)-en2(select)).*rsw(eta(select),period,0,a);
    en(~select) = en1(~select) + (en2(~select)-en1(~select)).*rsw(eta(~select),0,period,a);
    
end

function en = Ftwgb_stgb(ksi, vtheta, gamma, k, a)
% the piecewise rsw function for ease change points
% the ksi of change points is defined in vtheta 
% the energy of change points is defined in gamma
% k is the number of pieces that separate by change points
%

    % define positive and negative infinity
    %
    pInf = 999999;
    nInf = -999999;

    for i=2:k+1
 
        vthetaI = Ifn(ksi, vtheta(i-1), vtheta(i), 0, 0); %interval of rsw segments
      
        rightGammaI = Ifn(gamma(i)-gamma(i-1), nInf, 0, 1, 0); %interval of negative gradients
        rightGamma = gamma(i)+(gamma(i-1)-gamma(i)).*rsw(ksi, vtheta(i), vtheta(i-1), a);
      
        leftGammaI = Ifn(gamma(i)-gamma(i-1), 0, pInf, 0, 1); %interval of positive gradients
        leftGamma = gamma(i-1)+(gamma(i)-gamma(i-1)).*rsw(ksi, vtheta(i-1), vtheta(i), a);

        E(i,:) = vthetaI.*(rightGammaI.*rightGamma + leftGammaI.*leftGamma);
    end
    en = sum(E);
    
end

function I = Ifn(keppa, a, b, aOpen, bOpen)
% the indicator function for keppa
% that assigns the different segment for rsw function
% aOpen and bOpen are boonlean for open interval a and b, respectively.
% set to 1 if it is open interval
% set to 0 if it is close interval
%

    if aOpen ka = a<keppa; else ka = a<=keppa; end
    if bOpen kb = b>keppa; else kb = b>=keppa; end
    I = ka&kb;
    
end

function frsw = rsw(theta, thetaMin, thetaMax, a)
% This function computes the value of Read-Shockley-Wolf function at theta.
% The rsw function is normalized to be 1.0 at thetaMax and 0.0 at thetaMin.
%
% theta is the angle at which to compute the function
% thetaMin is the starting angle of the interval which the second
% derivative at this point is infinity.
% thetaMax is the end angle of the interval which the second derivative at
% this point is zero.
% a is a parameter for defining the shape of the RSW function.
%
    
    x = (theta-thetaMin) ./ (thetaMax-thetaMin);
    
    % For nummerical problem of sin of zero, the small value is assigned. 
    x(x<0.00001) = 0.00001;
    frsw = sin(pi/2.*x).*(1-a.*log(sin(pi/2.*x)));
    
end

% function en = rsw(theta, theta1, theta2, a)
% % en = rsw(theta,theta1,theta2,a)
% %
% % This function computes the value of Read-Shockley-Wolf function at theta.
% % The rsw function is normalized to be 1.0 at theta2 and 0.0 at theta1.
% % 
% % theta             angle at which to compute the function
% % theta1            the starting angle of the interval
% % theta2            the end angle of the interval
% % a                 parameter defining the shape of the RSW function
% %
% 
%     dtheta = theta2 - theta1  ;     % Interval of angles where defined
%     theta = (theta-theta1)./dtheta*pi/2 ;    % Normalized angle
%     % The rest is the rsw function evaluation
%     sins = sin(theta) ;
%     xlogx = zeros(size(sins));
% 
%     % Cut off at small sins to avoid 0*infinity problem.  The proper limit is 0.
%     select = sins >= 0.000001;
%     xlogx(select) = sins(select).*log(sins(select));
% 
%     en = sins - a*xlogx;
%     
% end

%
%
% BKR functions, calculating the Scaffolding Data for <100>, <110>, and <111>
% sets, are duplicated from the publication version, {Ref}.
%
%

function geom = distances_to_set(P,Q,whichaxes,dismax)
% geom = distances_to_set(P,Q,whichaxes,dismax)
%
% Calculates the geometry parameters for a given grain boundary relative to
% a given set of axes.
%
% P and Q are rotation matrices giving the orientations of the two grains.
% The grain boundary normal is fixed at [1,0,0].
%
% whichaxes is one of '100', '110', or '111'
%
% dismax is an optional parameter specifying the maximum distance to
% include. It defaults to slightly less than 1, which is the largest
% allowable distance before some anomalies start to appear.
%
% Result geom is a 4xn matrix, where the rows are distance, ksi, eta, and
% phi. It keeps all of the hits up to a distance of dismax.
%
% distance is 2*sin(delta/2) where delta is the angle of closest approach
% between a misorientation axis and one of the axes.  Note there are 24
% representations of the rotations and 3, 6, or 4 equivalent high-symmetry
% axes, so it calculates as many as 144 distances.  But only ones below
% the cutoff dismax are kept.
%
% Once it's picked the closest approximation to the boundary for a given
% axis and coset element, it finds the parameters ksi, eta, phi defining
% that idealized boundary (since the axis is defined, it's a 3-space).
%
% These are:
% phi, the angle between the rotation axis and the boundary plane normal
% (taken as the mean of the normals represented in the two actual grain
% orientations, which works when dismax is less than 1)
%
% ksi, the misorientation angle
%
% eta, a parameter giving the second axis of the boundary plane normal in
% terms of specified directions ('dirs') perpendicular to each
% high-symmetry axis.


    if ~exist('dismax','var') || isempty(dismax),
        dismax = 0.999999 ;  % Force the distance to be strictly less than one, allowing for roundoff
    end % Note if dismax >= 1, you're likely to get a warning about m1 being singular.

    switch whichaxes
        case {'110'}
            % Define 110 axes, normalize
            axes = [1  1  1  1  0  0 ;
                    1 -1  0  0  1  1 ;
                    0  0  1 -1  1 -1 ]/sqrt(2) ;

            % Define a crystal direction perpendicular to each rotation axis.
            % The formalism demands that this be an axis of at least two-fold
            % symmetry.
            dirs = [0  0  0  0  1  1 ;
                    0  0  1  1  0  0 ;
                    1  1  0  0  0  0 ] ;
        case {'111'}

            % Define 111 axes, normalize
            axes = [1  1 -1 -1 ;
                    1 -1  1 -1 ;
                    1 -1 -1  1 ]/sqrt(3) ;

            dirs = [1  1  1  1 ;
                   -1  1  1 -1 ;
                    0  0  0  0 ]/sqrt(2) ;

        case {'100'}
            % Define 100 axes, normalize
            axes = [1  0  0 ;
                    0  1  0 ;
                    0  0  1 ] ;

            dirs = [0  0  1 ;
                    1  0  0 ;
                    0  1  0 ] ;
        otherwise
            error('Undefined axis set')
    end

    naxes  = size(axes,2) ;
    period = pi*naxes/6 ;

    %  Define the symmetry operators

    rotX90  = [ 1  0  0 ;     %  Rotation by +90 degrees around X axis
                0  0 -1 ;
                0  1  0 ] ;

    rotY90  = [ 0  0  1 ;     %  Rotation by +90 degrees around Y axis 
                0  1  0 ;
               -1  0  0 ] ;   

    rotZ90  = [ 0 -1  0 ;      %  Rotation by +90 degrees around Z axis
                1  0  0 ;
                0  0  1 ] ;

    rotZ90m = [ 0  1  0 ;      %  Rotation by -90 degrees around Z axis
               -1  0  0 ;
                0  0  1 ] ;

    % Create 24 symmetry equivalent variants of Q
    % This is the coset appropriate for the rotation convention where Q'*P
    % is the misorientation represented in the grain frame.  If you're
    % getting odd results, e.g. misorientations that you know are CSL are
    % coming out entirely wrong, you may be using the opposite convention;
    % try replacing P and Q with P' and Q'.
    V = cell(24,1);

    V{1}  = Q;    
    V{2}  = V{1}*rotX90 ;        % Rotate the vectors three times around X by +90 degrees
    V{3}  = V{2}*rotX90 ;
    V{4}  = V{3}*rotX90 ;

    for j = 1:12                   % Rotate three times around Y by +90 degrees
      V{j+4} = V{j}*rotY90 ;
    end

    for j = 1:4
      V{j+16} = V{j} *rotZ90;    % Rotate three times around Z by +90 degrees
      V{j+20} = V{j} *rotZ90m;    % Rotate three times around Z by -90 degrees   
    end


    % Preallocate all parameter lists at their maximum possible sizes.
    % Redundant representations will be removed at the end.
    distances = zeros(1,24*naxes);
    phis      = zeros(1,24*naxes);
    ksis      = zeros(1,24*naxes);
    etas      = zeros(1,24*naxes);

    thisindex = 0;  % Number of hits found so far

    % Step through all combinations of symmetrically-equivalent axes and coset
    % elements V{j}.
    for i = 1:naxes

        ax   = axes(:,i) ;    % ax is a high-symmetry axis   

        dir = dirs(:,i) ;       %  This is the pivot vector used to partition 
                                %  the rotation around axis "i"
        dir2 = cross(ax,dir);   %  Completing the orthonormal coordinate set.
                                %  theta1 and theta2 are defined in the plane
                                %  spanned by (dir,dir2).


        for j = 1:24    % For each symmetry-related variant of the second grain

            Q     = V{j} ; 
            R     = Q'*P ;   %  This rotates any vector in cube P into a vector in cube Q

            q = mat2quat(R) ;   % Calculation from here on out is much easier with quaternions.
            axi = q(2:4)'/sqrt(sum(q(2:4).^2)); % Normalized rotation axis
            psi = 2*acos(q(1)); % Rotation angle

            dotp  = axi*ax ;

            % Compute rotational distance from boundary P/Q to the rotation set "i" 
            % This formula produces 2*sin(delta/2), where delta is the angle of
            % closest approach.
            dis   = 2*sqrt(abs(1 - dotp*dotp))*sin(psi/2) ;

            if dis < dismax
                thisindex = thisindex + 1;

                theta = 2*atan(dotp*tan(psi/2)) ; % angle of rotation about ax that most closely approximates R

                % Compute the normal of the best-fitting GB in grain 1
                n1    = P(1,:)' ;
                n2    = Q(1,:)' ;

                RA = quat2mat([cos(theta/2);sin(theta/2)*ax]);
                % RA is the rotation about ax that most closely approximates R

                % From this point on we're dealing with the idealized rotation RA, not
                % the original rotation R.
                m1    = n1 + RA'*n2 ;

                % The next problem comes up only for very large distances,
                % which are normally cut off
                if norm(m1) < 0.000001
                    disp('m1 is singular!!!')
                end

                m1    = m1/norm(m1)  ;  % Halfway between the two normal vectors from the two grains
                m2    = RA*m1 ;   % And the same represented in the other grain

                % Compute the inclination angle for the common rotation axis
                phi   = real(acos(abs(m1'*ax))) ; % "real" because of numerical problems when they're exactly parallel

                % Partition the total rotation angle "theta"
 
                if abs(ax'*m1) > 0.9999      % Check if the best-fitting GB is pure twist
                    theta1 = - theta/2 ;   % eta is meaningless for a twist boundary.
                    theta2 =   theta/2 ;
                else

                    theta1 = atan2(dir2'*m1,dir'*m1);
                    theta2 = atan2(dir2'*m2,dir'*m2);
                    % It's projecting m1 and m2 into the plane normal to ax and
                    % then determining the rotation angles of them relative to
                    % dir.

                end        

                % Reduce both angles to interval (-period/2,period/2],
                % semi-open with a small numerical error.
                theta2  = theta2 - round(theta2/period)*period ;
                theta1  = theta1 - round(theta1/period)*period ;

                % This implements the semi-open interval in order to avoid an
                % annoying numerical problem where certain representations are
                % double-counted.
                if abs(theta2+period/2)<0.000001,
                    theta2 = theta2 + period;
                end
                if abs(theta1+period/2)<0.000001,
                    theta1 = theta1 + period;
                end

                % Since this is only being run on fcc elements, which are
                % centrosymmetric, and all dir vectors are 2-fold axes, then
                % the operations of swapping theta1 and theta2, and of
                % multilying both by -1, are symmetries for the energy
                % function. This lets us fold everything into a small right
                % triangle in (ksi,eta) space:
                ksi     = abs(theta2 - theta1) ;
                eta     = abs(theta2 + theta1) ;

                % And store them in the vectors
                distances(thisindex) = dis;
                ksis(thisindex)      = ksi;
                etas(thisindex)      = eta;
                phis(thisindex)      = phi;
            end
        end      
    end   

    % Dump the excess pre-allocated ones and sort the rest in order of distance
    [distances,sortindex] = sort(distances(1:thisindex));
    ksis = ksis(sortindex);
    etas = etas(sortindex);
    phis = phis(sortindex);

    % Clean up redundancy.  Double-counting the same representation of one
    % boundary messes up the weighting functions in weightedmeanenergy.m

    % First round everything to 1e-6, so that negligible numerical
    % differences are dropped
    distances = 1e-6*round(distances*1e6);
    ksis = 1e-6*round(ksis*1e6);
    etas = 1e-6*round(etas*1e6);
    phis = 1e-6*round(phis*1e6);

    % And finally create the 4 x thisindex array of geometrical parameters
    geom = unique([distances',ksis',etas',phis'],'rows')';

end

function q = mat2quat(m)
% q = mat2quat(m)
%
% Auxiliary function converts a rotation matrix, assumed orthonormal, into
% a unit quaternion.
    t = m(1,1)+m(2,2)+m(3,3);
    e0 = sqrt(1+t)/2;
    if t > -0.999999999
        e = [m(2,3)-m(3,2);m(3,1)-m(1,3);m(1,2)-m(2,1)]/(4*e0);
    else
        e0 = 0;
        e3 = sqrt(-(m(1,1)+m(2,2))/2);
        if abs(e3) > 2e-8   % Check for singularity, allowing numerical error
            e = [m(1,3)/(2*e3) ; m(2,3)/(2*e3) ; e3];
        else
            e1 = sqrt((m(1,1)+1)/2);
            if e1 ~= 0
                e = [e1;m(2,1)/(2*e1);0];
            else
                e = [0;1;0];
            end
        end
    end
    
    q = [e0;-e];
end

function m = quat2mat(q)
% m = quat2mat(q)
%
% Auxiliary function converts a quaternion into a rotation matrix with no
% assumption about normalization.
    e0 = q(1);
    e1 = q(2);
    e2 = q(3);
    e3 = q(4);

    m = [e0^2+e1^2-e2^2-e3^2 , 2*(e1*e2-e0*e3) , 2*(e1*e3+e0*e2); ...
          2*(e1*e2+e0*e3) , e0^2-e1^2+e2^2-e3^2 , 2*(e2*e3-e0*e1); ...
          2*(e1*e3-e0*e2) , 2*(e2*e3+e0*e1) , e0^2-e1^2-e2^2+e3^2 ]...
          /(e0^2+e1^2+e2^2+e3^2);
end
