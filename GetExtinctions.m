
function exs = GetExtinctions( lambda )


lambda = double(lambda);  % Cast to double

num_lambda = length(lambda);

%Tabulated Molar Extinction Coefficient for Hemoglobin in Water  https://omlc.org/spectra/hemoglobin/moaveni.html
 % These values for the molar extinction coefficient e in [cm-1/(moles/liter)] were compiled by Scott Prahl (prahl@ece.ogi.edu) using data from
        %
        % J.M. Schmitt, "Optical Measurement of Blood Oxygenation by Implantable Telemetry," Technical Report G558-15, Stanford."
        % M.K. Moaveni, "A Multiple Scattering Field Theory Applied to Whole Blood," Ph.D. dissertation, Dept. of Electrical Engineering, University of Washington, 1970.
        % To convert this data to absorbance A, multiply by the molar concentration and the pathlength. For example, if x is the number of grams per liter and a 1 cm cuvette is being used, then the absorbance is given by
        %
        %         (e) [(1/cm)/(moles/liter)] (x) [g/liter] (1) [cm]
        %   A =  ---------------------------------------------------
        %                           66,500 [g/mole]
        %
        % using 66,500 as the gram molecular weight of hemoglobin.
        % To convert this data to absorption coefficient in (cm-1), multiply by the molar concentration and 2.303,
        %
        % µa = (2.303) e (x g/liter)/(66,500 g Hb/mole)

vLambdaHbOHb = ...
   [630 680 4280;
    640	440	3640;
    650	380	3420;
    660	320	3200;
    670	320	3080;
    680	320	2960;
    690	280	2560;
    700	320	2160;
    710	340	1840;
    720	360	1520;
    730	400	1500;
    740	440	1520;
    750	520	1620;
    760	600	1720;
    770	660	1420;
    780	720	1120;
    790	760	1020;
    800	800	920;
    810	860	880;
    820	920	840;
    830	980	840;
    840	1040 840;
    850	1060 800;
    860	1080 840;
    870	1120 840;
    880	1160 840;
    890	1180 860;
    900	1200 880;
    910	1220 920;
    920	1240 880;
    930	1240 800;
    940	1200 800;
    950	1200 720];

        vLambdaHbOHb(:,2) = vLambdaHbOHb(:,2) * 2.303;
        vLambdaHbOHb(:,3) = vLambdaHbOHb(:,3) * 2.303;

exs = zeros(num_lambda, 2);

if(size(lambda,2) > size(lambda,1))
    lambda=lambda';
end

% HbO, Hb
ind=find(lambda>=250 & lambda<=1000);
exs(ind,1)=interp1(vLambdaHbOHb(:,1),vLambdaHbOHb(:,2),lambda(ind));
exs(ind,2)=interp1(vLambdaHbOHb(:,1),vLambdaHbOHb(:,3),lambda(ind));


end