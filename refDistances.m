% Get distances, azimuths, elevation angles between obs sites and vent/ref

obs = struct('SABAMF',[ -15.75687       -71.82677	5036 ],...
             'SABA_BASE2', [-15.7366733333	-71.827369	5212 ],...
             'RADAR', [-15.750128      -71.836807	5140.2]);
    
ref = [-15.7863888, -71.8525, 5913];
vent   = [-15.786744, -71.855919, 5911];

sites = fieldnames(obs);
spheroid = referenceEllipsoid('WGS 84');
for i = 1:numel(sites)
    site = obs.(sites{i});
    [azV,elevV,dV]= geodetic2aer(vent(1),vent(2),vent(3),site(1),site(2),site(3),spheroid);
    [azR,elevR,dR]= geodetic2aer(ref(1),ref(2),ref(3),site(1),site(2),site(3),spheroid);
    
    fprintf('\n%s:\tAzimuth\tElev\tSlant Dist (m)\n',sites{i})
    fprintf('Vent:\t%.2f\t%.2f\t%.2f\n',azV,elevV,dV)
    fprintf('Ref:\t%.2f\t%.2f\t%.2f\n',azR,elevR,dR)
end