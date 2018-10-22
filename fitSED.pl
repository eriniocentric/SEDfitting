
$wav = 1e-6*sequence(10000);

$Pi = 3.1415;
$clight = 2.99792e8; # speed of light in  meters
$planck = 6.6260755e-34; #planck's constant in m^2 kg /s
$boltzmann = 1.3806503e-23; # m^2 kg s-2 K


$lambda_0 = 200e-6;#200 micron
$beta = 1.5;
$alpha = 2.0;




$lambda_c = 1e-6*(3/4)*((26.68+6.246*$alpha)**(-2)+(1.905e-4+7.243e-5*$alpha)*$temp)**(-1);#in micron

$normpl = 1e27*(1-exp(-($lambda_0/$lambda_c)**$beta))*($lambda_c)**(-3)*(exp($planck*$clight/($lambda_c*$boltzmann*$temp))-1);

$bb = $modelflux = (1-exp(-($lambda_0/$wav)**$beta))*($clight/$wav)**(3)*(exp($planck*$clight/($wav*$boltzmann*$temp))-1)**(-1);
$modelflux = (1-exp(-($lambda_0/$wav)**$beta))*($clight/$wav)**(3)*(exp($planck*$clight/($wav*$boltzmann*$temp))-1)**(-1) + $normpl*($wav)**($alpha)*exp(-($wav/$lambda_c)**2);

dev '/xs';
line log10($wav),log10($modelflux),{Axis=>'logxy'}; 
hold;
line log10($wav),log10($bb),{Color=>'red'};
