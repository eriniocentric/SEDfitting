
($z,$mass,$masserr,$AV,$AVerr) = rcols 'Software/Masses/GDDSVIzK1234_GC_calz-masses.dat',1,3,4,14,15
dev 'Software/Figures/AV_vs_mass_calzetti.ps/vcps',{AspectRatio=>1};
$sel = which($masserr< 0.2 & $z > 1.3);
$i = which($masserr< 0.2);
points $mass($i), $AV($i),{XRange=>[8.2,11.5],XTitle=>'LogM',Ytitle=>'A\dv\u Calzetti et al (2000)'};
hold;
points $mass($sel), $AV($sel),{XRange=>[8.2,11.5],Color=>red};
legend ['All GDDS','z>1.3'],8.5,1.5, {Colour=>['black','red'],Symbol=>[17,17],TextShift=>-1.7};
close_window;


($z,$mass,$masserr,$AV,$AVerr) = rcols 'Software/Masses/GDDSVIzK1234_GC-masses.dat',1,3,4,14,15
dev 'Software/Figures/AV_vs_mass_pei.ps/vcps',{AspectRatio=>1};
$sel = which($masserr< 0.2 & $z > 1.3);
$i = which($masserr< 0.2);
points $mass($i), $AV($i),{XRange=>[8.2,11.5],XTitle=>'LogM',Ytitle=>'A\dv\u SMC Law Pei (1992)'};
hold;
points $mass($sel), $AV($sel),{XRange=>[8.2,11.5],Color=>red};
legend ['All GDDS','z>1.3'],8.5,1.5, {Colour=>['black','red'],Symbol=>[17,17],TextShift=>-1.7};
close_window;
