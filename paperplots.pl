
use PDL;
use lib "$ENV{HOME}/Software";

use PDL::NiceSlice;
use PDL::Graphics::PGPLOT;
use PGPLOT;
use PDL::Graphics::PGPLOTOptions ('set_pgplot_options');
set_pgplot_options('CharSize' => 2.3,'HardCH'=> 1.5,'HardLW'=>3,'AspectRatio'=>1,'TightLabels'=>1);
use KGB::Cosmology;
use KGB::PegUtils;
use KGB::SpecUtils;
use KGB::Dust;
use PDL::Fit::Polynomial;



($id,$spflag,$spectralclass,$z,$massKarl,$massKarl_err,$K,$Conf,$weight,$sfr2000,$sfrOII,$restUB,$gini,$assym,$fac,$facburst,$burstratio,$burstratio_avg,$burstratio_err,$bbtemp,$temp_avg,$temp_err,$massVIzK,$massVIzKerr,$zmaxVIzK,$massVIzK12,$massVIzK12err,$zmaxVIzK12,$massVIzK1234,$massVIzK1234err,$zmaxVIzK1234,$massVIzK1234_cce,$massVIzK1234_cceerr,$zmaxVIzK1234_cce,$rest3micronflux,$rest3micronflux_err,$massVIzK_nb,$massVIzK_nberr,$massVIzK12_nb,$massVIzK12_nberr,$massVIzK12_nb,$massVIzK12_nberr,$t,$t_err,$tburst,$tburst_err,$stellarburstratio,$stellarburstratio_err,$blendflag,$IRACflag,$sfr2000sed,$sfrOIIsed,$sfr2000obs,$sfrOIIobs) = rcols("$ENV{HOME}/Software/IRexcess/irexcess_gc.txt");

$mass = $massVIzK1234_cce;
$ix = which($burstratio_avg < 1e-9);
$burstratio_avg($ix) .= 1e-9;
$dlogburstratio_err = 0.4343 * $burstratio_err / $burstratio_avg;

$ix = which($rest3micronflux < 1e-9);
$rest3micronflux($ix) .= 1e-9;

$ix = which($stellarburstratio < 1e-9);
$stellarburstratio($ix) .= 1e-9;
$stellarburstratio_err($ix) .= 1e-9;
badmask($stellarburstratio->inplace,1e-9);
badmask($stellarburstratio_err->inplace,1e-9);
$dlogstellarburstratio_err = 0.4343 * $stellarburstratio_err / $stellarburstratio;
 
$opt = {Device=> "$ENV{HOME}/Software/paper/figures/br_vs_sfr.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->env(-1.6,2.1,-4.,1,{Axis=>'LogXY',XTitle=>'SFR / M\dsun\u/yr',YTitle=>'F\d\gl=3\gmm\u (Greybody/PAH)/ F\d\gl=3\gmm\u (Stellar)',AspectRatio=>1});

#sptype cut
$idx1 = which($K < 20.6  & $spflag == 0 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7);
$idx2 = which($K < 20.6  & $spflag == 1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7);
$idx3 = which($K < 20.6  & $spflag == 2 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7);
$idx4 = which($K < 20.6  & $spflag == -1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7);

$win->errb(log10($sfrOII($idx1)), log10($burstratio_avg($idx1)), undef, $dlogburstratio_err($idx1), {Symbol=>18,Colour=>'blue',SymbolSize=>2});
$win->errb(log10($sfrOII($idx2)), log10($burstratio_avg($idx2)), undef, $dlogburstratio_err($idx2), {Symbol=>18,Colour=>'green',SymbolSize=>2});
$win->errb(log10($sfrOII($idx3)), log10($burstratio_avg($idx3)), undef, $dlogburstratio_err($idx3), {Symbol=>18,Colour=>'red',SymbolSize=>2});
$win->errb(log10($sfrOII($idx4)), log10($burstratio_avg($idx4)), undef, $dlogburstratio_err($idx4), {Symbol=>18,Colour=>'orange',SymbolSize=>2});

$sfr = 10**(pdl(-4,-3,-2,1,0,1,2,3)); 

$Lexcess = 500; #calculated in circumstellar/disksize.pl
$fluxratio = $Lexcess * $sfr * 1e6;
$win->line(log10($sfr),log10($fluxratio));

$win->legend(["Red Evolved","Intermediate","Star Forming","Mixed"],-1.4,0.7,{Colour=>['Red','Green','Blue','Black'],Symbol=>[17,16,18,17],Charsize=>1.2,SymbolSize=>1.5,TextShift=>-0.3});

$win->close();

$opt = {Device=> "$ENV{HOME}/Software/paper/figures/3micron_vs_sfr.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->env(-1.5,2.,29.,33.5,{Axis=>'LogXY',XTitle=>'SFR (M\d\(2281)\u/yr)',YTitle=>'L\d3\gmm\u (Greybody/PAH) (W/\A)',AspectRatio=>1});

$Lsun = 3.826e26; #in Watts
#write flux in solar units
#$rest3micronflux /= $Lsun; #in Lsun
#$rest3micronflux_err /= $Lsun;#in Lsun
$dlogrest3micronflux_err = 0.4343 * $rest3micronflux_err/ $rest3micronflux;

#sptype cut
$idx1 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == 0 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);
$idx2 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == 1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);
$idx3 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == 2 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);
$idx4 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == -1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);

$win->errb(log10($sfrOII($idx1)), log10($rest3micronflux($idx1)), undef, $dlogrest3micronflux_err($idx1), {Symbol=>18,Colour=>'blue',SymbolSize=>2});
$win->errb(log10($sfrOII($idx2)), log10($rest3micronflux($idx2)), undef, $dlogrest3micronflux_err($idx2), {Symbol=>18,Colour=>'green',SymbolSize=>2});
$win->errb(log10($sfrOII($idx3)), log10($rest3micronflux($idx3)), undef, $dlogrest3micronflux_err($idx3), {Symbol=>18,Colour=>'red',SymbolSize=>2});
$win->errb(log10($sfrOII($idx4)), log10($rest3micronflux($idx4)), undef, $dlogrest3micronflux_err($idx4), {Symbol=>18,Colour=>'orange',SymbolSize=>2});

$UL3micronflux = log10($rest3micronflux_err) + $dlogrest3micronflux_err;

$win->errb(log10($sfr2000($idx1)), $UL3micronflux($idx1),undef,undef,{Symbol=>31,Colour=>'blue',SymbolSize=>2 & $blendflag==0});
$win->errb(log10($sfr2000($idx2)), $UL3micronflux($idx2),undef,undef,{Symbol=>31,Colour=>'green',SymbolSize=>2 & $blendflag==0});
$win->errb(log10($sfr2000($idx3)), $UL3micronflux($idx3),undef,undef,{Symbol=>31,Colour=>'red',SymbolSize=>2 & $blendflag==0});
$win->errb(log10($sfr2000($idx4)), $UL3micronflux($idx4),undef,undef,{Symbol=>31,Colour=>'orange',SymbolSize=>2  & $blendflag==0});


#sptype cut
$idx1 = which($K < 20.6  & $rest3micronflux < 0.5*$rest3micronflux_err & $spflag == 0 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7);
$idx2 = which($K < 20.6  & $rest3micronflux < 0.5*$rest3micronflux_err & $spflag == 1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7);
$idx3 = which($K < 20.6  & $rest3micronflux < 0.5*$rest3micronflux_err & $spflag == 2 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7);
$idx4 = which($K < 20.6  & $rest3micronflux < 0.5*$rest3micronflux_err & $spflag == -1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7);

$win->errb(log10($sfrOII($idx1)), $UL3micronflux($idx1), undef,undef, {Symbol=>31,Colour=>'blue',SymbolSize=>2});
$win->errb(log10($sfrOII($idx2)), $UL3micronflux($idx2), undef,undef, {Symbol=>31,Colour=>'green',SymbolSize=>2});
$win->errb(log10($sfrOII($idx3)), $UL3micronflux($idx3), undef,undef, {Symbol=>31,Colour=>'red',SymbolSize=>2});
$win->errb(log10($sfrOII($idx4)), $UL3micronflux($idx4), undef,undef, {Symbol=>31,Colour=>'orange',SymbolSize=>2});

#plot fit

$sfr = 10**(pdl(-4,-3,-2,1,0,1,2,3)); 
$r3mu = 2.62e31 * $sfr**(0.351);

$win->line(log10($sfr),log10($r3mu));

#use order of magnitude approximation Lexcess in Lband = 0.642 * (0.15) * sfr * lifetime
$Lexcess = 350; #calculated in circumstellar/disksize.pl
$flux3micron = $Lsun / 3e4 * $sfr * 1e6 * $Lexcess; #in units of W/A
$win->line(log10($sfr),log10($flux3micron),{Color=>green,Linestyle=>dashed});
$win->legend(["Red Evolved","Intermediate","Star Forming","Mixed"],-1.1,33.15,{Colour=>['Red','Green','Blue','orange'],Symbol=>[17,16,18,17],Charsize=>1.1,SymbolSize=>1.5,TextShift=>-1.65});

$win->close();

$opt = {Device=> "$ENV{HOME}/Software/paper/figures/br_vs_ssfr_all.ps/vcps",AspectRatio=>0.8};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->env(log10(3e-4),1,-3,log10(7),{Axis=>'LogXY',XTitle=>'log SSFR (Gyr\u-1\d)',YTitle=>'L\d3\gmm\u (Greybody/PAH) / L\d3\gmm\u (Stellar)'});
$ageburst = $t-$tburst;
#sptype cut 

#sptype cut 
$idx1 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == 0 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7  & $blendflag==0 );
$idx2 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == 1 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7  & $blendflag==0);
$idx3 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == 2 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7   & $blendflag==0);
$idx4 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == -1 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7   & $blendflag==0);

$ssfr2000 = $sfr2000 / (10**$massVIzK1234_cce);

#$win->errb(log10($ssfr2000($idx1))+9., log10($burstratio_avg($idx1)),undef,$dlogburstratio_err($idx1),{Symbol=>17,Colour=>'blue',SymbolSize=>2});
#$win->errb(log10($ssfr2000($idx2))+9., log10($burstratio_avg($idx2)),undef,$dlogburstratio_err($idx2),{Symbol=>17,Colour=>'green',SymbolSize=>2});
#$win->errb(log10($ssfr2000($idx3))+9., log10($burstratio_avg($idx3)),undef,$dlogburstratio_err($idx3),{Symbol=>17,Colour=>'red',SymbolSize=>2});
#$win->errb(log10($ssfr2000($idx4))+9., log10($burstratio_avg($idx4)),undef,$dlogburstratio_err($idx4),{Symbol=>17,Colour=>'orange',SymbolSize=>2});

$ssfrOII = $sfrOII / (10**$massVIzK1234_cce);

#sptype cut
$idx1 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == 0 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7  );
$idx2 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == 1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7  );
$idx3 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == 2 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7  );
$idx4 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == -1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 );


$win->errb(log10($ssfrOII($idx4))+9., log10($burstratio_avg($idx4)), undef, $dlogburstratio_err($idx4), {Symbol=>17,Colour=>'orange',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx2))+9., log10($burstratio_avg($idx2)), undef, $dlogburstratio_err($idx2), {Symbol=>16,Colour=>'green',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx3))+9., log10($burstratio_avg($idx3)), undef, $dlogburstratio_err($idx3), {Symbol=>17,Colour=>'red',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx1))+9., log10($burstratio_avg($idx1)), undef, $dlogburstratio_err($idx1), {Symbol=>18,Colour=>'blue',SymbolSize=>2});

#sptype cut 
$idx1 = which($K < 20.6 &  $burstratio_avg < 0.5*$burstratio_err  & $spflag == 0 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7  & $blendflag==0);
$idx2 = which($K < 20.6 &  $burstratio_avg < 0.5*$burstratio_err  & $spflag == 1 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7  & $blendflag==0);
$idx3 = which($K < 20.6 &  $burstratio_avg < 0.5*$burstratio_err  & $spflag == 2 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7   & $blendflag==0);
$idx4 = which($K < 20.6 &  $burstratio_avg < 0.5*$burstratio_err  & $spflag == -1 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7   & $blendflag==0);

#compute upper limit measurement of burst to plot for bad burst fits
$upburstratio = $burstratio_avg + $burstratio_err;

$ssfr2000 = $sfr2000 / (10**$massVIzK1234_cce);


#$win->errb(log10($ssfr2000($idx1))+9., log10($upburstratio($idx1)),undef,undef,{Symbol=>31,Colour=>'blue',SymbolSize=>2});
#$win->errb(log10($ssfr2000($idx2))+9., log10($upburstratio($idx2)),undef,undef,{Symbol=>31,Colour=>'green',SymbolSize=>2});
#$win->errb(log10($ssfr2000($idx3))+9., log10($upburstratio($idx3)),undef,undef,{Symbol=>31,Colour=>'red',SymbolSize=>2});
#$win->errb(log10($ssfr2000($idx4))+9., log10($upburstratio($idx4)),undef,undef,{Symbol=>31,Colour=>'orange',SymbolSize=>2});

$ssfrOII = $sfrOII / (10**$massVIzK1234_cce);


#sptype cut
$idx1 = which($K < 20.6 &  $burstratio_avg < $burstratio_err & $spflag == 0 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7  & $blendflag==0 );
$idx2 = which($K < 20.6 &  $burstratio_avg < $burstratio_err & $spflag == 1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7   & $blendflag==0);
$idx3 = which($K < 20.6 &  $burstratio_avg < $burstratio_err & $spflag == 2 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7   & $blendflag==0);
$idx4 = which($K < 20.6 &  $burstratio_avg < $burstratio_err & $spflag == -1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7  & $blendflag==0);

$win->errb(log10($ssfrOII($idx4))+9., log10($upburstratio($idx4)), undef,undef, {Symbol=>31,Colour=>'orange',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx2))+9., log10($upburstratio($idx2)), undef,undef, {Symbol=>31,Colour=>'green',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx3))+9., log10($upburstratio($idx3)), undef,undef, {Symbol=>31,Colour=>'red',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx1))+9., log10($upburstratio($idx1)), undef,undef, {Symbol=>31,Colour=>'blue',SymbolSize=>2});

$ssfr = 10**(pdl(-4,-3,-2,1,0,1,2,3)); 
#$model = 0.33834126 * log10($ssfr) - 0.17239746;
$model = 0.350 * log10($ssfr) - 0.197;
$win->line(log10($ssfr),$model);

#toy model
$masstolightK = 0.3;#calculated in plotSED_cc.pl
$Lexcess = 350; #calculated in disksize.pl
#$toymodel =  $masstolightexcess**(-1) * ($ssfr * $masstolightK) * 0.01;
$toymodel = $Lexcess * $ssfr * 0.001;
#$win->line(log10($ssfr),log10($toymodel),{Color=>green,LineStyle=>dashed});
$win->legend(["Red Evolved","Intermediate","Star Forming","Mixed"],-3.1,0.5,{Colour=>['Red','Green','Blue','Orange'],Symbol=>[17,16,18,17],Charsize=>1.1,SymbolSize=>3,TextShift=>-1.65});

$win->close();



$opt = {Device=> "$ENV{HOME}/Software/paper/figures/br_vs_ssfrsed_all.ps/vcps",AspectRatio=>0.8};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->env(log10(3e-4),1,-3,log10(7),{Axis=>'LogXY',XTitle=>'log SSFR (Gyr\u-1\d)',YTitle=>'L\d3\gmm\u (Greybody/PAH) / L\d3\gmm\u (Stellar)'});
$ageburst = $t-$tburst;
#sptype cut 

#sptype cut 
$idx1 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == 0 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7  & $blendflag==0 );
$idx2 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == 1 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7  & $blendflag==0);
$idx3 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == 2 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7   & $blendflag==0);
$idx4 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == -1 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7   & $blendflag==0);

$ssfr2000 = $sfr2000sed / (10**$massVIzK1234_cce);

#$win->errb(log10($ssfr2000($idx1))+9., log10($burstratio_avg($idx1)),undef,$dlogburstratio_err($idx1),{Symbol=>17,Colour=>'blue',SymbolSize=>2});
#$win->errb(log10($ssfr2000($idx2))+9., log10($burstratio_avg($idx2)),undef,$dlogburstratio_err($idx2),{Symbol=>17,Colour=>'green',SymbolSize=>2});
#$win->errb(log10($ssfr2000($idx3))+9., log10($burstratio_avg($idx3)),undef,$dlogburstratio_err($idx3),{Symbol=>17,Colour=>'red',SymbolSize=>2});
#$win->errb(log10($ssfr2000($idx4))+9., log10($burstratio_avg($idx4)),undef,$dlogburstratio_err($idx4),{Symbol=>17,Colour=>'orange',SymbolSize=>2});

$ssfrOII = $sfrOIIsed / (10**$massVIzK1234_cce);

#sptype cut
$idx1 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == 0 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7  );
$idx2 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == 1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7  );
$idx3 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == 2 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7  );
$idx4 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $spflag == -1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 );


$win->errb(log10($ssfrOII($idx4))+9., log10($burstratio_avg($idx4)), undef, $dlogburstratio_err($idx4), {Symbol=>17,Colour=>'orange',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx2))+9., log10($burstratio_avg($idx2)), undef, $dlogburstratio_err($idx2), {Symbol=>16,Colour=>'green',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx3))+9., log10($burstratio_avg($idx3)), undef, $dlogburstratio_err($idx3), {Symbol=>17,Colour=>'red',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx1))+9., log10($burstratio_avg($idx1)), undef, $dlogburstratio_err($idx1), {Symbol=>18,Colour=>'blue',SymbolSize=>2});

#sptype cut 
$idx1 = which($K < 20.6 &  $burstratio_avg < 0.5*$burstratio_err  & $spflag == 0 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7  & $blendflag==0);
$idx2 = which($K < 20.6 &  $burstratio_avg < 0.5*$burstratio_err  & $spflag == 1 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7  & $blendflag==0);
$idx3 = which($K < 20.6 &  $burstratio_avg < 0.5*$burstratio_err  & $spflag == 2 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7   & $blendflag==0);
$idx4 = which($K < 20.6 &  $burstratio_avg < 0.5*$burstratio_err  & $spflag == -1 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7   & $blendflag==0);

#compute upper limit measurement of burst to plot for bad burst fits
$upburstratio = $burstratio_avg + $burstratio_err;

#$win->errb(log10($ssfr2000($idx1))+9., log10($upburstratio($idx1)),undef,undef,{Symbol=>31,Colour=>'blue',SymbolSize=>2});
#$win->errb(log10($ssfr2000($idx2))+9., log10($upburstratio($idx2)),undef,undef,{Symbol=>31,Colour=>'green',SymbolSize=>2});
#$win->errb(log10($ssfr2000($idx3))+9., log10($upburstratio($idx3)),undef,undef,{Symbol=>31,Colour=>'red',SymbolSize=>2});
#$win->errb(log10($ssfr2000($idx4))+9., log10($upburstratio($idx4)),undef,undef,{Symbol=>31,Colour=>'orange',SymbolSize=>2});

#sptype cut
$idx1 = which($K < 20.6 &  $burstratio_avg < $burstratio_err & $spflag == 0 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7  & $blendflag==0 );
$idx2 = which($K < 20.6 &  $burstratio_avg < $burstratio_err & $spflag == 1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7   & $blendflag==0);
$idx3 = which($K < 20.6 &  $burstratio_avg < $burstratio_err & $spflag == 2 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7   & $blendflag==0);
$idx4 = which($K < 20.6 &  $burstratio_avg < $burstratio_err & $spflag == -1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7  & $blendflag==0);

$win->errb(log10($ssfrOII($idx4))+9., log10($upburstratio($idx4)), undef,undef, {Symbol=>31,Colour=>'orange',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx2))+9., log10($upburstratio($idx2)), undef,undef, {Symbol=>31,Colour=>'green',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx3))+9., log10($upburstratio($idx3)), undef,undef, {Symbol=>31,Colour=>'red',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx1))+9., log10($upburstratio($idx1)), undef,undef, {Symbol=>31,Colour=>'blue',SymbolSize=>2});

$ssfr = 10**(pdl(-4,-3,-2,1,0,1,2,3)); 
#$model = 0.33834126 * log10($ssfr) - 0.17239746;
$model = 0.350 * log10($ssfr) - 0.197;
$win->line(log10($ssfr),$model);

#toy model
$masstolightK = 0.3;#calculated in plotSED_cc.pl
$Lexcess = 350; #calculated in disksize.pl
#$toymodel =  $masstolightexcess**(-1) * ($ssfr * $masstolightK) * 0.01;
$toymodel = $Lexcess * $ssfr * 0.001;
#$win->line(log10($ssfr),log10($toymodel),{Color=>green,LineStyle=>dashed});
$win->legend(["Red Evolved","Intermediate","Star Forming","Mixed"],-3.1,0.5,{Colour=>['Red','Green','Blue','Orange'],Symbol=>[17,16,18,17],Charsize=>1.1,SymbolSize=>3,TextShift=>-1.65});

$win->close();



$opt = {Device=> "$ENV{HOME}/Software/paper/figures/IRAC_vs_VIzK_UBcuts.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
#$win = PDL::Graphics::PGPLOT::Window->new({Device=>'/xs'});
$win->env(8.,12.,8.,12.,{Axis=>'LogXY',XTitle=>'VIzK log\d10\uM  (M\d\(2281)\u) ',YTitle=>'VIzK+IRAC log\d10\uM  (M\d\(2281)\u)',Justify=>1});

$mass1 = $massVIzK1234_cce;
$mass2 = $massVIzK;
$mass1_err = $massVIzK1234_cceerr;
$mass2_err = $massVIzKerr;

$idx = which($restUB < -0.1 & $massVIzKerr < 1 & $Conf > 1 & $Conf < 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'blue',Symbol=>18,SymbolSize=>1.5});
$idx = which($restUB > -0.1 & $restUB < 0.1 & 1 & $massVIzKerr < 1 & $Conf > 1 & $Conf < 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'green',Symbol=>16,SymbolSize=>1.5});
$idx = which($restUB > 0.1 & $massVIzKerr < 1 & $Conf > 1 & $Conf < 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'red',Symbol=>17,SymbolSize=>1.5});

$idx = which($restUB < -0.1 & $massVIzKerr < 1 & $Conf < 1 & $Conf > 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'blue',Symbol=>'star',SymbolSize=>1.5});
$idx = which($restUB > -0.1 & $restUB < 0.1 & $massVIzKerr < 1 & $Conf < 1 & $Conf > 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'green',Symbol=>'square',SymbolSize=>1.5});
$idx = which($restUB > 0.1 & $massVIzKerr < 1 & $Conf < 1 & $Conf > 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'red',Symbol=>'circle',SymbolSize=>1.5});
$win->line(sequence(20),sequence(20),{Colour=>'Black'});
$win->legend([' 0.1 < (U-B)\drest\u','-0.1 < (U-B)\drest\u < 0.1','-0.1 > (U-B)\drest\u'],10.3,8.8,{Colour=>['Red','Green','Blue'],Symbol=>[17,16,18],Charsize=>1.1,SymbolSize=>1.5,TextShift=>-0.3});
$win->close();


$opt = {Device=> "$ENV{HOME}/Software/paper/figures/IRAC_vs_VIzK_sptypes.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
#$win = PDL::Graphics::PGPLOT::Window->new({Device=>'/xs'});
$win->env(8.,12.,8.,12.,{Axis=>'LogXY',XTitle=>'VIzK log\d10\uM  (M\d\(2281)\u) ',YTitle=>'VIzK+IRAC log\d10\uM  (M\d\(2281)\u)',Justify=>1});

$mass1 = $massVIzK1234_cce;
$mass2 = $massVIzK;
$mass1_err = $massVIzK1234_cceerr;
$mass2_err = $massVIzKerr;

$idx = which($spflag == 0 & $massVIzKerr < 1 & $Conf > 1 & $Conf < 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'blue',Symbol=>18,SymbolSize=>1.5});
$idx = which($spflag == 1 & $massVIzKerr < 1 & $Conf > 1 & $Conf < 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'green',Symbol=>16,SymbolSize=>1.5});
$idx = which($spflag == 2 & $massVIzKerr < 1 & $Conf > 1 & $Conf < 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'red',Symbol=>17,SymbolSize=>1.5});
$idx = which($spflag == -1 & $massVIzKerr < 1 & $Conf > 1 & $Conf < 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'orange',Symbol=>17,SymbolSize=>1.5});

$idx = which($spflag == 0 & $massVIzKerr < 1 & $Conf < 1 & $Conf > 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'blue',Symbol=>'star',SymbolSize=>1.5});
$idx = which($spflag == 1 & $massVIzKerr < 1 & $Conf < 1 & $Conf > 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'green',Symbol=>'square',SymbolSize=>1.5});
$idx = which($spflag == 2 & $massVIzKerr < 1 & $Conf < 1 & $Conf > 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'red',Symbol=>'circle',SymbolSize=>1.5});
$idx = which($spflag == -1 & $massVIzKerr < 1 & $Conf < 1 & $Conf > 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'orange',Symbol=>'circle',SymbolSize=>1.5});


$win->line(sequence(20),sequence(20),{Colour=>'Orange'});
$win->legend(["Red Evolved","Intermediate","Star Forming","Mixed"],10.3,8.8,{Colour=>['Red','Green','Blue','Orange'],Symbol=>[17,16,18,17],Charsize=>1.1,SymbolSize=>1.5,TextShift=>-0.3});
$win->close();


$opt = {Device=> "$ENV{HOME}/Software/paper/figures/Burst_vs_noburst_UBcuts.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
#$win = PDL::Graphics::PGPLOT::Window->new({Device=>'/xs'});
$win->env(8.,12.,8.,12.,{Axis=>'LogXY',XTitle=>'One component log\d10\uM  (M\d\(2281)\u) ',YTitle=>'Two component log\d10\uM  (M\d\(2281)\u)',Justify=>1});

$mass1 = $massVIzK12;
$mass2 = $massVIzK12_nb;
$mass1_err = $massVIzK12err;
$mass2_err = $massVIzK12_nberr;

$idx = which($restUB < -0.1 & $massVIzKerr < 1 & $Conf > 1 & $Conf < 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'blue',Symbol=>18,SymbolSize=>1.5});
$idx = which($restUB > -0.1 & $restUB < 0.1 & 1 & $massVIzKerr < 1 & $Conf > 1 & $Conf < 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'green',Symbol=>16,SymbolSize=>1.5});
$idx = which($restUB > 0.1 & $massVIzKerr < 1 & $Conf > 1 & $Conf < 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'red',Symbol=>17,SymbolSize=>1.5});

$idx = which($restUB < -0.1 & $massVIzKerr < 1 & $Conf < 1 & $Conf > 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'blue',Symbol=>'star',SymbolSize=>1.5});
$idx = which($restUB > -0.1 & $restUB < 0.1 & $massVIzKerr < 1 & $Conf < 1 & $Conf > 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'green',Symbol=>'square',SymbolSize=>1.5});
$idx = which($restUB > 0.1 & $massVIzKerr < 1 & $Conf < 1 & $Conf > 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'red',Symbol=>'circle',SymbolSize=>1.5});
$win->line(sequence(20),sequence(20),{Colour=>'Orange'});
$win->legend([' 0.1 < (U-B)\drest\u','-0.1 < (U-B)\drest\u < 0.1','-0.1 > (U-B)\drest\u'],10.3,8.8,{Colour=>['Red','Green','Blue'],Symbol=>[17,16,18],Charsize=>1.1,SymbolSize=>1.5,TextShift=>-0.3});
$win->close();


$opt = {Device=> "$ENV{HOME}/Software/paper/figures/Burst_vs_noburst_sptypes.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
#$win = PDL::Graphics::PGPLOT::Window->new({Device=>'/xs'});
$win->env(8.,12.,8.,12.,{Axis=>'LogXY',XTitle=>'One component log\d10\uM  (M\d\(2281)\u) ',YTitle=>'Two component log\d10\uM  (M\d\(2281)\u)',Justify=>1});

$mass1 = $massVIzK12;
$mass2 = $massVIzK12_nb;
$mass1_err = $massVIzK12err;
$mass2_err = $massVIzK12_nberr;

$idx = which($spflag == 0 & $massVIzKerr < 1 & $Conf > 1 & $Conf < 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'blue',Symbol=>18,SymbolSize=>1.5});
$idx = which($spflag == 1 & $massVIzKerr < 1 & $Conf > 1 & $Conf < 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'green',Symbol=>16,SymbolSize=>1.5});
$idx = which($spflag == 2 & $massVIzKerr < 1 & $Conf > 1 & $Conf < 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'red',Symbol=>17,SymbolSize=>1.5});
$idx = which($spflag == -1 & $massVIzKerr < 1 & $Conf > 1 & $Conf < 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'orange',Symbol=>17,SymbolSize=>1.5});

$idx = which($spflag == 0 & $massVIzKerr < 1 & $Conf < 1 & $Conf > 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'blue',Symbol=>'star',SymbolSize=>1.5});
$idx = which($spflag == 1 & $massVIzKerr < 1 & $Conf < 1 & $Conf > 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'green',Symbol=>'square',SymbolSize=>1.5});
$idx = which($spflag == 2 & $massVIzKerr < 1 & $Conf < 1 & $Conf > 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'red',Symbol=>'circle',SymbolSize=>1.5});
$idx = which($spflag == -1 & $massVIzKerr < 1 & $Conf < 1 & $Conf > 7);
$win->errb($mass2($idx),$mass1($idx),$mass2_err($idx),$mass1_err($idx),{Col=>'orange',Symbol=>'circle',SymbolSize=>1.5});


$win->line(sequence(20),sequence(20),{Colour=>'Orange'});
$win->legend(["Red Evolved","Intermediate","Star Forming","Mixed"],10.3,8.8,{Colour=>['Red','Green','Blue','Orange'],Symbol=>[17,16,18,17],Charsize=>1.1,SymbolSize=>1.5,TextShift=>-0.3});
$win->close();

#chi square paper plots


$opt = {Device=> "$ENV{HOME}/Software/paper/figures/reducedchi2_GB_vs_noGB_sptypes.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
#$win = PDL::Graphics::PGPLOT::Window->new({Device=>'/xs'});
$win->env(0,28,0,28,{Xtitle=>'\gx\u2\d - Stellar Only', YTitle=>'\gx\u2\d - with Greybody/PAH'});

($id,$chi2,$chi2err,$chi2orig,$chi2origerr,$nadata,$chi2cce,$chi2cceerr,$chi2origcce,$chi2origcceerr,$ndata,$spflag,$spectralclass,$z,$massKarl,$massKarl_err,$K,$Conf,$weight,$sfr2000,$sfrOII,$restUB,$gini,$assym,$fac,$facburst,$burstratio,$burstratio_avg,$burstratio_err,$bbtemp,$temp_avg,$temp_err,$massVIzK,$massVIzKerr,$zmaxVIzK,$massVIzK12,$massVIzK12err,$zmaxVIzK12,$massVIzK1234,$massVIzK1234err,$zmaxVIzK1234,$massVIzK1234_cce,$massVIzK1234_cceerr,$zmaxVIzK1234_cce,$rest4micronflux,$rest4micronflux_err,$massVIzK_nb,$massVIzK_nberr,$massVIzK12_nb,$massVIzK12_nberr,$massVIzK12_nb,$massVIzK12_nberr,$t,$t_err,$tburst,$tburst_err) = rcols("$ENV{HOME}/Software/chi2/chi2excess.dat");

#sptype cut
$idx1 = which($K < 20.6 & $spflag == 0 & ($Conf > 1 | $z<9) & $z>0.01 & $Conf < 7 );
$idx2 = which($K < 20.6 & $spflag == 1 & ($Conf > 1 | $z<9) & $z>0.01 & $Conf < 7 );
$idx3 = which($K < 20.6 & $spflag == 2 & ($Conf > 1 | $z<9) & $z>0.01 & $Conf < 7 ) ;
$idx4 = which($K < 20.6 & $spflag == -1 & ($Conf > 1 | $z<9) & $z>0.01 & $Conf < 7 );

#convert monte carlo chi squares to reduced chi squares
$chi2cce /= $ndata;
$chi2 /= $ndata;
$chi2cceerr /= $ndata;
$chi2err /= $ndata;

$win->errb($chi2($idx4), $chi2cce($idx4),$chi2err($idx4),$chi2cceerr($idx4),{Symbol=>17,Colour=>'orange',SymbolSize=>1.3,LineWidth=>1});
$win->errb($chi2($idx2), $chi2cce($idx2),$chi2err($idx2),$chi2cceerr($idx2),{Symbol=>16,Colour=>'green',SymbolSize=>1.3,LineWidth=>1});
$win->errb($chi2($idx3), $chi2cce($idx3),$chi2err($idx3),$chi2cceerr($idx3),{Symbol=>17,Colour=>'red',SymbolSize=>1.3,LineWidth=>1});
$win->errb($chi2($idx1), $chi2cce($idx1),$chi2err($idx1),$chi2cceerr($idx1),{Symbol=>18,Colour=>'blue',SymbolSize=>1.3,LineWidth=>1});

$win->line(sequence(200),sequence(200),{Colour=>'Black',LineStyle=>Dashed,Linewidth=>1.5});
$win->legend(["Red Evolved","Intermediate","Star Forming","Mixed"],4,25,{Colour=>['Red','Green','Blue','Orange'],Symbol=>[17,16,18,17],Charsize=>1.1,SymbolSize=>1.5,TextShift=>-1.7});
$win->close();

$opt = {Device=> "$ENV{HOME}/Software/paper/figures/chi2_GB_vs_noGB_sptypes.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
#$win = PDL::Graphics::PGPLOT::Window->new({Device=>'/xs'});
$win->env(0,175,0,175,{Xtitle=>'\gx\u2\d - Stellar', YTitle=>'\gx\u2\d - Stellar+GB'});

$chi2origcce *= $ndata;
$chi2orig *= $ndata;
$chi2origcceerr *= $ndata;
$chi2origerr *= $ndata;

$win->errb($chi2($idx1), $chi2cce($idx1),$chi2err($idx1),$chi2cceerr($idx1),{Symbol=>18,Colour=>'blue',SymbolSize=>1.7});
$win->errb($chi2($idx2), $chi2cce($idx2),$chi2err($idx2),$chi2cceerr($idx2),{Symbol=>16,Colour=>'green',SymbolSize=>1.7});
$win->errb($chi2($idx3), $chi2cce($idx3),$chi2err($idx3),$chi2cceerr($idx3),{Symbol=>17,Colour=>'red',SymbolSize=>1.7});
$win->errb($chi2($idx4), $chi2cce($idx4),$chi2err($idx4),$chi2cceerr($idx4),{Symbol=>17,Colour=>'orange',SymbolSize=>1.7});
$win->line(sequence(200),sequence(200),{Colour=>'Orange'});
$win->legend(["Red Evolved","Intermediate","Star Forming","Mixed"],5,165,{Colour=>['Red','Green','Blue','Orange'],Symbol=>[17,16,18,17],Charsize=>1.1,SymbolSize=>1.5,TextShift=>-0.3});
$win->close();



$opt = {Device=> "$ENV{HOME}/Software/paper/figures/chi2_GB_vs_noGB_sptypeszoom.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
#$win = PDL::Graphics::PGPLOT::Window->new({Device=>'/xs'});
$win->env(0,55,0,55,{Xtitle=>'\gx\u2\d - Stellar', YTitle=>'\gx\u2\d - Stellar+GB'});

$win->errb($chi2($idx1), $chi2cce($idx1),$chi2err($idx1),$chi2cceerr($idx1),{Symbol=>18,Colour=>'blue',SymbolSize=>1.7});
$win->errb($chi2($idx2), $chi2cce($idx2),$chi2err($idx2),$chi2cceerr($idx2),{Symbol=>16,Colour=>'green',SymbolSize=>1.7});
$win->errb($chi2($idx3), $chi2cce($idx3),$chi2err($idx3),$chi2cceerr($idx3),{Symbol=>17,Colour=>'red',SymbolSize=>1.7});
$win->errb($chi2($idx4), $chi2cce($idx4),$chi2err($idx4),$chi2cceerr($idx4),{Symbol=>17,Colour=>'orange',SymbolSize=>1.7});
$win->line(sequence(200),sequence(200),{Colour=>'Orange'});
$win->legend(["Red Evolved","Intermediate","Star Forming","Mixed"],4,50,{Colour=>['Red','Green','Blue','Orange'],Symbol=>[17,16,18,17],Charsize=>1.1,SymbolSize=>1.5,TextShift=>-0.3});
$win->close();


$opt = {Device=> "$ENV{HOME}/Software/paper/figures/logchi2_GB_vs_noGB_sptypes.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
#$win = PDL::Graphics::PGPLOT::Window->new({Device=>'/xs'});
$win->env(0,3,0,3,{Xtitle=>'log(\gx\u2\d) - Stellar', YTitle=>'log(\gx\u2\d) - Stellar+GB'});

$logchi2orig = log10($chi2orig);
$logchi2origcce = log10($chi2origcce);
$logchi2origerr = 0.4343 * $chi2origerr / $chi2orig;
$logchi2origcceerr = 0.4343 * $chi2origcceerr / $chi2origcce;

$win->errb($logchi2orig($idx1), $logchi2origcce($idx1),$logchi2origerr($idx1),$logchi2origcceerr($idx1),{Symbol=>18,Colour=>'blue',SymbolSize=>1.7});
$win->errb($logchi2orig($idx2), $logchi2origcce($idx2),$logchi2origerr($idx2),$logchi2origcceerr($idx2),{Symbol=>16,Colour=>'green',SymbolSize=>1.7});
$win->errb($logchi2orig($idx3), $logchi2origcce($idx3),$logchi2origerr($idx3),$logchi2origcceerr($idx3),{Symbol=>17,Colour=>'red',SymbolSize=>1.7});
$win->errb($logchi2orig($idx4), $logchi2origcce($idx4),$logchi2origerr($idx4),$logchi2origcceerr($idx4),{Symbol=>17,Colour=>'orange',SymbolSize=>1.7});
$win->line(sequence(10)-5,sequence(10)-5,{Colour=>'Orange'});
$win->legend(["Red Evolved","Intermediate","Star Forming","Mixed"],0.2,2.7,{Colour=>['Red','Green','Blue','Orange'],Symbol=>[17,16,18,17],Charsize=>1.1,SymbolSize=>1.5,TextShift=>-0.3});
$win->close();


($id,$spflag,$spectralclass,$z,$massKarl,$massKarl_err,$K,$Conf,$weight,$sfr2000,$sfrOII,$restUB,$gini,$assym,$fac,$facburst,$burstratio,$burstratio_avg,$burstratio_err,$bbtemp,$temp_avg,$temp_err,$massVIzK,$massVIzKerr,$zmaxVIzK,$massVIzK12,$massVIzK12err,$zmaxVIzK12,$massVIzK1234,$massVIzK1234err,$zmaxVIzK1234,$massVIzK1234_cce,$massVIzK1234_cceerr,$zmaxVIzK1234_cce,$rest3micronflux,$rest3micronflux_err,$massVIzK_nb,$massVIzK_nberr,$massVIzK12_nb,$massVIzK12_nberr,$massVIzK12_nb,$massVIzK12_nberr,$t,$t_err,$tburst,$tburst_err,$stellarburstratio,$stellarburstratio_err,$blendflag,$IRACflag) = rcols("$ENV{HOME}/Software/IRexcess/irexcess_gc.txt");

$opt = {Device=> "$ENV{HOME}/Software/paper/figures/br_vs_ssfr_masscuts.ps/vcps",AspectRatio=>0.8};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->env(log10(3e-4),1,-3,log10(7),{Axis=>'LogXY',XTitle=>'log SSFR (Gyr\u-1\d)',YTitle=>'L\d3\gmm\u (Greybody/PAH) / L\d3\gmm\u (Stellar)'});
$ageburst = $t-$tburst;
#sptype cut 

#sptype cut 
$idx1 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $massVIzK1234_cce < 10.2 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7  & $blendflag==0 );
$idx2 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $massVIzK1234_cce > 10.2 & $massVIzK1234_cce < 10.8 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7  & $blendflag==0);
$idx3 = which($K < 20.6 &  $burstratio_avg > 0.5*$burstratio_err  & $massVIzK1234_cce > 10.8 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7   & $blendflag==0);

$ssfrOII = $sfrOII / (10**$massVIzK1234_cce);


$win->errb(log10($ssfrOII($idx2))+9., log10($burstratio_avg($idx2)), undef, $dlogburstratio_err($idx2), {Symbol=>16,Colour=>'green',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx3))+9., log10($burstratio_avg($idx3)), undef, $dlogburstratio_err($idx3), {Symbol=>17,Colour=>'red',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx1))+9., log10($burstratio_avg($idx1)), undef, $dlogburstratio_err($idx1), {Symbol=>18,Colour=>'blue',SymbolSize=>2});

#compute upper limit measurement of burst to plot for bad burst fits
$upburstratio = $burstratio_avg + $burstratio_err;

#sptype cut
$idx1 = which($K < 20.6 &  $burstratio_avg < $burstratio_err & $massVIzK1234_cce < 10.2 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7  & $blendflag==0 );
$idx2 = which($K < 20.6 &  $burstratio_avg < $burstratio_err & $massVIzK1234_cce > 10.2 & $massVIzK1234_cce < 10.8 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7   & $blendflag==0);
$idx3 = which($K < 20.6 &  $burstratio_avg < $burstratio_err & $massVIzK1234_cce > 10.8 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7   & $blendflag==0);

$win->errb(log10($ssfrOII($idx2))+9., log10($upburstratio($idx2)), undef,undef, {Symbol=>31,Colour=>'green',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx3))+9., log10($upburstratio($idx3)), undef,undef, {Symbol=>31,Colour=>'red',SymbolSize=>2});
$win->errb(log10($ssfrOII($idx1))+9., log10($upburstratio($idx1)), undef,undef, {Symbol=>31,Colour=>'blue',SymbolSize=>2});

$ssfr = 10**(pdl(-4,-3,-2,1,0,1,2,3)); 
#$model = 0.33834126 * log10($ssfr) - 0.17239746;
$model = 0.350 * log10($ssfr) - 0.197;
$win->line(log10($ssfr),$model);

#toy model
$masstolightK = 0.3;#calculated in plotSED_cc.pl
$Lexcess = 350; #calculated in disksize.pl
#$toymodel =  $masstolightexcess**(-1) * ($ssfr * $masstolightK) * 0.01;
$toymodel = $Lexcess * $ssfr * 0.001;
#$win->line(log10($ssfr),log10($toymodel),{Color=>green,LineStyle=>dashed});
$win->legend(["log M > 10.8","10.2 < log M < 10.8","logM < 10.2"],-3.1,0.5,{Colour=>['Red','Green','Blue','Orange'],Symbol=>[17,16,18],Charsize=>1.1,SymbolSize=>3,TextShift=>-1.65});

$win->close();



$opt = {Device=> "$ENV{HOME}/Software/paper/figures/3micron_vs_z.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->env(0.5,2.,29.,33.5,{Axis=>'LogY',XTitle=>'z',YTitle=>'L\d3\gmm\u (Greybody/PAH) (W/\A)',AspectRatio=>1});

$Lsun = 3.826e26; #in Watts
#write flux in solar units
#$rest3micronflux /= $Lsun; #in Lsun
#$rest3micronflux_err /= $Lsun;#in Lsun
$dlogrest3micronflux_err = 0.4343 * $rest3micronflux_err/ $rest3micronflux;

#sptype cut
$idx1 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == 0 &  ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);
$idx2 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == 1 &  ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);
$idx3 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == 2 &  ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);
$idx4 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == -1 &  ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);

$win->errb($z($idx1), log10($rest3micronflux($idx1)), undef, $dlogrest3micronflux_err($idx1), {Symbol=>18,Colour=>'blue',SymbolSize=>2});
$win->errb($z($idx2), log10($rest3micronflux($idx2)), undef, $dlogrest3micronflux_err($idx2), {Symbol=>18,Colour=>'green',SymbolSize=>2});
$win->errb($z($idx3), log10($rest3micronflux($idx3)), undef, $dlogrest3micronflux_err($idx3), {Symbol=>18,Colour=>'red',SymbolSize=>2});
$win->errb($z($idx4), log10($rest3micronflux($idx4)), undef, $dlogrest3micronflux_err($idx4), {Symbol=>18,Colour=>'orange',SymbolSize=>2});

$UL3micronflux = log10($rest3micronflux_err) + $dlogrest3micronflux_err;

#sptype cut
$idx1 = which($K < 20.6  & $rest3micronflux < 0.5*$rest3micronflux_err & $spflag == 0 &  ($Conf > 1 | $z<9) & $Conf < 7);
$idx2 = which($K < 20.6  & $rest3micronflux < 0.5*$rest3micronflux_err & $spflag == 1 &  ($Conf > 1 | $z<9) & $Conf < 7);
$idx3 = which($K < 20.6  & $rest3micronflux < 0.5*$rest3micronflux_err & $spflag == 2 &  ($Conf > 1 | $z<9) & $Conf < 7);
$idx4 = which($K < 20.6  & $rest3micronflux < 0.5*$rest3micronflux_err & $spflag == -1 &  ($Conf > 1 | $z<9) & $Conf < 7);

$win->errb($z($idx1), $UL3micronflux($idx1), undef,undef, {Symbol=>31,Colour=>'blue',SymbolSize=>2});
$win->errb($z($idx2), $UL3micronflux($idx2), undef,undef, {Symbol=>31,Colour=>'green',SymbolSize=>2});
$win->errb($z($idx3), $UL3micronflux($idx3), undef,undef, {Symbol=>31,Colour=>'red',SymbolSize=>2});
$win->errb($z($idx4), $UL3micronflux($idx4), undef,undef, {Symbol=>31,Colour=>'orange',SymbolSize=>2});

$win->legend(["Red Evolved","Intermediate","Star Forming","Mixed"],0.8,33.15,{Colour=>['Red','Green','Blue','orange'],Symbol=>[17,16,18,17],Charsize=>1.1,SymbolSize=>1.5,TextShift=>-1.65});

$win->close();


$opt = {Device=> "$ENV{HOME}/Software/paper/figures/3micron_vs_mass.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->env(9.6,11.7,29.,33.5,{Axis=>'LogY',XTitle=>'log M',YTitle=>'L\d3\gmm\u (Greybody/PAH) (W/\A)',AspectRatio=>1});

$Lsun = 3.826e26; #in Watts
#write flux in solar units
#$rest3micronflux /= $Lsun; #in Lsun
#$rest3micronflux_err /= $Lsun;#in Lsun
$dlogrest3micronflux_err = 0.4343 * $rest3micronflux_err/ $rest3micronflux;

#sptype cut
$idx1 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == 0 &  ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);
$idx2 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == 1 &  ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);
$idx3 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == 2 &  ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);
$idx4 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == -1 &  ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);

$win->errb($mass($idx1), log10($rest3micronflux($idx1)), undef, $dlogrest3micronflux_err($idx1), {Symbol=>18,Colour=>'blue',SymbolSize=>2});
$win->errb($mass($idx2), log10($rest3micronflux($idx2)), undef, $dlogrest3micronflux_err($idx2), {Symbol=>18,Colour=>'green',SymbolSize=>2});
$win->errb($mass($idx3), log10($rest3micronflux($idx3)), undef, $dlogrest3micronflux_err($idx3), {Symbol=>18,Colour=>'red',SymbolSize=>2});
$win->errb($mass($idx4), log10($rest3micronflux($idx4)), undef, $dlogrest3micronflux_err($idx4), {Symbol=>18,Colour=>'orange',SymbolSize=>2});

$UL3micronflux = log10($rest3micronflux_err) + $dlogrest3micronflux_err;

#sptype cut
$idx1 = which($K < 20.6  & $rest3micronflux < 0.5*$rest3micronflux_err & $spflag == 0 &  ($Conf > 1 | $z<9) & $Conf < 7);
$idx2 = which($K < 20.6  & $rest3micronflux < 0.5*$rest3micronflux_err & $spflag == 1 &  ($Conf > 1 | $z<9) & $Conf < 7);
$idx3 = which($K < 20.6  & $rest3micronflux < 0.5*$rest3micronflux_err & $spflag == 2 &  ($Conf > 1 | $z<9) & $Conf < 7);
$idx4 = which($K < 20.6  & $rest3micronflux < 0.5*$rest3micronflux_err & $spflag == -1 &  ($Conf > 1 | $z<9) & $Conf < 7);

$win->errb($mass($idx1), $UL3micronflux($idx1), undef,undef, {Symbol=>31,Colour=>'blue',SymbolSize=>2});
$win->errb($mass($idx2), $UL3micronflux($idx2), undef,undef, {Symbol=>31,Colour=>'green',SymbolSize=>2});
$win->errb($mass($idx3), $UL3micronflux($idx3), undef,undef, {Symbol=>31,Colour=>'red',SymbolSize=>2});
$win->errb($mass($idx4), $UL3micronflux($idx4), undef,undef, {Symbol=>31,Colour=>'orange',SymbolSize=>2});

$win->legend(["Red Evolved","Intermediate","Star Forming","Mixed"],0.8,33.15,{Colour=>['Red','Green','Blue','orange'],Symbol=>[17,16,18,17],Charsize=>1.1,SymbolSize=>1.5,TextShift=>-1.65});

$win->close();



$opt = {Device=> "$ENV{HOME}/Software/paper/figures/3micron_divmass_vs_sfr.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->env(-1.5,2.,18.8,22.1,{Axis=>'LogXY',XTitle=>'SFR (M\d\(2281)\u/yr)',YTitle=>'L\d3\gmm\u (Greybody/PAH) (W/\A)',AspectRatio=>1});

$Lsun = 3.826e26; #in Watts
#write flux in solar units
#$rest3micronflux /= $Lsun; #in Lsun
#$rest3micronflux_err /= $Lsun;#in Lsun
$dlogrest3micronflux_err = 0.4343 * $rest3micronflux_err/ $rest3micronflux;

#sptype cut
$idx1 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == 0 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);
$idx2 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == 1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);
$idx3 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == 2 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);
$idx4 = which($K < 20.6  & $rest3micronflux > 0.5*$rest3micronflux_err & $spflag == -1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0);

$win->errb(log10($sfrOIIsed($idx1)), log10($rest3micronflux($idx1))-$mass($idx1), undef, $dlogrest3micronflux_err($idx1), {Symbol=>18,Colour=>'blue',SymbolSize=>2});
$win->errb(log10($sfrOIIsed($idx2)), log10($rest3micronflux($idx2))-$mass($idx2), undef, $dlogrest3micronflux_err($idx2), {Symbol=>18,Colour=>'green',SymbolSize=>2});
$win->errb(log10($sfrOIIsed($idx3)), log10($rest3micronflux($idx3))-$mass($idx3), undef, $dlogrest3micronflux_err($idx3), {Symbol=>18,Colour=>'red',SymbolSize=>2});
$win->errb(log10($sfrOIIsed($idx4)), log10($rest3micronflux($idx4))-$mass($idx4), undef, $dlogrest3micronflux_err($idx4), {Symbol=>18,Colour=>'orange',SymbolSize=>2});

$UL3micronflux = log10($rest3micronflux_err) + $dlogrest3micronflux_err;

$win->errb(log10($sfr2000sed($idx1)), $UL3micronflux($idx1),undef,undef,{Symbol=>31,Colour=>'blue',SymbolSize=>2 & $blendflag==0});
$win->errb(log10($sfr2000sed($idx2)), $UL3micronflux($idx2),undef,undef,{Symbol=>31,Colour=>'green',SymbolSize=>2 & $blendflag==0});
$win->errb(log10($sfr2000sed($idx3)), $UL3micronflux($idx3),undef,undef,{Symbol=>31,Colour=>'red',SymbolSize=>2 & $blendflag==0});
$win->errb(log10($sfr2000sed($idx4)), $UL3micronflux($idx4),undef,undef,{Symbol=>31,Colour=>'orange',SymbolSize=>2  & $blendflag==0});

#plot fit

$sfr = 10**(pdl(-4,-3,-2,1,0,1,2,3)); 
$r3mu = 2.62e31 * $sfr**(0.351);

$win->line(log10($sfr),log10($r3mu));

#use order of magnitude approximation Lexcess in Lband = 0.642 * (0.15) * sfr * lifetime
$Lexcess = 350; #calculated in circumstellar/disksize.pl
$flux3micron = $Lsun / 3e4 * $sfr * 1e6 * $Lexcess; #in units of W/A
$win->line(log10($sfr),log10($flux3micron),{Color=>green,Linestyle=>dashed});
$win->legend(["Red Evolved","Intermediate","Star Forming","Mixed"],-1.1,33.15,{Colour=>['Red','Green','Blue','orange'],Symbol=>[17,16,18,17],Charsize=>1.1,SymbolSize=>1.5,TextShift=>-1.65});

$win->close();
