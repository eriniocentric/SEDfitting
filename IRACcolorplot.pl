
use lib "$ENV{HOME}/Software";
use PDL;
use PDL::NiceSlice;
use KGB::Cosmology;
use KGB::PegUtils;
use KGB::SpecUtils;
use KGB::Dust;
use PDL::Graphics::PGPLOT;
use PDL::GSL::RNG;
use PGPLOT;
use PDL::Graphics::PGPLOTOptions ('set_pgplot_options');
set_pgplot_options('CharSize' => 2.3,'HardCH'=> 1.7,'HardLW'=>3,'AspectRatio'=>1,'TightLabels'=>1);

$flux0 = 3631; # AB zero-point in Jy

($id,$ch1,$ch2,$ch3,$ch3err,$ch4,$ch4err,$mips,$f1,$f2,$f3,$f4) = rcols("/Users/mentuch/Software/Catalogs/GDDSMIPS_summary_fr3.txt",0,29,32,35,37,38,40,41,47,48,49,50);

$cols = sequence(40)+50;


($spflag,$spectralclass,$z,$massKarl,$massKarl_err, $Conf,$weight,$sfr2000,$sfrOII,$restUB,$gini,$assym,$fac,$facburst,$burstratio,$burstratio_avg,$burstratio_err,$bbtemp,$temp_avg,$temp_err,$massVIzK,$massVIzKerr,$zmaxVIzK,$massVIzK12,$massVIzK12err,$zmaxVIzK12,$massVIzK1234,$massVIzK1234err,$zmaxVIzK1234,$massVIzK1234_cce,$massVIzK1234_cceerr,$zmaxVIzK1234_cce,$rest4micronflux,$rest4micronflux_err,$massVIzK_nb,$massVIzK_nberr,$massVIzK12_nb,$massVIzK12_nberr,$massVIzK12_nb,$massVIzK12_nberr) = rcols("$ENV{HOME}/Software/Catalogs/GDDSMIPS_summary_fr3.txt",50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90);

$flux1 = $flux0*10**(-0.4*$ch1);
$flux2 = $flux0*10**(-0.4*$ch2);
$flux3 = $flux0*10**(-0.4*$ch3);
$flux4 = $flux0*10**(-0.4*$ch4);

$ch3to1 = $flux3 / $flux1;
$ch4to2 = $flux4 / $flux2;

$mipsto4 = (1e-6) * $mips / $flux4;
$mipsUL = (1e-6) *  70. * ones(nelem($flux4)) / $flux4;
dev "$ENV{HOME}/Software/paper/figures/IRACMIPScolorplot.ps/vcps",{AspectRatio=>1};

env 0,5.1,0,35,{Col=>'Black',XRange=>[0,5],YRange=>[0,50],YTitle=>'S\d24\u/S\d8.0\u',XTitle=>'S\d8.0\u/S\d4.5\u'};  

#identify IR exess objects
#$ix = which($f4 > 1.3 & $mips > 70 & $ch4err < 1);
#points $ch4to2($ix),$mipsto4($ix),{Color=>'black',Symbol=>'circle',SymbolSize=>2.3};

#$ix = which($f4 > 1.3 & $ch4err < 1 & $mips < 70);
#points $ch4to2($ix), $mipsUL($ix), {Symbol=>'circle',color=>'black',SymbolSize=>2.3};


$ix = which($spflag == -1 & $ch4err < 0.3 & $mips > 70);
points $ch4to2($ix),$mipsto4($ix),{Symbolsize=>2,Symbol=>17,Col=>'orange'};

$ix = which($spflag == 0 & $ch4err < 0.3 & $mips > 70);
points $ch4to2($ix),$mipsto4($ix),{Symbolsize=>2,Symbol=>18,Col=>'blue'};

$ix = which($spflag == 1 & $ch4err < 0.3 & $mips > 70);
points $ch4to2($ix),$mipsto4($ix),{Symbolsize=>2,Symbol=>16,Col=>'green'};

$ix = which($spflag == 2 & $ch4err < 0.3 & $mips > 70);
points $ch4to2($ix),$mipsto4($ix),{Symbolsize=>2,Symbol=>17,Col=>'red'};


#plot upperlimits
$ix = which($ch4err < 0.3 & $mips < 70);
$mipsUL = (1e-6) *  70. * ones(nelem($flux4)) / $flux4;
points $ch4to2($ix), $mipsUL($ix)-0.6, {Symbol=>'31',SymbolSize=>1.5};


$ix = which($spflag == -1 & $ch4err < 0.3 & $mips < 70);
points $ch4to2($ix),$mipsUL($ix),{Symbolsize=>2,Symbol=>17,Col=>'orange'};
$ix = which($spflag == 0 & $ch4err < 0.3 & $mips < 70);
points $ch4to2($ix),$mipsUL($ix),{Symbolsize=>2,Symbol=>18,Col=>'blue'};
$ix = which($spflag == 1 & $ch4err < 0.3 & $mips < 70);
points $ch4to2($ix),$mipsUL($ix),{Symbolsize=>2,Symbol=>16,Col=>'green'};
$ix = which($spflag == 2 & $ch4err < 0.3 & $mips < 70);
points $ch4to2($ix),$mipsUL($ix),{Symbolsize=>2,Symbol=>17,Col=>'red'};

#$ix = which($f4 > 1.7 & $ch4err < 1 & $mips < 70);
#points $ch4to2($ix), $mipsUL($ix), {Symbol=>'circle',color=>'black',SymbolSize=>2.3};
legend(["Red Evolved","Intermediate","Star Forming","Mixed","SB DOGs","AGN DOGs","Normal DOGs","ULIRGs"],3.1,32,{Colour=>['Red','Green','Blue','orange','Black','Black','Black','Black'],Symbol=>[17,16,18,17,'cross','diamond','circle','triangle'],Charsize=>1.,SymbolSize=>1.5,TextShift=>-0.35});
#legend(["MIPS non-detection","IRAC "],1.,1,{Charsize=>0,Symbol=>['31','circle'],Color=>['black','black'],TextShift=>-0.3});

#overplot SBSMGs range from Pope et al. (2008)
($x,$y) = rcols("$ENV{HOME}/Software/Catalogs/PopePlots/SFDOGs.txt");
points $x,$y,{Colour=>'black',Symbol=>'cross'};
#overplot AGN SMGs from Pope et al 
($x,$y) = rcols("$ENV{HOME}/Software/Catalogs/PopePlots/AGNDOGs.txt");
points $x,$y,{Colour=>'black',Symbol=>'diamond'};
#overplot all other DOGs from Pope et al 
($x,$y) = rcols("$ENV{HOME}/Software/Catalogs/PopePlots/DOGs.txt");
points $x,$y,{Colour=>'black',Symbol=>'circle',Symbolsize=>1.5};
#overplot ULIRGS at z>1 from Sajina et al 2007
($x,$y) = rcols("$ENV{HOME}/Software/Catalogs/PopePlots/ULIRGs.txt");
points $x,$y,{Colour=>'black',Symbol=>'triangle'};

line pdl(2,2),pdl(0,50),{Linestyle=>'dashed'};

#now plot M82 SED evol track

($z,$S8,$S24) = rcols("$ENV{HOME}/Software/Catalogs/PopePlots/SBtrack.dat");
line $S8,$S24,{Color=>'black',Linestyle=>3,Color=>13};

$i = which( $z == 1);
text $z(($i)),$S8(($i)),$S24(($i))-1,{Charsize=>0.9,color=>13};
$i = which( $z == 2);
text $z(($i)),$S8(($i))-0.02,$S24(($i))+1,{Charsize=>0.9,color=>13};
$i = which( $z == 3);
text $z(($i)),$S8(($i))+.04,$S24(($i))-0.5,{Charsize=>0.9,color=>13};
$i = which( $z == 4);
text "M82",1.4,3,{Charsize=>0.9,color=>13};


#now plot MK231 AGN SED evol track

($z,$S8,$S24) = rcols("$ENV{HOME}/Software/Catalogs/PopePlots/AGNtrack.dat");
line $S8,$S24,{Color=>'black',Linestyle=>3,Color=>purple};

$i = which( $z == 1);
text $z(($i)),$S8(($i))+0.1,$S24(($i))+0.5,{Charsize=>0.9,Color=>purple};
$i = which( $z == 2);
text $z(($i)),$S8(($i)),$S24(($i))+.5,{Charsize=>0.9,color=>purple};
$i = which( $z == 3);
text $z(($i)),$S8(($i))+0.3,$S24(($i))+.7,{Charsize=>0.9,color=>purple};
$i = which( $z == 4);
text "Mrk 231",4.1,2.5,{Charsize=>0.9,color=>purple};

#plot stellar evol tracks

($z,$i,$j,$h,$k,$ch1,$ch2,$ch3,$ch4,$mips) = rcols("$ENV{HOME}/Software/Catalogs/sedEVOLBLUE_AV=1.dat");

#line 10**(-0.4*($ch4-$ch2)),10**(-0.4*($mips-$ch4)),{Color=>blue,Linestyle=>dotted};

($z,$i,$j,$h,$k,$ch1,$ch2,$ch3,$ch4,$mips) = rcols("$ENV{HOME}/Software/Catalogs/sedEVOLRED_AV=0.0.dat");
#line 10**(-0.4*($ch4-$ch2)),10**(-0.4*($mips-$ch4)),{Color=>red,Linestyle=>5};

close_window();


($id,$z,$Conf,$V,$Verr,$R,$Rerr,$i,$ierr,$K,$Kerr,$ch1,$ch1err,$ch2,$ch2err,$ch3,$ch3err,$ch4,$ch4err,$spflag,$f1,$f2,$f3,$f4,$blendflag) = rcols("/Users/mentuch/Software/Catalogs/GDDS_fr2.txt",0,1,5,11,13,14,16,17,19,26,28,29,31,32,34,35,37,38,40,41,45,46,47,48,49);


#convert mags to AB
$V += -0.0243;
$i += 0.4138;

$flux1 = $flux0*10**(-0.4*$ch1);
$flux2 = $flux0*10**(-0.4*$ch2);
$flux3 = $flux0*10**(-0.4*$ch3);
$flux4 = $flux0*10**(-0.4*$ch4);

$ch3to1 = $flux3 / $flux1;
$ch4to2 = $flux4 / $flux2;

dev "$ENV{HOME}/Software/paper/figures/IRAC-opticalcolorplot_zcuts.ps/vcps";

$opt_ch1_colour = $i - $ch1;
env 17.7,24,-0.1,4.5,{Color=>'black',XTitle=>'[3.6 \gmm]',YTitle=>'I - [3.6 \gmm]'};
#$ix = which($f4 > 1.3 & $ch4err < 1 & $ch3err < 1 & $z>0.5 & $z<2);
#points $ch1($ix),$opt_ch1_colour($ix),{YRange=>[-1.3,5.0],XRange=>[17,25],Color=>'black',SymbolSize=>2.3,Symbol=>'circle'};
#hold;
$ix = which($ch4err < 1 & $ch3err < 1 & $K < 20.6 & $Conf > 1  & $Conf < 7 & $blendflag == 1 );

$ix = which($ch4err < 1 & $ch3err < 1 & $z < 0.7 & $K < 20.6 & ($Conf > 1 | $z<9) & $z>0.01 & $Conf < 7 & $Conf > 1 & $blendflag == 0);
points $ch1($ix),$opt_ch1_colour($ix),{Colour=>'blue',Symbol=>18,SymbolSize=>1.6};
$ix = which($ch4err < 1 & $ch3err < 1 & $z > 0.7 & $z < 1.0 & $K < 20.6 & ($Conf > 1 | $z<9) & $z>0.01 & $Conf < 7 & $blendflag == 0);
points $ch1($ix),$opt_ch1_colour($ix),{Symbol=>16,SymbolSize=>1.6,Colour=>'cyan'};
$ix = which($ch4err < 1 & $ch3err < 1 & $z > 1.0 & $z < 1.5  & $K < 20.6 & ($Conf > 1 | $z<9) & $z>0.01 & $Conf < 7 & $blendflag == 0);
points $ch1($ix),$opt_ch1_colour($ix),{Symbol=>13,SymbolSize=>1.6,Colour=>'orange'};
$ix = which($ch4err < 1 & $ch3err < 1 & $z > 1.5 & $z < 2  & $K < 20.6 & ($Conf > 1 | $z<9) & $z>0.01 & $Conf < 7 & $blendflag == 0);
points $ch1($ix),$opt_ch1_colour($ix),{Colour=>'purple',Symbolsize=>1.6,Symbol=>17};

($ze,$mass,$ie,$je,$he,$ke,$ch1e,$ch2e,$ch3e,$ch4e) = rcols("$ENV{HOME}/Software/Catalogs/isochronesblue.dat");

$iplot = which($ze == 0.7 & $mass > 10.2);
#line $ch1e($iplot),$ie($iplot)-$ch1e($iplot),{Color=>blue,Linestyle=>dotted};
$iplot = which($ze == 0.9 & $mass > 10.2 );
#line $ch1e($iplot),$ie($iplot)-$ch1e($iplot),{Color=>green,Linestyle=>dotted};
$iplot = which($ze == 1.25 & $mass > 10.2 );
#line $ch1e($iplot),$ie($iplot)-$ch1e($iplot),{Color=>orange,Linestyle=>dotted};
$iplot = which($ze == 1.75  & $mass > 10.2);
#line $ch1e($iplot),$ie($iplot)-$ch1e($iplot),{Color=>red,Linestyle=>dotted};

($ze,$mass,$ie,$je,$he,$ke,$ch1e,$ch2e,$ch3e,$ch4e) = rcols("$ENV{HOME}/Software/Catalogs/isochronesred.dat");

$iplot = which($ze == 0.7 );
line $ch1e($iplot),$ie($iplot)-$ch1e($iplot),{Color=>blue,Linestyle=>dashed};
text "z = 0.7",22.8,1.68,{Color=>blue,Charsize=>1};
$iplot = which($ze == 0.9 );
line $ch1e($iplot),$ie($iplot)-$ch1e($iplot),{Color=>cyan,Linestyle=>dashed};
text "z = 0.9",22.8,2.28,{Color=>cyan,Charsize=>1};
$iplot = which($ze == 1.25 );
line $ch1e($iplot),$ie($iplot)-$ch1e($iplot),{Color=>orange,Linestyle=>dashed};
text "z = 1.25",22.8,3.35,{Color=>orange,Charsize=>1};
$iplot = which($ze == 1.75 );
line $ch1e($iplot),$ie($iplot)-$ch1e($iplot),{Color=>purple,Linestyle=>dashed};
text "z = 1.75",22.8,4.25,{Color=>purple,Charsize=>1};

legend(["0.5 < z < 0.7","0.7 < z < 1.0","1.0 < z < 1.5","1.5 < z < 2.0"],18.5,0.7,{Colour=>['blue','cyan','Orange','purple'],Symbol=>[18,16,13,17],Charsize=>1.15,SymbolSize=>4,TextShift=>-1.7});

release;
close_window();


dev "$ENV{HOME}/Software/paper/figures/IRACcolorplot_sptypes.ps/vcps";

env -1.3,1.3,-.65,0.35,{Color=>'black',XTitle=>'[5.8 \gmm] - [8.0 \gmm]',YTitle=>'[3.6 \gmm] - [4.5 \gmm]'};
		   
$ch1ch2 = $ch1-$ch2;
$ch3ch4 = $ch3-$ch4;

$ix = which($ch4err < 0.3 & $ch3err < 0.3 & $K < 20.6 & $spflag == 0 & ($Conf > 1 | $z<9) & $z>0.01 & $Conf < 7  & $blendflag == 0);
points $ch3ch4($ix),$ch1ch2($ix),{Colour=>'blue',Symbol=>18,SymbolSize=>2};
$ix = which($ch4err < 0.3 & $ch3err < 0.3 & $K < 20.6 & $spflag == 1 & ($Conf > 1 | $z<9) & $z>0.01 & $Conf < 7 & $blendflag == 0 );
points $ch3ch4($ix),$ch1ch2($ix),{Symbol=>16,SymbolSize=>2,Colour=>'green'};
$ix = which($ch4err < 0.3 & $ch3err < 0.3 & $K < 20.6 & $spflag == 2 & ($Conf > 1 | $z<9) & $z>0.01 & $Conf < 7 & $blendflag == 0 );
points $ch3ch4($ix),$ch1ch2($ix),{Symbol=>17,SymbolSize=>2,Colour=>'red'};
$ix = which($ch4err < 0.3 & $ch3err < 0.3 & $K < 20.6 & $spflag == -1 & ($Conf > 1 | $z<9) & $z>0.01 & $Conf < 7 & $blendflag == 0 );
points $ch3ch4($ix),$ch1ch2($ix),{Colour=>'orange',Symbolsize=>2,Symbol=>17};

#plot mid-IR selection criteria from Stern et al 2005
poly pdl(-0.061,-0.061,0.937,1.328),pdl(1.031,-0.156,0.025,1.028),{Linestyle=>'Dashed',Fill=>0};
text "AGN (Stern et al. 2005)", 0,0.29,{Charsize=>0.8};
($z,$i,$j,$h,$k,$ch1,$ch2,$ch3,$ch4,$mips) = rcols("$ENV{HOME}/Software/Catalogs/sedEVOLBLUE_AV=1.dat");
line $ch3-$ch4,$ch1-$ch2,{Color=>blue,Linestyle=>dotted};

($z,$i,$j,$h,$k,$ch1,$ch2,$ch3,$ch4,$mips) = rcols("$ENV{HOME}/Software/Catalogs/sedEVOLRED_AV=0.0.dat");
line $ch3-$ch4,$ch1-$ch2,{Color=>red,Linestyle=>5};

legend(["Red Evolved","Intermediate","Star Forming","Mixed"],0.2,-.45,{Colour=>['Red','Green','Blue','Orange'],Symbol=>[17,16,18,17],Charsize=>1.1,SymbolSize=>1.5,TextFraction=>0.8});

text "z = 0.5", -0.95,-0.35,{Charsize=>0.9};
text "z = 2", -0.8,0.22,{Charsize=>0.9};
close_window();


points $ch3ch4($ix),$ch1ch2($ix),{Symbol=>16,SymbolSize=>2,Colour=>'green'};
$ix = which($ch4err < 0.3 & $ch3err < 0.3 & $K < 20.6 & $spflag == 2 & ($Conf > 1 | $z<9) & $z>0.01 & $Conf < 7 & $blendflag == 0 );
points $ch3ch4($ix),$ch1ch2($ix),{Symbol=>17,SymbolSize=>2,Colour=>'red'};
$ix = which($ch4err < 0.3 & $ch3err < 0.3 & $K < 20.6 & $spflag == -1 & ($Conf > 1 | $z<9) & $z>0.01 & $Conf < 7 & $blendflag == 0 );
points $ch3ch4($ix),$ch1ch2($ix),{Colour=>'orange',Symbolsize=>2,Symbol=>17};

#plot mid-IR selection criteria from Stern et al 2005
poly pdl(-0.061,-0.061,0.937,1.328),pdl(1.031,-0.156,0.025,1.028),{Linestyle=>'Dashed',Fill=>0};
text "AGN (Stern et al. 2005)", 0,0.29,{Charsize=>0.8};
($z,$i,$j,$h,$k,$ch1,$ch2,$ch3,$ch4,$mips) = rcols("$ENV{HOME}/Software/Catalogs/sedEVOLBLUE_AV=1.dat");
line $ch3-$ch4,$ch1-$ch2,{Color=>blue,Linestyle=>dotted};

($z,$i,$j,$h,$k,$ch1,$ch2,$ch3,$ch4,$mips) = rcols("$ENV{HOME}/Software/Catalogs/sedEVOLRED_AV=0.0.dat");
line $ch3-$ch4,$ch1-$ch2,{Color=>red,Linestyle=>5};

legend(["Red Evolved","Intermediate","Star Forming","Mixed"],0.2,-.45,{Colour=>['Red','Green','Blue','Orange'],Symbol=>[17,16,18,17],Charsize=>1.1,SymbolSize=>1.5,TextFraction=>0.8});

text "z = 0.5", -0.95,-0.35,{Charsize=>0.9};
text "z = 2", -0.8,0.22,{Charsize=>0.9};
close_window();
