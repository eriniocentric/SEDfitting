## this plots the fluxratios between observed and model fluxes blueshifted to the rest frame
## this file uses the models calculated from GDDSVIzK12 

use PDL;
use lib "$ENV{HOME}/Software";

use PDL::NiceSlice;
use PDL::Graphics::PGPLOT;
use PGPLOT;
use PDL::Graphics::PGPLOTOptions ('set_pgplot_options');
set_pgplot_options('CharSize' => 1,'HardCH'=> 1.5,'HardLW'=>3,'AspectRatio'=>1);

use KGB::Cosmology;
use KGB::PegUtils;
use KGB::SpecUtils;
use KGB::Dust;

$ENV{DATADIR} = "$ENV{HOME}/Software/filters";

$OmegaM = 0.3; $H0 = 70;
$soften = 0.03;
$Law = "SMC";
$nebfac = 2.0; # Extra factor of extinction for neb lines
$pc = 3.086E16; # One parsec in meters
$flux0 = 3631; # AB zero-point in Jy
$clight = 2.99792e8; # speed of light in  meters
$planck = 6.6260755e-34; #planck's constant in m^2 kg /s
$boltzmann = 1.3806503e-23; # m^2 kg s-2 K
$Pi = 3.1415;
# Work out filter effective wavelengths for plots
@filters = (
#    'sys_mosaic_U',
#    'sys_mosaic_B',
    'sys_mosaic_V',
#    'sys_mosaic_R',
    'sys_mosaic_I',
    'sys_mosaic_Z',
#   'sys_circi_J',
#    'sys_circi_H',
    'sys_circi_K',
    'ch1',
    'ch2',
    'ch3',
    'ch4',
    );
@vega2ab = (
#    0.5479, #U
#   -0.1219, #B
    -0.0245, #V
#    0.1799, #R
    0.4137, #I
    0.4965, #Z
#    0.8806, #J
#    1.3281, #H
    1.8106, #K
    0,
    0,
    0,
    0,
   );

@v=(); foreach my $f (@filters) {   my ($w, $t) = get_filter($f . ".dat"); push @v, sum($w*$t)/sum($t); }
$wav0 = pdl(@v);
print "Filter eff. wavelengths = ", $wav0,"\n";

open IN, "$ENV{HOME}/Software/Catalogs/GDDSphotometry_all.txt" or die "Input file not found\n";
open OUTMAGS, ">$ENV{HOME}/Software/SINGS/restframemags_VIzK12.txt" or die "Can't open output file\n";


while(<IN>) {
   next if /^#/;

   # Read  catalog format
   @v = split;
   $id = $v[0]; $zsp = $v[1]; $zph = $v[2]; $Conf = $v[5]; 
   
   $v = pdl(@v);
   $z = $zsp;
   $weight = $v[7];
   $flag = ( ($Conf <=1 | $zsp>9) & $zsp>0.01  ); # Good photo-z only (and not star)
   $K = $v[26];
   $z=$zsp;
   if ($flag == 1) {
       $z = $zph;
   }							      
      
   next if $z<0.5 or $z>2;
#   print "######################### Object $id  ##### z = $zsp #########################\n";

   $dmags = $v->dice([11,17,20,26,29,32,35,38]) + pdl(@vega2ab); $dmags_v = $v->dice([13,19,22,28,31,34,37,40])**2 + $soften**2; #
								 # Data
								 # and
								 # variance
								 # (softened)
   $dflux = $flux0*10**(-0.4*$dmags); $dflux_v = $dmags_v * ($dflux/1.08574)**2;  # Data and errors in Jy

   $K = $dmags(3);
   #nan flux with bad errors
   $ix = which($dmags_v > 0.3); 
   
   next if($dmags_v((6,7)) > 0.3);
   $dflux($ix) .= nan;

   open BLEND, "$ENV{HOME}/Software/Catalogs/GDDS_blend.cat";
   while (<BLEND>) {
	 next unless m/^$id/;
	 @v4 = split;
	 $blendflag = $v4[1];
       }
   close BLEND;
   next if ($blendflag == 1);
open IN2, "$ENV{HOME}/Software/Masses/GDDSVIzK12/SEDparameters.dat"; #using VIzK1234 + greybody + pah

   while (<IN2>) {
       
       next unless m/^$id/;
       @v2 = split;
       $id2 = $v2[0]; $zsp2 = $v2[1]; $specfile = $v2[2]; $itime = $v2[3]; $fac = $v2[4]; $AVV =$v2[5]; 
       
       open IN4, "$ENV{HOME}/Software/Catalogs/GDDSSummary_051219.txt" or die "Input file not found\n";
       while (<IN4>) {
	   next unless m/^$id/;
	   @v3 = split;
	   $spectralclass = $v3[42]; $agn = $v3[29];
       }
       close IN4;
       print "Getting SED from $specfile\n";
       ($t2,$dummy,$wav,$spec,$emspec) = read_peg_spec($specfile,{SPLIT_EMISSION=>1}); # Note #spec returned in W/A
       

       $spec *= 1E10; $emspec *= 1E10; # Convert to 1E10 mass galaxy to avoid rounding
       
       $spec2 = $spec->dice_axis(1,$itime); 
       $emspec2 = $emspec->dice_axis(1,$itime); 
       
# Redden the spectrum in the same way expand_massgrid2.p does
       $attn = peidust($Law, $wav)/peidust($Law, 5500);
       
       $spec2 = $spec2 * 10**( -0.4 * $attn * $AVV ) + $emspec2 * 10**(-0.4 * $attn * $nebfac * $AVV);

       
	open FILE2, "$ENV{HOME}/Software/Masses/GDDSVIzK12-masses.dat" or die "Died on FILE2";
	while(<FILE2>) {
	    next unless m/$id/;

	    @v2 = split;
	    $mass = $v2[3]; $masserr = $v2[4];

	}
       close FILE2;
 
       
#extract parameters from spectrum file name
       chomp $specfile;
       print "Processing file $specfile    ";
	 
       $specfile =~ /massgrid5.Z(.*)_T(.*)_tb(.*)_logf(.*).pegspec/ ;
       $metals = $1; $Tau = $2; $tburst = $3; $fburst = $4;
       $fburst = "-9999" if $fburst eq "-none";
       
       print "Z=$metals T=$Tau tb=$tburst f=$fburst\n";  
       
#get appropriate age
       # Time steps in Myr to use for long SFH calculation (24 elements)
       $time1 =  pdl qw/ 0 25 50 75 100 200 300 400 500 600 700 800 900
1000 1200 1400 1600 1800 2000 2500.000 3000 4000 5000 7000 9000 11000 13000/;
# Extra time steps to use near vicinity of burst 	       
       $time2 = pdl qw/-25 0 25 50 75 100 200 300 400 500 600 700 800 900 1000/;

       $time = uniq(qsort(append($time1, $time2+$tburst)));
       $t = $time(($itime));
 
       $spectotal = $fac*$spec2;
  
       #compute rest frame mags
       $ABflag=1;
      
       $U = mag($wav,$spectotal,'sys_mosaic_U')->at(0);
       $B = mag($wav,$spectotal,'sys_mosaic_B')->at(0);
       $V = mag($wav,$spectotal,'sys_mosaic_V')->at(0);
       $R = mag($wav,$spectotal,'sys_mosaic_R')->at(0);
       $I = mag($wav,$spectotal,'sys_mosaic_I')->at(0);
       $z = mag($wav,$spectotal,'sys_mosaic_z')->at(0);
       $J = mag($wav,$spectotal,'2mass_J')->at(0);
       $H = mag($wav,$spectotal,'2mass_H')->at(0);
       $K = mag($wav,$spectotal,'2mass_K')->at(0);
       $L = mag($wav,$spectotal,'ch1')->at(0);
       $M = mag($wav,$spectotal,'ch2')->at(0);
       $N = mag($wav,$spectotal,'ch3')->at(0);
       $O = mag($wav,$spectotal,'ch4')->at(0);

#retrieve SFR from steph-sfr.dat if available
       $sfr2000 = -99.99;
       $sfrOII  = -99.99;
       
       open SFR, "$ENV{HOME}/Software/Catalogs/steph-sfr.dat" or die "Input file not found\n";
       while (<SFR>) {
	   next unless m/^$id/; 
	   @v3 = split;
	   $id3 = $v3[0]; $zsp3 = $v3[1]; $sfr2000 = $v3[6]; $sfrOII = $v3[7];
       }   
       close SFR;
       $attnOII = peidust($Law, 3627)/peidust($Law, 5500)->at(0);

       $sfrOII /= 10**(-0.4 * $attnOII * $nebfac * $AVV);
       $sfr2000 /= 10**(-0.4 * $attnOII * $nebfac * $AVV);
       print "$sfrOII,$sfr2000 \n";

       $spflag = -1; $col = 'black'; $sym = '17';
       
       if ($spectralclass =~ m/001/) {$spflag = 2; $col = 'red'; $sym = '17';}
       if ($spectralclass =~ m/010/) {$spflag = 1; $col = 'green'; $sym = '16';}
       if ($spectralclass =~ m/100/) {$spflag = 0; $col = 'blue'; $sym = '18'};

       printf OUTMAGS "%10s %2.0g  %5.3g  %5.3g %5.3g %5.3g   %8.6g  %8.6g  %8.6g %8.6g   %8.6g %8.6g %8.6g %8.6g  %8.6g %8.6g  %8.6g   %8.6g %8.6g  %8.6g  %8.6g %8.6g  %8.6g\n", $id,$spflag,$Conf,$K,$sfr2000,$sfrOII,$AVV,$U,$B,$V,$R,$I,$z,$J,$H,$K,$L,$M,$N,$O;
         
   }
   close IN2;
}
close OUT;
close OUTFLUX;
close OUTMAGS;
