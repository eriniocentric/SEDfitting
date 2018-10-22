#@id_list = ("SA02-1777","SA12-7045","SA12-8668","SA15-6695");
#$id = "SA12-6633";
#@id_list = ("SA02-1255","SA02-1842","SA12-5836","SA12-7972",,"SA12-8025","SA12-5869","SA12-6072","SA12-5592","SA12-8895","SA15-4367","SA15-7543","SA15-5005");

#@id_list = ("SA02-0744","SA02-0756","SA02-0782","SA02-1280","SA02-1400","SA02-1636","SA02-1702","SA02-1722","SA02-1741","SA02-1785","SA02-2134","SA02-2171","SA02-2530","SA12-5957","SA12-5965","SA12-6896","SA12-6974","SA12-7099","SA12-7595","SA12-7660","SA12-8139","SA12-9012","SA15-3841","SA15-3853","SA15-4231","SA15-4367","SA15-4634","SA15-4762","SA15-5731","SA15-6396","SA15-6718","SA15-6968","SA15-7501","SA15-9333");

@id_list = ("SA22-0435");

#list of objects with ACS image info
#@id_list = ("SA02-0744","SA02-0756","SA02-1280","SA02-1400","SA02-1636","SA02-1702","SA02-1722","SA02-2134","SA02-2171","SA02-2530","SA12-5957","SA12-6896","SA12-7099","SA12-7660","SA15-5731","SA15-6718","SA15-6968","SA15-7501","SA15-9333");

use lib "$ENV{HOME}/Software";
use PDL;
use PDL::NiceSlice;
use KGB::Cosmology;
use KGB::PegUtils;
use KGB::SpecUtils;
use KGB::Dust;
use PDL::Graphics::PGPLOT::Window;
use PDL::Graphics::PGPLOTOptions ('set_pgplot_options');
set_pgplot_options('CharSize' => 1.0,'HardCH'=> 1.6,'HardLW'=>2.5,AspectRatio=>1);

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

# Work out filter effective wavelengths for plots
@filters = (
#    'sys_mosaic_U',
#    'sys_mosaic_B',
    'sys_mosaic_V',
#    'sys_mosaic_R',
    'sys_mosaic_I',
    'sys_mosaic_Z',
#    'sys_circi_J',
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


foreach $id (@id_list) {
   
    open IN, "$ENV{HOME}/Software/Catalogs/GDDSphotometry_all.txt" or die "Input file not found\n";
    
#first will do it for VIzK1234 mags
    
    while(<IN>) {
	next unless m/^$id/;
	
	# Read  catalog format
	s/\*//;
	@v = split;

	$id1 = $v[0]; $zsp = $v[1]; $zph = $v[2]; $Conf = $v[5];
	
	$flag = ( ($Conf <=1 | $zsp>9) & $zsp>0.01  ); # Good photo-z only (and not star)
    								   
	$z=$zsp;
	if ($flag == 1) {
	    $z = $zph;
	}			
#	$massVIzK = $v[43]; $massVIzKerr = $v[44];$massVIzK1234 = $v[60]; $massVIzK1234err = $v[61];
 #get masses directly from output file instead     
	
	open FILE1, "$ENV{HOME}/Software/Masses/GDDSVIzK1234-masses.dat" or die "Died on FILE1";
	while(<FILE1>) {
	    next unless m/$id/;
	    @v1 = split;
	    print @v1;
	    $massVIzK1234 = $v1[3]; $massVIzK1234err = $v1[4];
	    print "massVIzK1234 is $massVIzK1234 \n";
	}
	close FILE1;
	open FILE2, "$ENV{HOME}/Software/Masses/GDDSVIzK1234_greybody_plusPAH_850K-masses.dat" or die "Died on FILE2";
	while(<FILE2>) {
	    next unless m/$id/;
	    @v2 = split;
	    print @v2;
	    $massstarburst = $v2[3]; $massstarbursterr = $v2[4];
	    print "massstarburst is $massstarburst\n";
	}
	close FILE2;
	print "Masses are $massVIzK1234 and $massstarburst \n";
	$v = pdl(@v);
	
	print "######################### Object $id1  ##### z = $z #####CONF = $Conf #######\n";
	
	$dmags = $v->dice([11,17,20,26,29,32,35,38]) + pdl(@vega2ab); $dmags_v = $v->dice([13,19,22,28,31,34,37,40])**2 + $soften**2; #
	# Data
	# and
	# variance
	# (softened)

	$dflux = $flux0*10**(-0.4*$dmags); $dflux_v = $dmags_v * ($dflux/1.08574)**2;  # Data and errors in Jy
	
   $ix = which($dmags_v<1);  # This keeps all VIzK points, and throws out only a couple of really bad spitzer points
    }
    close IN;

    open IN2, "$ENV{HOME}/Software/Masses/GDDSVIzK1234_greybody_plusPAH_850K/SEDparameters.dat";
    while (<IN2>) {
	next unless m/^$id/;
	@v = split;
	$id2 = $v[0]; $zsp = $v[1]; $specfile = $v[2]; $itime = $v[3]; $fac = $v[4]; $facburst = $v[5]; $AVV =$v[6];
    }
    close IN2;

    print "Getting SED from $specfile\n";
    ($t2,$dummy, $wav,$spec,$emspec) = read_peg_spec($specfile,{SPLIT_EMISSION=>1}); # Note #spec returned in W/A
    $spec *= 1E10; $emspec *= 1E10; # Convert to 1E10 mass galaxy to avoid rounding
    
    $spec2 = $spec->dice_axis(1,$itime); 
    $emspec2 = $emspec->dice_axis(1,$itime); 
    
# Redden the spectrum in the same way expand_massgrid2.p does
    $attn = peidust($Law, $wav)/peidust($Law, 5500);
    
    $spec2 = $spec2 * 10**( -0.4 * $attn * $AVV ) + $emspec2 * 10**(-0.4 * $attn * $nebfac * $AVV);
    
#    $fileburst = "$ENV{HOME}/Software/starburst/SBgrid/1kpc/OB40/L1050Av35.4_1d4.plot";
  
#    print "Opening starburst file: $fileburst\n";
    
#    ($wav2,$specburst2) = rcols($fileburst,{Lines=>"-1:3"}); #wav in microns, spec in Jy
    
#   $wav2 *= 1e4; #convert to angstroms
#    $specburst2 *= 1e-26*($clight / ($wav2 * 1e-10)) / ($wav2); #convert from jy to W/A
#    $specburst2 *= 4*$Pi * (50. *1E6*$pc)**2; #write in W/m^2/A 
    
    #now fit the additional starburst component. Must interpolate burst spectrum  to match PEGASE SED

    $wav2 = 1e-9 * sequence(20000) + 1e-7; #wav in meters 
    $specburst2 = 2 * $planck * ($clight ** 2) * $wav2 ** (-6) * ( exp( $clight * $planck /($wav2 * $boltzmann * 850))) ** (-1); #in W/m^3

    $specburst2 *= 1e10; #in W/m^2/A 

    $wav2 *= 1e10; #convert to angstroms


# add 3.3 um PAH feature following Lorentian Profile parametrized by Verstraete et al. (2001)

   $PAH = 2.27 * $specburst2->index(which( $wav2 == 40000)) / (1 + ((1 / ($wav2 * 1e-8) - 3039.1)/19.4)**2); 

    $specburst2 += $PAH;

    $specburst = interpol($wav,$wav2,$specburst2);
    
    $spectotal = $fac*$spec2 + $facburst*$specburst;
    
    ($wz, $sz) = redshift_spectra($wav, $spectotal, pdl($z));
    
# Coerce 1D
    die "Something wrong \n" if (nelem($sz) ne $sz->getdim(0)) or (nelem($wz) ne $wz->getdim(0)); # Should never happen! 
    $wz = $wz->clump(-1)->copy; $sz = $sz->clump(-1)->copy;
    
    $DL = lumdist($z);
    
    $sz /= 4*$Pi * ($DL*1E6*$pc)**2; 
   
# $sz now matches that used to produce the colors for this model
     
# Output fit spectra to file
#$win = PDL::Graphics::PGPLOT::Window->new(Device=>"$ENV{HOME}/Software/Figures/$id-SED_VIzK1234.ps/cps"); 
    $win = PDL::Graphics::PGPLOT::Window->new(Device=>"$ENV{HOME}/Software/IRexcess/Figures/cce/$id-SEDcombined.ps/cps"); 
#    $win = PDL::Graphics::PGPLOT::Window->new(Device=>'/xs');
    $plotlogflux = log10($dflux($ix));
    $dplotlogflux = 0.4343 * sqrt($dflux_v($ix)) / $dflux($ix);
    $win->env(log10(3000),log10(100000),max($plotlogflux)-2.5,max($plotlogflux)+0.7,{Axis=>'LogXY',XTitle=>'Wavelength (\A)',YTitle=>'Flux (Jy)'});
    $win->points(log10($wav0($ix)), $plotlogflux,{Col=>Red,SymbolSize=>1.5});
    $win->errb(log10($wav0($ix)), $plotlogflux, $dplotlogflux,{Col=>Red});
    $win->line(log10($wz), log10( $sz * ($wz/10000)**2 / 3e-16 )); # Plotting in Jy
    $win->text('log\d10\uM = '.$massstarburst.'\(2233)'.$massstarbursterr." Greybody + PAH", log10(3300),max($plotlogflux)-2.3,{Charsize=>1.2});
    $win->text('log\d10\uM = '.$massVIzK1234.'\(2233)'.$massVIzK1234err.' Stellar only', log10(3400),max($plotlogflux)-2.1,{Charsize=>1.2});
    $win->text("$id: z = $zsp",log10(3400),max($plotlogflux)-1.9,{Charsize=>1.2});
#    $win->close();
    
#now we can do it for PEGASE fit without cce component
    
    open IN3, "$ENV{HOME}/Software/Masses/GDDSVIzK1234/SEDparameters.dat";
    
    while (<IN3>) { 
	next unless m/^$id/;
	@v = split;
	$id = $v[0]; $zsp = $v[1]; $specfile = $v[2]; $itime = $v[3]; $fac = $v[4]; $AVV =$v[5];
    }
    
    print "Getting SED from $specfile\n";
    ($t2,$dummy, $wav,$spec,$emspec) = read_peg_spec($specfile,{SPLIT_EMISSION=>1}); # Note #spec returned in W/A
    $spec *= 1E10; $emspec *= 1E10; # Convert to 1E10 mass galaxy to avoid rounding
    
    $spec2 = $spec->dice_axis(1,$itime); 
    $emspec2 = $emspec->dice_axis(1,$itime); 
    
# Redden the spectrum in the same way expand_massgrid2.p does
    $attn = peidust($Law, $wav)/peidust($Law, 5500);
    
    $spec2 = $spec2 * 10**( -0.4 * $attn * $AVV ) + $emspec2 * 10**(-0.4 * $attn * $nebfac * $AVV);
    
    ($wz, $sz) = redshift_spectra($wav, $spec2, pdl($z));
    
# Coerce 1D
    die "Something wrong \n" if (nelem($sz) ne $sz->getdim(0)) or (nelem($wz) ne $wz->getdim(0)); # Should never happen! 
    $wz = $wz->clump(-1)->copy; $sz = $sz->clump(-1)->copy;

    $DL = lumdist($z);
    
    $sz /= 4*$Pi * ($DL*1E6*$pc)**2; 
# $sz now matches that used to produce the colors for this model
    $sz *= $fac; # Normalization
    
    $win->line(log10($wz), log10( $sz * ($wz/10000)**2 / 3e-16 ),{LineStyle=>'Dashed'});
#    $win->legend(["VIz'K + IRAC + Starburst SED fit","VIz'K + IRAC SED fit"],log10(10000),max($plotlogflux)-2.3,{LineStyle=>['Solid','Dashed'],TextFraction=>0.75});
    $win->close();

    $winfits = PDL::Graphics::PGPLOT::Window->new(Device=>"$ENV{HOME}/Software/IRexcess/Figures/$id-fits.ps/cps",AspectRatio=>1); 
    $fitsimage = rfits( "$ENV{HOME}/Software/IRexcess/fits_images/".$id."_F814W.fits");
    #$win->env(0,1,0,1,{PlotPosition => [0.7,0.9,0.2,0.4]});
    $winfits->fits_imag($fitsimage,{DrawWedge=>0,Charsize=>1.3});
    $winfits->close();
    
# Output fit spectra to file
#    $win1 = PDL::Graphics::PGPLOT::Window->new(Device=>"$ENV{HOME}/Software/Figures/$id-SED_VIzK.ps/cps"); 
#    $plotlogflux = log10($dflux($ix));
#    $dplotlogflux = 0.4343 * sqrt($dflux_v($ix)) / $dflux($ix);
#    $win1->env(log10(3000),log10(100000),max($plotlogflux)-2.5,max($plotlogflux)+0.7,{Axis=>'LogXY',XTitle=>'Wavelength (\A)',YTitle=>'Flux (Jy)',Title=>"$id"});
#    $win1->points(log10($wav0($ix)), $plotlogflux,{Col=>Red});
#    $win1->errb(log10($wav0($ix)), $plotlogflux, $dplotlogflux,{Col=>Red});
#    $win1->line(log10($wz), log10( $sz * ($wz/10000)**2 / 3e-16 )); # Plotting in Jy
#    $win1->close();
}
