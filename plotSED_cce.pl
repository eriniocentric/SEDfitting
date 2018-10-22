#@id_list = ("SA02-1777","SA12-7045","SA12-8668","SA15-6695");
#$id = "SA12-6633";
#@id_list = ("SA02-1255","SA02-1842","SA12-5836","SA12-7972",,"SA12-8025","SA12-5869","SA12-6072","SA12-5592","SA12-8895","SA15-4367","SA15-7543","SA15-5005");

#@id_list = ("SA02-0744","SA02-0756","SA02-0782","SA02-1280","SA02-1400","SA02-1636","SA02-1702","SA02-1722","SA02-1741","SA02-1785","SA02-2134","SA02-2171","SA02-2530","SA12-5957","SA12-5965","SA12-6896","SA12-6974","SA12-7099","SA12-7595","SA12-7660","SA12-8139","SA12-9012","SA15-3841","SA15-3853","SA15-4231","SA15-4367","SA15-4634","SA15-4762","SA15-5731","SA15-6396","SA15-6718","SA15-6968","SA15-7501","SA15-9333");


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
set_pgplot_options('CharSize' => 1.0,'HardCH'=> 1.8,'HardLW'=>2.5,AspectRatio=>0.618);

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
#open to store mass to light ratio
$ML = zeroes(229);
$iml =0;

#foreach $id (@id_list) {

open IN, "$ENV{HOME}/Software/Catalogs/GDDSphotometry_all.txt" or die "Input file not found\n";

#first will do it for VIzK12 mags

while(<IN>) {
    next if /^#/;
	
	# Read  catalog format
	
	@v = split;
    
    $id = $v[0]; $zsp = $v[1]; $zph = $v[2]; $Conf = $v[5];
   
    $flag = ( ($Conf <=1 | $zsp>9) & $zsp>0.01  ); # Good photo-z only (and not star)
    
    $z=$zsp;
    if ($flag == 1) {
	$z = $zph;
    }		
    next if $z < 0.5 or $z > 2;
#	$massVIzK = $v[43]; $massVIzKerr = $v[44];$massVIzK1234 = $v[60]; $massVIzK1234err = $v[61];
    #get masses directly from output file instead     
    
    open FILE1, "$ENV{HOME}/Software/Masses/GDDSVIzK12-masses.dat" or die "Died on FILE1";
    while(<FILE1>) {
	next unless m/$id/;
	@v1 = split;
	print @v1;
	$massVIzK12 = $v1[3]; $massVIzK12err = $v1[4];
	print "massVIzK12 is $massVIzK12 \n";
    }
    close FILE1;
    open FILE2, "$ENV{HOME}/Software/Masses/GDDSVIzK1234_GC-masses.dat" or die "Died on FILE2";
    while(<FILE2>) {
	next unless m/$id/;
	@v2 = split;
	print @v2;
	$massstarburst = $v2[3]; $massstarbursterr = $v2[4];
	$facburst_avg = $v2[5]; $facburst_err = $v2[6]; 
	$fac_avg = $v2[28]; $fac_err = $v2[29];
	print "massstarburst is $massstarburst\n";
    }
	close FILE2;

       $burstratio_avg = ($facburst_avg*$specburst->index( which( $wav == 30050))) / ($fac_avg * $spec2->index( which ($wav == 30050)));
       $burstratio_avg = $burstratio_avg->at(0);
       $burstratio_err = 0;

       
       if ($burstratio_avg > 0) {
	   $burstratio_err = $burstratio_avg * sqrt(($facburst_err/$facburst_avg)**2 + ($fac_err/$fac_avg)**2);
       }
    
    print "Masses are $massVIzK12 and $massstarburst \n";
    $v = pdl(@v);
    
    print "######################### Object $id  ##### z = $z #####CONF = $Conf #######\n";
    
    $dmags = $v->dice([11,17,20,26,29,32,35,38]) + pdl(@vega2ab); $dmags_v = $v->dice([13,19,22,28,31,34,37,40])**2 + $soften**2; #
	# Data
    # and
    # variance
	# (softened)
    
    $dflux = $flux0*10**(-0.4*$dmags); $dflux_v = $dmags_v * ($dflux/1.08574)**2;  # Data and errors in Jy
    $dflux_v(4:7) *= 2; # multiply spizter errors by 2
    $ix = which($dmags_v<0.3);  # This keeps all VIzK points, and throws out only a couple of really bad spitzer points
   
    next if ($dmags_v((7)) > 0.3);
    
   open BLEND, "$ENV{HOME}/Software/Catalogs/GDDS_blend.cat";
   while (<BLEND>) {
	 next unless m/^$id/;
	 @v4 = split;
	 $blendflag = $v4[1];
       }
    close BLEND;
    next if ($blendflag == 1);

    open IN2, "$ENV{HOME}/Software/Masses/GDDSVIzK1234_GC/SEDparameters.dat";
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
    

#get greybody mags in flux space, fit the additional starburst component

    ($wav2,$specburst2) = rcols("$ENV{HOME}/Software/IRexcess/fit_ism_sed.dat",0,1); #wav in microns, starts at 1.03 microns, flux in lambda*F_lambda
    $wav2 *= 1e4; #convert to angstroms
    $specburst2 *=1E45; #round up to match W/A of PEGASE models
    
    $specburst2 /= $wav2; #spec now in F_lambda
    
    $specburst = interpol($wav,$wav2,$specburst2);
    
    $spectotal = $fac*$spec2 + $facburst*$specburst;
    
    $absR = mag($wav,($spectotal / ( 4*$Pi*(10.*$pc)**2)),'sys_mosaic_R')->at(0);

    $absK = mag($wav,($spectotal / ( 4*$Pi*(10.*$pc)**2)),'sys_circi_K')->at(0);

    ($wz, $sz) = redshift_spectra($wav, $spectotal, pdl($z));
    
    #also do for PAH modified greybody
    ($wzburst,$szburst) = redshift_spectra($wav,$facburst*$specburst,pdl($z));
    #also do for stellar comp separately
    ($wzst,$szst) = redshift_spectra($wav,$fac*$spec2,pdl($z));
		
# Coerce 1D
    die "Something wrong \n" if (nelem($sz) ne $sz->getdim(0)) or (nelem($wz) ne $wz->getdim(0)); # Should never happen! 
    $wz = $wz->clump(-1)->copy; $sz = $sz->clump(-1)->copy;
    
    $DL = lumdist($z);
    
    $sz /= 4*$Pi * ($DL*1E6*$pc)**2; 

# Coerce 1D
    die "Something wrong \n" if (nelem($szburst) ne $szburst->getdim(0)) or (nelem($wzburst) ne $wzburst->getdim(0)); # Should never happen! 
    $wzburst = $wzburst->clump(-1)->copy; $szburst = $szburst->clump(-1)->copy;
    
    $DL = lumdist($z);
    
    $szburst /= 4*$Pi * ($DL*1E6*$pc)**2; 

   # Coerce 1D
    die "Something wrong \n" if (nelem($szst) ne $szst->getdim(0)) or (nelem($wzst) ne $wzst->getdim(0)); # Should never happen! 
    $wzst = $wzst->clump(-1)->copy; $szst = $szst->clump(-1)->copy;
    
    $DL = lumdist($z);
    
    $szst /= 4*$Pi * ($DL*1E6*$pc)**2; 
    
     
#calculate model mags
    $ABflag = 1;
    $count = 0;
    $mags = zeroes(8);
    foreach $filt (@filters) {
	$mfilt = mag($wz,$sz,$filt)->((0),:);  # Note unit dim is removed
	badmask($mfilt, -99999, $mfilt);
	$nnn = sum($mfilt < -999);
	$mags(($count)) .= $mfilt;
	$count++;
    }
    
 # get spectral type from GDDSsummary

    
    open IN4, "$ENV{HOME}/Software/Catalogs/GDDSSummary_051219.txt" or die "Input file not found\n";
    while (<IN4>) {
	next unless m/^$id/;
	@v3 = split;
	$spectralclass = $v3[42]; $agn = $v3[29];
	$spectralclass =~ s/\"//;
	$spectralclass =~ s/\"//;
	$spectralclass =~ s/\*//;
    }
    close IN4;
    
    $flux = $flux0 * 10**(-0.4*$mags);	 # Convert to flux space
# Output fit spectra to file
#$win = PDL::Graphics::PGPLOT::Window->new(Device=>"$ENV{HOME}/Software/Figures/$id-SED_VIzK1234.ps/cps"); 
    $win = PDL::Graphics::PGPLOT::Window->new(Device=>"$ENV{HOME}/Software/paper/figures/cce/$id-SEDcombined.ps/vcps"); 
#    $win = PDL::Graphics::PGPLOT::Window->new(Device=>'/xs');
    $plotlogflux = log10($dflux($ix));
    $dplotlogflux = 0.4343 * sqrt($dflux_v($ix)) / $dflux($ix);
    $win->env(log10(3000),log10(100000),max($plotlogflux)-2.5,max($plotlogflux)+0.7,{Axis=>'LogXY',XTitle=>'\gl (\A)',YTitle=>'F\d\gn\u (Jy)'});
    $win->points(log10($wav0($ix)), $plotlogflux,{Col=>Red,Symbol=>17,SymbolSize=>4});
    $win->errb(log10($wav0($ix)), $plotlogflux, $dplotlogflux,{Col=>red});
    $win->line(log10($wz), log10( $sz * ($wz/10000)**2 / 3e-16 )); # Plotting in Jy
#plot PAH modified greybody
    $win->line(log10($wzburst), log10( $szburst * ($wzburst/10000)**2 / 3e-16 ),{LineStyle=>'dotted',Color=>'blue'}); # Plotting in Jy
#plot stellar comp seperately
    $win->line(log10($wzst), log10( $szst * ($wzst/10000)**2 / 3e-16 ),{LineStyle=>'dashed',Color=>'green'}); # Plotting in Jy

    $obs33 = ones(2) * 3 * (1+$z) * 1e4;
    #for reference plot 3 microns
    $win->line(log10($obs33),pdl(-1,-10),{Color=>red,Linestyle=>solid});
    $win->text('\gl\drest\u=3\gmm',log10($obs33((0)))-0.23,max($plotlogflux)-2.3,{Colour=>red,Charsize=>1.2});

#    $win->points(log10($wav0($ix)),log10($flux($ix)),{Col=>blue}); #plotting model as mags
  #  $win->text('log\d10\uM = '.$massstarburst.'\(2233)'.$massstarbursterr." Stellar + Greybody", 4.1,max($plotlogflux)-2.3,{Charsize=>1.2});
  #  $win->text('log\d10\uM = '.$massVIzK12.'\(2233)'.$massVIzK12err.' Stellar only', 4.1,max($plotlogflux)-2.1,{Charsize=>1.2});
    $win->text("$id: z = $zsp",log10(3400),max($plotlogflux)+0.1,{Charsize=>1.4});
    $win->text('F\d3\gmm\u(GB/PAH)/F\d3\gmm\u(stellar) ='.sprintf( "%5.2f",$burstratio_avg)." +/-".sprintf( "%5.2f",$burstratio_err),log10(3400),max($plotlogflux)+0.35,{Charsize=>1.4});
    $win->text("SpType = $spectralclass",log10(3400),max($plotlogflux)-0.1,{Charsize=>1.4});
    $win->close();
    
#now we can do it for PEGASE fit without cce component to VIzK12 mags
#$ML($iml) .= 10**($massstarburst) * 10**(0.4*($absR - 4.76)); 
$ML($iml) .= 10**($massstarburst) * 10**(0.4*($absK - 5.14)); 

print $ML;
$iml++;
  }
