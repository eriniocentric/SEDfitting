#adapted program to handle large images by breaking them into panels of 50x50 pixels

#uses long wavelength images remaped to 10"/pix platescale

$ENV{DATADIR} = "$ENV{HOME}/Research/filters";

use lib "$ENV{HOME}/Software";
use PDL;
use PDL::NiceSlice;
use KGB::PegUtils;
use KGB::SpecUtils;
use KGB::Dust;
use PDL::Graphics::PGPLOT;
use PDL::GSL::RNG;
use PDL::Graphics::LUT;

$PDL::BIGPDL = 1;
$PEGASE = 1;
$plotBB = 1; #flag to plot modified BB results

#### this can be set up to be automated later
#$name = 'ngc4125';
#$dist = 23.9;#Mpc from Tonry et al 2001
#$AMW = pdl(0.082, 0.063, 0.063, 0.051, 0.037, 0.017, 0.011, 0.007, 0.003); #NGC 4125 from NED

$plsc = 10.;#plate scale for FWHM=6 images

$name = 'm51';
$dist = 9.9;
$AMW = pdl(0.152, 0.117, 0.094, 0.094, 0.068, 0.032, 0.020, 0.013, 0.005); #BVRRIJHKL band foreground galactic extinction for M51

$extinctioncorrection = 10**(-0.4 * ($AMW));#divid datacube by this to correction for galactic foreground extinction


$workdir = "$ENV{HOME}/Research/VNGS/".$name."/SEDfitting";
#make SED fitting directory if it doesn't exist
system "mkdir $workdir" unless (-e $workdir);

srand(42);
$rng = PDL::GSL::RNG->new('taus');

@filters = ('kpnoB',
	    'kpnoV',
	    'kpnoHA1563',
	    'kpnoR',
	    'kpnoI',
	    '2mass_J',
	    '2mass_H',
	    '2mass_K',
	    'ch1',
	    'ch2',
	    'ch3',
	    'ch4',
	    'mips24',
	    'pacs_blue',
	    'pacs_red',
	    'spire_1ext',
	    'spire_2ext',
	    'spire_3ext'
	   );

@v=(); foreach my $f (@filters) {   my ($w, $t) = get_filter($f . ".dat"); push @v, sum($w*$t)/sum($t); }
$wavfilter = pdl(@v);
print "Filter eff. wavelengths = ", $wavfilter,"\n";
$nphot = nelem($wavfilter);

#$starphot = pdl(0,1,3,4,5,6,7); #select V through K to fit stellar models to photometry, but not Ha or B (ngc 4125 seems to have a B excess -Xray emission?
$starphot = pdl(0,1,3,4,5,6,7);
$dustphot = pdl(11,12,13,14,15,16); #select 8 through 350 micron to fit dust models to photometry

$Msolar = 1.98892e33; #mass of the sun in g
$Pi = 3.1415;
$clight = 2.99792e8; # speed of light in  meters
$planck = 6.6260755e-34; #planck's constant in m^2 kg /s
$boltzmann = 1.3806503e-23; # m^2 kg s-2 K
$sigma = 5.67e-8; # W/m^2/K^4

$H0 = 70;
$OmegaM = 0.3;
$nebfac = 1.0/0.44; # Extra factor of extinction for neb lines

$massH = 1.67e-24; #mass of hydrogen in grams
$Lsolar = 3.82e26; #luminosity of the sun in W
$massSun = 1.98e33; #mass of Sun in g
$massHsolar = $massH / $massSun; #mass of hydrogen atom in solar masses

$ABflag = 1;
$flux0 = 3631; # AB zero-point in Jy
$pc = 3.086E16; # One parsec in meters



$extendedcorrection = ones($nphot);

#to convert to flux densities from point source to extended, remain in RSRF flux weighted intentensities 
#need to divide by K4p term

$extendedcorrection(-3:-1) .= pdl(1.0113,1.0087,1.0065); #should be dividing this
#$extendedcorrection(-3:-1) .=  pdl(0.9939,0.9898,0.9773);#old one mistake!


#open up data to fit 

$imB = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_B_500.fits",{hdrcpy=>1};
$imV = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_V_500.fits",{hdrcpy=>1};
$imHa = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_Ha_500.fits",{hdrcpy=>1};
$imR = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_R_500.fits",{hdrcpy=>1};
$imI = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_I_500.fits",{hdrcpy=>1};
$imJ = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_J_500.fits",{hdrcpy=>1};
$imH = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_H_500.fits",{hdrcpy=>1};
$imK = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_K_500.fits",{hdrcpy=>1};
$im3pt6 = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_3pt6_500.fits",{hdrcpy=>1};
$im4pt5 = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_4pt5_500.fits",{hdrcpy=>1};
$im5pt6 = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_5pt6_500.fits",{hdrcpy=>1};
$im8 = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_8_500.fits",{hdrcpy=>1};
$im24 = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_24_500.fits",{hdrcpy=>1};
$im70 = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_70_500.fits",{hdrcpy=>1};
$im160 = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_160_500.fits",{hdrcpy=>1};
$im250 = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_250_500.fits",{hdrcpy=>1};
$im350 = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_350_500.fits",{hdrcpy=>1};
$im500 = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_500_500.fits",{hdrcpy=>1};

#temporary fix to 70 micron data as described in email from Helene Roussel
$im70 /= 1.119;
$im160 /= 1.174;

($xdims,$ydims) = dims( $imB);

#open up error images for SPIRE, all other bands will include calibration and skybackground
#uncertainties only
$imBerr = zeroes($imB);
$imVerr = zeroes($imV);
$imHaerr = zeroes($imHa);
$imRerr = zeroes($imR);
$imIerr = zeroes($imI);
$imJerr = zeroes($imJ);
$imHerr = zeroes($imH);
$imKerr = zeroes($imK);
$im3pt6err = zeroes($im3pt6);
$im4pt5err = zeroes($im4pt5);
$im5pt6err = zeroes($im5pt6);
$im8err = zeroes($im8);
$im24err = zeroes($im24);
$im70err = zeroes($im70);
$im160err = zeroes($im160);
if (-e "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_70err_500.fits") {
  $im70err = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_70err_500.fits",{hdrcpy=>1}; 
  $im160err = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_160err_500.fits",{hdrcpy=>1};
  $im70err /= 1.119;
  $im160err /= 1.174;
}
$im250err = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_250err_500.fits",{hdrcpy=>1};
$im350err = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_350err_500.fits",{hdrcpy=>1};
$im500err = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/".$name."_500err_500.fits",{hdrcpy=>1};


#if only want to plot up part of the image which contains the galaxy
($dum,$mean,$stddev,$med) = rcols "$ENV{HOME}/Research/VNGS/".$name."/images/processed/noise_500.dat";

#create data and dataerr cubes

$datacube = $imB->glue(2,$imV,$imHa,$imR,$imI,$imJ,$imH,$imK,$im3pt6,$im4pt5,$im5pt6,$im8,$im24,$im70,$im160,$im250,$im350,$im500); #dims are XDIMSxYDIMSxPHOTNUM
$datacube = $datacube->mv(2,0);#data cube of dims PHOTNUMXXDIMSxYDIMS
$datacube /= $extendedcorrection; #apply Herschel correction to convert to point source flux densities to extended
$datacube(0:8) /= $extinctioncorrection; #apply correction for foreground galactic extinction


$dataerr = $imBerr->glue(2,$imVerr,$imHaerr,$imRerr,$imIerr,$imJerr,$imHerr,$imKerr,$im3pt6err,$im4pt5err,$im5pt6err,$im8err,$im24err,$im70err,$im160err,$im250err,$im350err,$im500err); #dims are XDIMSxYDIMSxPHOTNUM
$dataerr = $dataerr->mv(2,0);#data cube of dims PHOTNUMxXDIMSxYDIMS
$dataerr /= $extendedcorrection; #apply Herschel correction to convert to RSRF weigth flux densities
$dataerr(0:8) /= $extinctioncorrection; #apply correction for foreground galactic extinction

#add stddev of sky background as additional error 
#Most images have photometric accuracy within 5% (broad band) or 10% (narrow band). (SINGS handbook) #they are not background subtracted!
#JHK 2MASS calibration uncertainties are 2-3% ( Jarrett et al 2003)
#calibration uncertainties in IRAC are 2-3% (IRAC Handbook)
#cabliration uncertainties in MIPS 24 is 2% (SINGS delivery pdf)
#calibration uncertainties are 10% and 20% for PACS 70,s160 um (Poglitsch 2010 Table 5) and 7% for all SPIRE bands (SPIRE OM S5.2.10)
#$calerr = pdl(0.05,0.05,0.1,0.05,0.05,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.1,0.2,0.07,0.07,0.07) * $datacube;
#updated PACS errors to match that from updated PACS handbook
$calerr = pdl(0.05,0.05,0.1,0.05,0.05,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.05,0.07,0.07,0.07) * $datacube;

$dataerr = sqrt( $dataerr**2 + $stddev**2 + $calerr**2); #uncertainties added in quadrature


$SN = $datacube > 3*$dataerr; #create NPHOTxXPIXxYPIX mask where pixels have S/N> 1

$comSN1 = sumover($SN(0:7)); #sumover S/N images, will include pixels that have good S/N in at least 3 stellar SED bands
$comSN2 = sumover($SN(10:16)); #and good S/N in three IR bands

$galaxymask = $comSN1 > 4 & $comSN2 > 2;
#$galaxymask = ones($xdims,$ydims);

#$invstarmask = rfits "$ENV{HOME}/Research/VNGS/".$name."/images/processed/stars.fits";
#$starmask = $invstarmask == 0;
#$galaxymask *= $starmask;
$datacube *= $galaxymask->dummy(0);
$dataerr *= $galaxymask->dummy(0);#multiply by galaxy mask

########OPEN STAR MODEL FLUXES AND PARAMETERS


if ($PEGASE) {
  $starsgrid = rfits "$ENV{HOME}/Research/VNGS/code/stars-massgrid-PEGASE-10Mpc.fits";
}
else {
  $starsgrid = rfits "$ENV{HOME}/Research/VNGS/code/stars-massgrid-FSPS-10Mpc.fits";
}

$modelnum = $starsgrid((0),:);
$logt = $starsgrid((1),:);
$logMs = $starsgrid((2),:);
$metals = $starsgrid((3),:);
$opticaldepth = $starsgrid((4),:);#in Myr
$Tau = $starsgrid((5),:);#in Myr
$tburst = $starsgrid((6),:);#in Myr
$logf = $starsgrid((7),:);
 
$nstarmodel = nelem $logt;

$galaxyage = 10**$logt;# in yr
$galaxyage *= 1e-6;#in Myr

if ($PEGASE) {

#$modelsel = which(($metals > 0.01 & $galaxyage > $tburst & $galaxyage > 5e3 & $opticaldepth < 1 & $logf != 0));
$modelsel = which($metals > 0.01 & $galaxyage > $tburst & ($galaxyage-$tburst) < 1e3 & $galaxyage > 2e3 & $opticaldepth < 0.9 & $logf < -0.3);
}
else {
#  $modelsel = which($metals > 0.01 & $galaxyage > $tburst & (($galaxyage-$tburst) < 1e3 ) );
 $modelsel = which($metals > 0.01 & $galaxyage > $tburst & ($galaxyage-$tburst) < 1e3 & $galaxyage > 2e3  & $opticaldepth < 1.6 & $logf < -0.3);#this is for M51
# $modelsel = which($metals > 0.01 & $galaxyage > $tburst & $galaxyage > 2e3 );
}

$starfluxes = pdl($starsgrid(8:-1,$modelsel)/($dist/10)**2);#convert to galaxy distance since models are at arbitrary distance of 10 Mpc

#######OPEN DUST MODEL FLUXES AND PARAMETERS


$dustgrid = rfits "$ENV{HOME}/Research/VNGS/code/dust-massgrid-10Mpc.fits";

$Umin = $dustgrid((0),:);
$Umax = $dustgrid((1),:);
$qpah = $dustgrid((2),:);
$ipah = $dustgrid((3),:);
$gamma = $dustgrid((4),:);

$dustfluxes = $dustgrid(5:-1,:)->copy;#dims are PHOTNUMxNDUSTMODEL and convert to galaxy distance
$dustfluxes /= ($dist/10)**2; #dims are PHOTNUMxNDUSTMODEL and convert to galaxy distance
$dustfluxes = $dustfluxes->dummy(1,$ydims);
$dustfluxes = $dustfluxes->dummy(1,$xdims); #dims are PHOTNUMxXDIMSxYDIMSxNDUSTMODEL


if (1) {#to skip fitting

system "rm -rf $workdir/results/*";
system "mkdir $workdir/results";
system "mkdir $workdir/results/images";
system "mkdir $workdir/results/images/ps";

$subpanelx = 20.;
$subpanely = 20.;

$nsubpanelx = ceil($xdims/$subpanelx);
$nsubpanely = ceil($ydims/$subpanely);

$starfluxes = $starfluxes->dummy(1,$subpanely);
$starfluxes = $starfluxes->dummy(1,$subpanelx); #dims are PHOTNUMxXDIMSxYDIMSxNSELMODEL

$Ntrials = 1;

$Uminmap = zeroes($Ntrials,$xdims,$ydims);
$Umaxmap = zeroes($Ntrials,$xdims,$ydims);
$qpahmap = zeroes($Ntrials,$xdims,$ydims);
$gammamap = zeroes($Ntrials,$xdims,$ydims);
$dustchi2map = zeroes($Ntrials,$xdims,$ydims);
$dustmap = zeroes($Ntrials,$xdims,$ydims);
$dustmodelindex = zeroes($Ntrials,$xdims,$ydims);

$metalsmap     = zeroes($Ntrials,$xdims,$ydims);
$modelnummap   = zeroes($Ntrials,$xdims,$ydims);
$opticaldepthmap         = zeroes($Ntrials,$xdims,$ydims);
$logtmap = zeroes($Ntrials,$xdims,$ydims);
$starchi2map  = zeroes($Ntrials,$xdims,$ydims);
$starmassmap   = zeroes($Ntrials,$xdims,$ydims);
$tburstmap = zeroes($Ntrials,$xdims,$ydims);
$fburstmap = zeroes($Ntrials,$xdims,$ydims);
$Taumap = zeroes($Ntrials,$xdims,$ydims);

$starfluximage = zeroes($nphot,$xdims,$ydims);
$dustfluximage = zeroes($nphot,$xdims,$ydims);

#start MC trials
$k = 0;

for ($k == 0; $k<$Ntrials; $k++) {
      
  $datacube2 = $datacube->copy;
  if ($k>0) { # 1st trial = orig. photometry, subsequent trials = M-C the fluxes + errors
    $datacube2 = $datacube + $dataerr*( grandom($datacube))/3;#tighten the gaussian distribution 
  }
  
  $xpanel = 0;
  
  for ($xpanel = 0; $xpanel < $nsubpanelx; $xpanel++) { 
    $xpanel_end = $subpanelx-1;
    $ypanel = 0;
    for ($ypanel = 0; $ypanel < $nsubpanely; $ypanel++) {
      $ypanel_end = $subpanely-1;
    
      print "On trial $k for panel $xpanel,$ypanel of $nsubpanelx,$nsubpanely panels \n";
      $xs = $xpanel * $subpanelx;
      $xe = (($xpanel+1) * $subpanelx) -1;
      $ys = $ypanel * $subpanely;
      $ye = (($ypanel+1) * $subpanely) -1;
      
      if ($xe > $xdims-1) {
	$xpanel_end = $xdims-$xs-1;
	$xe = $xdims-1;
      }
      if ($ye > $ydims-1) {
	$ypanel_end = $ydims-$ys-1;	
	$ye = $ydims-1;
      }

      $normstar = sumover($starfluxes($starphot,0:$xpanel_end,0:$ypanel_end) * $datacube2($starphot,$xs:$xe,$ys:$ye) / $dataerr($starphot,$xs:$xe,$ys:$ye)**2) / sumover($starfluxes($starphot,0:$xpanel_end,0:$ypanel_end)**2 / $dataerr($starphot,$xs:$xe,$ys:$ye)**2);
      $chi2 = sumover( ($datacube2($starphot,$xs:$xe,$ys:$ye) - $normstar->dummy(0) * $starfluxes($starphot,0:$xpanel_end,0:$ypanel_end))**2 / $dataerr($starphot,$xs:$xe,$ys:$ye)**2 );	
      $normstar = $normstar->mv(2,0);
      $chi2 = $chi2->mv(2,0);  
      
      $imin = minimum_ind($chi2);
      $modelindex = $imin;
      #die;
      #get model fluxes and renormalize to get stellar mass in each pixel/sr
      $wav=0;
      for ( $wav==0; $wav < $nphot; $wav++) {
	$mwvi= $wav+8; #model wave index, need to add 8 to account for first seven columns in starsgrid not being photometry
	$starflux = pdl($starsgrid(($mwvi)));
	
	$starfluximage(($wav),$xs:$xe,$ys:$ye) .= $starflux($modelsel)->($modelindex)/ ($dist/10)**2;
      }
      
      $starfluximage(13:-1) .= 0; #zero flux for wav>30 micron
      
      $mstar = sumover($starfluximage($starphot,$xs:$xe,$ys:$ye) * $datacube2($starphot,$xs:$xe,$ys:$ye) / $dataerr($starphot,$xs:$xe,$ys:$ye)**2) / sumover($starfluximage($starphot,$xs:$xe,$ys:$ye)**2 / $dataerr($starphot,$xs:$xe,$ys:$ye)**2); 
      badmask($mstar->inplace,0);
      
      $starfluximage(:,$xs:$xe,$ys:$ye) *= $mstar->dummy(0); #normalize fluxes
      
      $logtmap(($k),$xs:$xe,$ys:$ye) .= $logt($modelsel)->($modelindex);
      $metalsmap(($k),$xs:$xe,$ys:$ye) .= $metals($modelsel)->($modelindex);
      $opticaldepthmap(($k),$xs:$xe,$ys:$ye) .= $opticaldepth($modelsel)->($modelindex);
      $Taumap(($k),$xs:$xe,$ys:$ye) .= $Tau($modelsel)->($modelindex);
      $tburstmap(($k),$xs:$xe,$ys:$ye) .= $tburst($modelsel)->($modelindex);
      $fburstmap(($k),$xs:$xe,$ys:$ye) .= 10**($logf($modelsel)->($modelindex));
      
      $starmassmap(($k),$xs:$xe,$ys:$ye) .= $mstar*(10**($logMs($modelsel)->($modelindex)));
      $chi2sort = qsort($chi2);
      $starchi2map(($k),$xs:$xe,$ys:$ye) .= $chi2sort((0));
      $modelnummap(($k),$xs:$xe,$ys:$ye) .= $modelnum($modelsel)->($modelindex);
      
      #now fit dust to data - modelled stars
      
      $datacubedust = $datacube2-$starfluximage;
      dev '/xs';
      ctab lut_data 'rainbow';
      fits_imag $starfluximage(5),{Justify=>1};
      $normdust = sumover($dustfluxes($dustphot,$xs:$xe,$ys:$ye) * $datacubedust($dustphot,$xs:$xe,$ys:$ye) / $dataerr($dustphot,$xs:$xe,$ys:$ye)**2) / sumover($dustfluxes($dustphot,$xs:$xe,$ys:$ye)**2 / $dataerr($dustphot,$xs:$xe,$ys:$ye)**2);
      
      #calculate chisquare statistic, number of measurements N = 13 photometric bands
      
      $chi2 = sumover( ($datacubedust($dustphot,$xs:$xe,$ys:$ye) - $normdust->dummy(0) * $dustfluxes($dustphot,$xs:$xe,$ys:$ye))**2 / $dataerr($dustphot,$xs:$xe,$ys:$ye)**2 );	 
      badmask($normdust->inplace,0);			
      badmask($chi2->inplace,999999);
      $normdust = $normdust->mv(2,0);
      $chi2 = $chi2->mv(2,0);
      
      $imin = minimum_ind($chi2);
      
      #get model fluxes and renormalize to get dust mass in each pixel/sr
      $wav=0;
      
      for ( $wav==0; $wav < $nphot; $wav++) {
	$mwvi= $wav+5; #model wave index, need to add 5 to account for first five columns in modelgrid not being photometry
	$dustflux = $dustgrid(($mwvi));
	
	$dustfluximage(($wav),$xs:$xe,$ys:$ye) .= $dustflux($imin)/($dist/10)**2;
	
      }
      badmask($dustfluximage->inplace, 0);
       
      $mdust = sumover($dustfluximage($dustphot,$xs:$xe,$ys:$ye) * $datacubedust($dustphot,$xs:$xe,$ys:$ye) / $dataerr($dustphot,$xs:$xe,$ys:$ye)**2) / sumover($dustfluximage($dustphot,$xs:$xe,$ys:$ye)**2 / $dataerr($dustphot,$xs:$xe,$ys:$ye)**2); #dims are XDIMSxYDIMS
      badmask($mdust->inplace,0);
      
      $dustfluximage(:,$xs:$xe,$ys:$ye) *= $mdust->dummy(0); #normalize fluxes

      $Uminmap(($k),$xs:$xe,$ys:$ye) .= $Umin($imin);
      $Umaxmap(($k),$xs:$xe,$ys:$ye) .= $Umax($imin);
      $qpahmap(($k),$xs:$xe,$ys:$ye) .= $qpah($imin);
      $gammamap(($k),$xs:$xe,$ys:$ye) .= $gamma($imin);
      $dustmap(($k),$xs:$xe,$ys:$ye) .= $mdust;
      $chi2sort = qsort($chi2);
      $dustchi2map(($k),$xs:$xe,$ys:$ye) .= $chi2sort((0));#store best fit model chi2 values
      $dustmodelindex(($k),$xs:$xe,$ys:$ye) .= $imin; #store best fit model index
      
    }#close y panel loop
  }#close x panel loop
  
  #normalize dustflux image and starflux image
  #$starfluximage *= $starmassmap(($k))->dummy(0);
  #$dustfluximage *= $dustmap(($k))->dummy(0);
  
  #plot up trial results
  dev '/xs',2,2,{AspectRatio=>1,WindowWidth=>8};
  ctab(lut_data('rainbow'));
  fits_imag $starmassmap(($k)),1e6,{ITF=>'log',Justify=>1,Title=>"Trial $k - Stellar Mass"};
  fits_imag $logtmap(($k)),{Justify=>1,Title=>'log t\dstar\u'};
  fits_imag $Uminmap(($k)),{ITF=>'log',Justify=>1,Title=>"Umin"};
  fits_imag $dustmap(($k)),0,{ITF=>'log',Justify=>1,Title=>"Dust Mass"};

  
  if ($k == 0) {
    
    #plot up image, model image and residuals in each filter
    $totalmodelflux = zeroes($datacube);

    $totalmodelflux(0:12) .= $starfluximage(0:12); #insert model flux images for stellar models from B to 24 micron
    $totalmodelflux(5:-1) .= $totalmodelflux(5:-1) + $dustfluximage(5:-1);#add dust model for J-IR bands 
    
    $residualimage = abs($totalmodelflux-$datacube)/$dataerr;
    badmask($residualimage->inplace,0);
    $wav=0;

    for ( $wav==0; $wav < $nphot; $wav++) {
      $band = $bands[$wav];
      dev "Research/VNGS/".$name."/SEDfitting/results/images/ps/$band.ps/cps",3,1;
      ctab(lut_data('rainbow'));
      fits_imag $datacube($wav),0,{Justify=>1};
      fits_imag $totalmodelflux($wav),0,{Justify=>1};
      fits_imag $residualimage($wav),0,3,{Justify=>1};
      close_window;
      $difimage = $datacube(($wav)) - $totalmodelflux(($wav));
      wfits $difimage, "Research/VNGS/".$name."/SEDfitting/results/images/difimage_$band.fits";
      wfits $residualimage(($wav)),"Research/VNGS/".$name."/SEDfitting/results/images/residual_$band.fits";
      wfits $totalmodelflux(($wav)),"Research/VNGS/".$name."/SEDfitting/results/images/model_$band.fits";
      
    }
    
    
    
  
  #dev '/xs',2,2,{AspectRatio=>1,WindowWidth=>10};
    dev "Research/VNGS/".$name."/SEDfitting/results/exampleSED.ps/cps",2,2;
  
  #select indices to plot
    @title = ('Centre','~10" offset','~20"offset','~30"offset');
    #do something fancy here with rvals to determine pixel value at each of these offsets      
    $xplot = floor($xdims/2)+sequence(4);
   if ($name =~ /m51/) {
     $xplot = pdl(22,16,17,17);
     $yplot = pdl(40,50,29,66);
     @title = ('NGC 5194 Center','NGC 5194 Inner Arm','NGC 5194 Interarm','NGC 5195');
   }
    #$yplot = floor($ydims/2)*ones(4);
  
    $i=0;
    for ($i=0;$i<nelem($xplot);$i++) {
      
      points log10($wavfilter(0:-2)*1e-4),log10($datacube(0:-2,$xplot(($i)),$yplot(($i)))),{XTitle=>'\gl (\gmm)',YTitle=>'I\d\gn\u (MJy/sr)',Col=>black,XRange=>[-1,3],YRange=>[-2,4],Axis=>'logxy',SymbolSize=>1.6};
      hold;
      errb log10($wavfilter(0:-2)*1e-4),log10($datacube(0:-2,$xplot(($i)),$yplot(($i)))),0.4343*3*$dataerr(0:-2,$xplot(($i)),$yplot(($i)))/$datacube(0:-2,$xplot(($i)),$yplot(($i))),{Col=>black};
      $text = $title[$i];
      text $text,-0.8,3.5,{Charsize=>1.4};
      
      #plot up best fitting stars model
         
      if ($PEGASE) { 
	open STARFILE, "$ENV{HOME}/Research/PEGASE.2/modelspectra/00index" or die "Can not open 00index file\n"; # List of models
	
	$imodel = $modelnummap(0,($xplot($i)),($yplot($i)))->at(0);
	$counter = 0;
	while ($file = <STARFILE>) {
	  if ($counter == $imodel) {
	    $starfile = $file;
	    $starfile =~ s/^\s+//;#remove leading white space
	    $starfile =~ s/\s+$//;#remove trailing white space
	  }
	  $counter++;
	  
	}
	close STARFILE;
	
	($t,$Ms, $wav, $starspec,$emspec) = read_peg_spec("$ENV{HOME}/Research/PEGASE.2/modelspectra/$starfile",{SPLIT_EMISSION=>1}); # Note #spec returned in W/A/Msolar
	$logage = $logtmap(0,($xplot($i)),($yplot($i)))->at(0);
	$age = 1e-6*10**$logage;
	
	$selt = qsorti( abs($t - $age))->at(0);
	
	$A_v = $opticaldepthmap(0,($xplot($i)),($yplot($i)))->at(0);    
	
	# Redden the spectrum 
	$attn = peidust('MW', $wav)/peidust('MW', 5500);
	
	$spec2 = $starspec * 10**(-0.4 * ($attn * $A_v)) + $emspec * 10**(-0.4 * $attn * $nebfac * $A_v);#apply factor of 2 additional extinction for nebular lines
	$spec2 = $spec2(:,($selt));
	$Mdiv = $Ms(($selt));

	$spec2 /= 4*$Pi * ($dist*1E6*$pc)**2;

	$starmass = $starmassmap(0,($xplot($i)),($yplot($i)))->at(0);

	$spec = $spec2 *($starmass/($Mdiv));
      }
      
      else {
	open STARFILE, "$ENV{HOME}/Software/FSPS/OUTPUTS/00index" or die "Can not# open 00index file\n"; # List of models
	
	$imodel = $modelnummap(0,($xplot($i)),($yplot($i)))->at(0);
	$counter = 0;
	while ($file = <STARFILE>) {
	    if ($counter == $imodel) {
	      $starfile = $file;
	      $starfile =~ s/^\s+//;#remove leading white space
	      $starfile =~ s/\s+$//;#remove trailing white space
	    }
	    $counter++;
	    
	  }
	  close STARFILE;
	
	($wav,$spec,$t,$Ms, $logLbol, $logSFR) = read_spec($starfile);#spectrum returned in W/A and wave in A
	
	$age = $logtmap(0,($xplot($i)),($yplot($i)))->at(0);
	$selt = qsorti(abs( $t - $age))->at(0);
	
	$logMdiv = $Ms(($selt));
	$spec = $spec(:,($selt));
	
	$spec /= 4*$Pi * ($dist*1E6*$pc)**2;
	
	$starmass = $starmassmap(0,($xplot($i)),($yplot($i)))->at(0);
	
	$spec *=($starmass/(10**($logMdiv)));
      }
      
      line log10($wav*1e-4),log10($spec*($wav/1e4)**2/3e-16),{Col=>purple};
      $starmass *= 1e6*($plsc/3600)**2/3282.85;
      text 'M\dstar\u = '.sprintf("%5.2g",$starmass).' M\d\(2281)\u',-0.8,3.2,{Charsize=>1.4};
      
      points log10($wavfilter($starphot)*1e-4),log10($starfluximage($starphot,$xplot(($i)),$yplot(($i)))),{SymbolSize=>3,Symbol=>'triangle',Color=>purple};
      
      #select best model fit for pixel
      $idustmodel = $dustmodelindex((0),($xplot($i)),($yplot($i)))->at(0);
      $Umini = $Umin($idustmodel)->at(0);
      $Umaxi = $Umax($idustmodel)->at(0);
      $ipahi = $ipah($idustmodel)->at(0);
      $gammai = $gamma($idustmodel)->at(0); 
      
      if ($Umini < 10) {
	$ISMdustfile = "$ENV{HOME}/Research/DustModels/DL07spec/U".sprintf("%3.2f",$Umini)."/U".sprintf("%3.2f",$Umini)."_".sprintf("%3.2f",$Umini)."_MW3.1_".sprintf("%1s0",$ipahi).".txt";
	$PDRdustfile = "$ENV{HOME}/Research/DustModels/DL07spec/U".sprintf("%3.2f",$Umini)."/U".sprintf("%3.2f",$Umini)."_1e6_MW3.1_".sprintf("%1s0",$ipahi).".txt";
      } else {
	$ISMdustfile = "$ENV{HOME}/Research/DustModels/DL07spec/U".sprintf("%3.1f",$Umini)."/U".sprintf("%3.1f",$Umini)."_".sprintf("%3.1f",$Umini)."_MW3.1_".sprintf("%1s0",$ipahi).".txt";
	$PDRdustfile = "$ENV{HOME}/Research/DustModels/DL07spec/U".sprintf("%3.1f",$Umini)."/U".sprintf("%3.1f",$Umini)."_1e6_MW3.1_".sprintf("%1s0",$ipahi).".txt";
      }

      ($wavdust,$ISMspecdust,$jnu) = rcols $ISMdustfile, {Lines=>'-1:62'};#wav in micron,flux units of ergs/s/hydrogen nucleon
      #convert to W/A / 1E10 Mdust(Msolar)
      $wavdust *= 1e4; #lambda in A
      
      ($wavdust2,$PDRspecdust,$jnu) = rcols $PDRdustfile, {Lines=>'-1:62'};#wav in micron,flux units of ergs/s/hydrogen nucleon
      #convert to W/A / 1E10 Mdust(Msolar)

      $mass = $dustmap(0,($xplot($i)),($yplot($i)))->at(0);

      $ISMspecdust /= 0.01*$massHsolar; #erg/s/Mdust
      $ISMspecdust *= 1e-7 / $wavdust; # in W/A/Mdust (Msolar)
      $ISMspecdust /= 4*$Pi * ($dist*1E6*$pc)**2; #W/m^2/A     
      $ISMspecdust *= $mass;
      
      $PDRspecdust /= 0.01*$massHsolar; #erg/s/Mdust
      $PDRspecdust *= 1e-7 / $wavdust; # in W/A/Mdust (Msolar)
      $PDRspecdust /= 4*$Pi * ($dist*1E6*$pc)**2; #W/m^2/A     
      $PDRspecdust *= $mass;
      
      $specdust =  (1-$gammai) * $ISMspecdust + $gammai * $PDRspecdust; #units of ergs/s/Hnucleon

      line log10($wavdust*1e-4),log10($specdust* ($wavdust/1e4)**2 / 3e-16),{Col=>red};

      line log10($wavdust*1e-4),log10((1-$gammi)*$ISMspecdust* ($wavdust/1e4)**2 / 3e-16),{LineStyle=>dashed,Col=>orange};
      line log10($wavdust*1e-4),log10($gammai*$PDRspecdust* ($wavdust/1e4)**2 / 3e-16),{LineStyle=>dotted,Col=>orange};

      $mass *= 1e6*($plsc/3600)**2/3282.85;
      text 'M\ddust\u = '.sprintf("%5.2g",$mass).' M\d\(2281)\u',-0.8,2.9,{Charsize=>1.4};
      points log10($wavfilter($dustphot)*1e-4),log10($dustfluximage($dustphot,$xplot(($i)),$yplot(($i)))),{SymbolSize=>3,Symbol=>'triangle',Color=>red};
      
      #also plot modified BB result
      if ($plotBB) {
	$betamean = rfits "Research/VNGS/".$name."/SEDfitting/results-IR/beta.fits";
	$betarms = rfits "Research/VNGS/".$name."/SEDfitting/results-IR/betaerr_MC.fits";
	$tempmean = rfits "Research/VNGS/".$name."/SEDfitting/results-IR/temp.fits";
	$temprms = rfits "Research/VNGS/".$name."/SEDfitting/results-IR/temperr_MC.fits";
	$uncert = rfits "Research/VNGS/".$name."/SEDfitting/results-IR/temperr_chi2.fits";
	$massdustmean = rfits "Research/VNGS/".$name."/SEDfitting/results-IR/dustmass.fits";
	$massdustrms = rfits "Research/VNGS/".$name."/SEDfitting/results-IR/dustmasserr.fits";
	$chi2mean = rfits "Research/VNGS/".$name."/SEDfitting/results-IR/chi2.fits";
	
	$photsel = pdl(13,14,15,16);#fitted BB for 70 to 500um obs
	$bbtemp = $tempmean(($xplot($i)),($yplot($i)));
	$bbtemperr = $temprms(($xplot($i)),($yplot($i)));
	$beta = $betamean(($xplot($i)),($yplot($i)));
	$betaerr  = $betarms(($xplot($i)),($yplot($i)));
	$dflux = $datacube(:,($xplot($i)),($yplot($i)));
	$dfluxerr = $dataerr(:,($xplot($i)),($yplot($i)));
	$mdust = $massdustmean(($xplot($i)),($yplot($i)));
	$mdusterr = $massdustrms(($xplot($i)),($yplot($i)));
	
	$mdust *= 1e6*($plsc/3600)**2/3282.85;
	$mdusterr *= 1e6*($plsc/3600)**2/3282.85;
	#plot some pixels 
  
	$wav2 = 1e-8 * sequence(100000) + 1e-7; #wav in meters 
	$blackbody = 2 * $planck * ($clight**2) * $wav2**(-5) * ( exp( $clight * $planck /($wav2 * $boltzmann * $bbtemp))-1)**(-1); #in W/m^2/sr/m
	$dw = $wav2(1:-1)- $wav2(0:-2);
	
	$integrate = sumover($dw*$blackbody(1:-1));#should add upto $sigma*T^4/$Pi (and it does)
	
	$blackbody /= 1e10; #in W/m^2/A/sr      
	
	$greybody = 1e-10*$blackbody/(($wav2)**($beta));
	
	$wav2 *= 1e10; #convert to angstroms
	
	
	@bbflux = ();
	@gbflux = ();
	foreach $filt (@filters) {
    
	  $greymag = mag($wav2,$greybody,$filt);
	  $bbmag = mag($wav2,$blackbody,$filt);
	  print $greymag,$bbmag;
	  $fluxgrey = $flux0 * 10**(-0.4*$greymag); #Jy
	  $fluxbb = $flux0 * 10**(-0.4*$bbmag); #Jy
	  #      $fluxfilt = flux($wav2,$flux,$filt); #W/A/Msolar (emission per steradian per solar mass
	  
	  badmask($fluxgrey, -99999, $fluxgrey);
	  badmask($fluxbb, -99999, $fluxbb);
	  $nnn = sum($fluxbb < -999);
	  
	  push @bbflux, $fluxbb;
	  push @gbflux, $fluxgrey;
	} # Loop over filter
	$modelfluxes = pdl(@gbflux)->flat;
	$bbfluxes = pdl(@bbflux)->flat;    
	
	$norm = sumover($modelfluxes($photsel) * $dflux($photsel) / $dfluxerr($photsel)**2) /sumover($modelfluxes($photsel)**2 / $dfluxerr($photsel)**2);
	
	$chi2 = sumover( ($dflux($photsel) - $norm->dummy(0) * $modelfluxes($photsel))**2 / $dfluxerr($photsel)**2 ); 					 
#	errb log10($wavfilter*1e-4),log10($dflux),0.4343*sqrt($dfluxerr**2)/$dflux,{Col=>red};
#	$text = $title[$i];
#	text $text,1.8,1.5;
	text "T=".sprintf("%3.1f",$bbtemp).'+/-'.sprintf("%3.1f",$bbtemperr)."K",1.7,0.5,{Charsize=>1.2};
	text 'M\ddust\u='.sprintf("%5.2g",$mdust).' M\d\(2281)\u',1.7,0.1,{Charsize=>1.2};
	line log10($wav2* 1e-4),log10($norm*$greybody*(($wav2/1e4)**2)/3e-16),{LineStyle=>dashed,col=>green};
	points log10($wavfilter($photsel)*1e-4),log10($norm*$modelfluxes($photsel)),{Symbolsize=>1.6,Symbol=>'triangle',col=>green};
      }

      points log10($wavfilter(0:-2)*1e-4),log10($datacube(0:-2,$xplot(($i)),$yplot(($i)))),{XTitle=>'\gl (\gmm)',YTitle=>'I\d\gn\u (MJy/sr)',Col=>black,XRange=>[-1,3],YRange=>[-2,4],Axis=>'logxy',SymbolSize=>2};
 release;
      
    }
    close_window;
  }#end of k=0 uncertainty and plotting bit
  
}#end of trials  
($logtmean,$logtrms) = statsover($logtmap);
($starmassmean,$starmassrms) = statsover($starmassmap);
($tburstmean,$tburstrms) = statsover($tburstmap);
($fburstmean,$fburstrms) = statsover($fburstmap);
($Taumean,$Taurms) = statsover($Taumap);
($metalsmean,$metalsrms) = statsover($metalsmap);
($opticaldepthmean,$opticaldepthrms) = statsover($opticaldepthmap);
($Uminmean,$Uminrms) = statsover($Uminmap);
($Umaxmean,$Umaxrms) = statsover($Umaxmap);
($qpahmean,$qpahrms) = statsover($qpahmap);
($gammamean,$gammarms) = statsover($gammamap);
($dustmassmean,$dustmassrms) = statsover($dustmap);
($dustchi2mean,$dustchi2rms) = statsover($dustchi2map);
($starchi2mean,$starchi2rms) = statsover($starchi2map);
$timeburstmap = 10**$logtmap - $tburstmap*1e6;
($timeburstmean,$timeburstrms) = statsover($timeburstmap);

$logtmean    *= $galaxymask;
$logtrms     *= $galaxymask;
$metalsmean     *= $galaxymask;
$metalsrms      *= $galaxymask;
$opticaldepthmean         *= $galaxymask;
$opticaldepthrms          *= $galaxymask;
$tburstmean     *= $galaxymask;
$tburstrms      *= $galaxymask;
$fburstmean     *= $galaxymask;
$fburstrms      *= $galaxymask;
$Taumean        *= $galaxymask;
$Taurms         *= $galaxymask;
$starmassmean   *= $galaxymask;
$starmassrms    *= $galaxymask;
$Uminmean       *= $galaxymask;
$Uminrms        *= $galaxymask;
$qpahmean       *= $galaxymask;
$qpahrms        *= $galaxymask;
$gammamean      *= $galaxymask;
$gammarms       *= $galaxymask;
$dustmassmean   *= $galaxymask;
$dustmassrms    *= $galaxymask;

wfits $logtmean,"Research/VNGS/".$name."/SEDfitting/results/t_star_mean.fits";
wfits $logtrms,"Research/VNGS/".$name."/SEDfitting/results/t_star_rms.fits";
wfits $metalsmean,"Research/VNGS/".$name."/SEDfitting/results/Z_mean.fits";
wfits $metalsrms,"Research/VNGS/".$name."/SEDfitting/results/Z_rms.fits";
wfits $opticaldepthmean,"Research/VNGS/".$name."/SEDfitting/results/opticaldepth_mean.fits";
wfits $opticaldepthrms,"Research/VNGS/".$name."/SEDfitting/results/opticaldepth_rms.fits";
wfits $tburstmean,"Research/VNGS/".$name."/SEDfitting/results/tburst_mean.fits";
wfits $tburstrms,"Research/VNGS/".$name."/SEDfitting/results/tburst_rms.fits";

wfits $timeburstmean,"Research/VNGS/".$name."/SEDfitting/results/timeburst_mean.fits";
wfits $timeburstrms,"Research/VNGS/".$name."/SEDfitting/results/timeburst_rms.fits";

wfits $fburstmean,"Research/VNGS/".$name."/SEDfitting/results/fburst_mean.fits";
wfits $fburstrms,"Research/VNGS/".$name."/SEDfitting/results/fburst_rms.fits";
wfits $Taumean,"Research/VNGS/".$name."/SEDfitting/results/Tau_mean.fits";
wfits $Taurms,"Research/VNGS/".$name."/SEDfitting/results/Tau_rms.fits";
wfits $starmassmean,"Research/VNGS/".$name."/SEDfitting/results/Mstar_mean.fits";
wfits $starmassrms,"Research/VNGS/".$name."/SEDfitting/results/Mstar_rms.fits";
wfits $Uminmean,"Research/VNGS/".$name."/SEDfitting/results/Umin_mean.fits";
wfits $Uminrms,"Research/VNGS/".$name."/SEDfitting/results/Umin_rms.fits";
wfits $qpahmean,"Research/VNGS/".$name."/SEDfitting/results/qpah_mean.fits";
wfits $qpahrms,"Research/VNGS/".$name."/SEDfitting/results/qpah_rms.fits";
wfits $gammamean,"Research/VNGS/".$name."/SEDfitting/results/gamma_mean.fits";
wfits $gammarms,"Research/VNGS/".$name."/SEDfitting/results/gamma_rms.fits";
wfits $dustmassmean,"Research/VNGS/".$name."/SEDfitting/results/Mdust_mean.fits";
wfits $dustmassrms,"Research/VNGS/".$name."/SEDfitting/results/Mdust_rms.fits";
wfits $dustchi2mean,"Research/VNGS/".$name."/SEDfitting/results/chi2_dust.fits";
wfits $starchi2mean,"Research/VNGS/".$name."/SEDfitting/results/chi2_star.fits";

wfits $galaxymask,"Research/VNGS/".$name."/SEDfitting/results/galaxymask.fits";

}#close if statement to skip fitting

$logtmean    = rfits "Research/VNGS/".$name."/SEDfitting/results/t_star_mean.fits",{hdrcpy=>1};
$logtrms     = rfits "Research/VNGS/".$name."/SEDfitting/results/t_star_rms.fits",{hdrcpy=>1};
$metalsmean     = rfits "Research/VNGS/".$name."/SEDfitting/results/Z_mean.fits",{hdrcpy=>1};
$metalsrms      = rfits "Research/VNGS/".$name."/SEDfitting/results/Z_rms.fits",{hdrcpy=>1};
$opticaldepthmean         = rfits "Research/VNGS/".$name."/SEDfitting/results/opticaldepth_mean.fits",{hdrcpy=>1};
$opticaldepthrms          = rfits "Research/VNGS/".$name."/SEDfitting/results/opticaldepth_rms.fits",{hdrcpy=>1};
$tburstmean     = rfits "Research/VNGS/".$name."/SEDfitting/results/tburst_mean.fits",{hdrcpy=>1};
$tburstrms      = rfits "Research/VNGS/".$name."/SEDfitting/results/tburst_rms.fits",{hdrcpy=>1};
$timeburstmean     = rfits "Research/VNGS/".$name."/SEDfitting/results/timeburst_mean.fits",{hdrcpy=>1};
$timeburstrms      = rfits "Research/VNGS/".$name."/SEDfitting/results/timeburst_rms.fits",{hdrcpy=>1};
$fburstmean     = rfits "Research/VNGS/".$name."/SEDfitting/results/fburst_mean.fits",{hdrcpy=>1};
$fburstrms      = rfits "Research/VNGS/".$name."/SEDfitting/results/fburst_rms.fits",{hdrcpy=>1};
$Taumean        = rfits "Research/VNGS/".$name."/SEDfitting/results/Tau_mean.fits",{hdrcpy=>1};
$Taurms         = rfits "Research/VNGS/".$name."/SEDfitting/results/Tau_rms.fits",{hdrcpy=>1};
$starmassmean   = rfits "Research/VNGS/".$name."/SEDfitting/results/Mstar_mean.fits",{hdrcpy=>1};
$starmassrms    = rfits "Research/VNGS/".$name."/SEDfitting/results/Mstar_rms.fits",{hdrcpy=>1};
$Uminmean       = rfits "Research/VNGS/".$name."/SEDfitting/results/Umin_mean.fits",{hdrcpy=>1};
$Uminrms        = rfits "Research/VNGS/".$name."/SEDfitting/results/Umin_rms.fits",{hdrcpy=>1};
$dustmassmean   = rfits "Research/VNGS/".$name."/SEDfitting/results/Mdust_mean.fits",{hdrcpy=>1};
$dustmassrms    = rfits "Research/VNGS/".$name."/SEDfitting/results/Mdust_rms.fits",{hdrcpy=>1};
$qpahmean       = rfits "Research/VNGS/".$name."/SEDfitting/results/qpah_mean.fits",{hdrcpy=>1};
$qpahrms        = rfits "Research/VNGS/".$name."/SEDfitting/results/qpah_rms.fits",{hdrcpy=>1};
$gammamean      = rfits "Research/VNGS/".$name."/SEDfitting/results/gamma_mean.fits",{hdrcpy=>1};
$gammarms       = rfits "Research/VNGS/".$name."/SEDfitting/results/gamma_rms.fits",{hdrcpy=>1};


$Uminres = badmask(($Uminrms/$Uminmean)->inplace,0);
$dustmassres = badmask(($dustmassrms/$dustmassmean)->inplace,0);
$starmassres = badmask(($starmassrms/$starmassmean)->inplace,0);

$dusttotalmass = sum($dustmassmean)/sum($galaxymask);
$dusttotalmasserr =sqrt( sum( ($dustmassrms**2)))/sum($galaxymask);

$startotalmass = sum($starmassmean)/sum($galaxymask);
$startotalmasserr =sqrt( sum( ($starmassrms**2)))/sum($galaxymask);

$dusttostarmass = badmask(($dustmassmean/$starmassmean)->inplace,0);

printf "Stellar Total Mass is %10.4g +/- %10.4g", $startotalmass,$startotalmasserr;
printf "Dust Total Mass is %10.4g +/- %10.4g", $dusttotalmass,$dusttotalmasserr;


#dev '/xs',2,2,{WindowWidth=>8,AspectRatio=>1};
dev "$ENV{HOME}/Research/VNGS/".$name."/SEDfitting/results/massplots.ps/vcps",2,3;
ctab(lut_data('rainbow'));
fits_imag $starmassmean,0,{ITF=>'log',Justify=>1,Title=>'Mstar (Msolar/pix)'};
fits_imag $starmassres,0,1,{Justify=>1,Title=>'Mstar residuals'};
fits_imag $dustmassmean,0,{ITF=>'log',Justify=>1,Title=>'Mdust (Msolar/pix)'};
fits_imag $dustmassres,0,1,{Justify=>1,Title=>'Mdust residuals'};
fits_imag $dusttostarmass,0,{ITF=>'sqrt',Justify=>1,Title=>'Mdust/Mstar'};     
fits_imag $opticaldepthmean,{Justify=>1,Title=>'optical depth'};
fits_imag $logtmean,{Justify=>1,Title=>'log10 Age_star'};  
fits_imag $metalsmean,{Justify=>1,Title=>'Z'};
fits_imag $timeburstmean,0,{Justify=>1,Title=>'burst age'};
fits_imag $fburstmean,{Justify=>1,Title=>'burst mass fraction'};
fits_imag $Uminmean,0,{ITF=>'sqrt',Justify=>1,Title=>'Umin'};
fits_imag $qpahmean,0,5,{Justify=>1,Title=>'QPAH'};
close_window;


# Stuff for buffer manipulation
BEGIN{
 $BUFFNROWS = 5000000; # This should be bigger than we need 
 $counter = -9999; # Ptr to where we are in the buffer
}

# Put stuff on big output array
sub wbuff {
   my @x = @_;
   my ($j, $val);
   # Init buffer
   if (ref($buff) eq "") {
       $buff = zeroes(float, scalar(@x),$BUFFNROWS); 
       $counter = 0; # Ptr to the first element free
   }
   my $n = nelem($x[0]); # Size (note 1st arg must be vector for this to work)
   my $ct2 = $counter + $n - 1;
   for ($j=0; $j<scalar(@x); $j++) {
       $buff(($j), $counter:$ct2) .= $x[$j]; # Handles scalars and vectors!!!
   }
   $counter += $n;
   die "Output buffer size exceeded" if $counter>$buff->getdim(1);
}



