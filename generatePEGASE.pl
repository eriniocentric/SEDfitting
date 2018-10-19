#!/usr/local/bin/perl

#create SEDs for consisting of SED models from PEGASE code
#placed at arbitrary distance of 10 Mpc

use lib "$ENV{HOME}/Software";

use PDL;
use PDL::NiceSlice;
use PDL::Math;
use KGB::PegUtils;
use KGB::SpecUtils;
use KGB::Cosmology;
use KGB::Dust;

$PDL::BIGPDL = 1;

$ENV{DATADIR} = "$ENV{HOME}/Research/filters";

@filters = ('uprime',
	    'gprime',
	    'rprime',
	    'iprime',
	    'zprime',
	    'UKIRT_J',
	    'UKIRT_H',
	    'UKIRT_K',
	    'WISE-1',
	    'WISE-2',
	    'WISE-3',
	    'WISE-4',
	    'pacs_blue',
	    'pacs_red',
	    'spire_1',
	    'spire_2',
	    'spire_3'
	   );


$ABflag = 1; #Use AB System
$KGB::SpecUtils::FastFlux=0; # Take shortcuts for speed?


($z1,$z2) = @ARGV;
die "Usage: expand_massgrid_uds.p z1 z2\n" if $z2 eq "";
print "Will generate expanded PEGASE flux and mass grid from z $z1 to $z2\n\n";

use PDL::Graphics::PGPLOT; 

$H0 = 70;
$OmegaM = 0.3;
$nebfac = 1.0/0.44; # Extra factor of extinction for neb lines
$pc = 3.086E16; # One parsec in meters
$clight = 3e8;#m/s
$flux0 = 3631; # AB zero-point in Jy

$massH = 1.67e-24; #mass of hydrogen in grams
$Lsolar = 3.82e26; #luminosity of the sun in W
$massSun = 1.98e33; #mass of Sun in g
$massHsolar = $massH / $massSun; #mass of hydrogen atom in solar masses

for ($zz=$z1; $zz<=$z2+0.00005; $zz+=0.01) {

  printf "\n\n\n################ Processing Z %10.3f ################\n\n", $zz;
  
  $buff = ""; # Buffer for output results
  
  open FILES, "$ENV{HOME}/Research/PEGASE.2/modelspectra/00index" or die "Can not open 00index file\n"; # List of models
  
  $modelnum = long(0);
  
  while(<FILES>) {
    $file = $_; chomp $file;
    print "Processing file $file    ";
    
    $file =~ /massgrid5.Z(.*)_T(.*)_tb(.*)_logf(.*).pegspec/ ;
    
    $metals = $1; $Tau = $2; $tburst = $3; $logf = $4;
    print "\nZ=$metals T=$Tau tb=$tburst logf=$logf\n";  
    $logf == -9999 if $logf =~ m/'none'/;
    ($t,$Ms, $wav, $starspec,$emspec) = read_peg_spec("$ENV{HOME}/Research/PEGASE.2/modelspectra/".$file,{SPLIT_EMISSION=>1}); # Note #spec returned in W/A/Msolar, t returned in Myr
    
    #add higher wavelength resolution to properly interpolate stellar continuum at lamba > 6.4 microns, matching the 200 A resolution below 6.4 micron
    
    $wav_hr = uniq(qsort(append($wav,64100+200*sequence(5000))));#increase wavelength resolution from 6-30 micron
    
    $starspec_hr = zeroes(nelem($wav_hr),nelem($t));
    $emspec_hr = zeroes(nelem($wav_hr),nelem($t));
    
    $i=0;
    for ($i==0; $i < nelem($t); $i++) {
      $starspec_hr(:,(($i))) .= 10**interpol($wav_hr,$wav,log10($starspec(:,(($i)))));
      $emspec_hr(:,(($i))) .= interpol($wav_hr,$wav,$emspec(:,(($i))));
    }
    
    $logMs = log10($Ms); #match logMs used in generate_FSPS
    
    $logt = log10(1e6*$t); #to match logt in yr used in generate_FSPS
    #($wav, $starspec,$emspec,$totalspec) = rcols("$ENV{HOME}/Research/PopStarModels/$file"); # Note wav in A, spec in Lsolar
    #$dwav = rotate($wav,-1)-$wav; $dwav(-1).=$dwav(-2);
    
    #redden stars - perhaps do later for PAH/NIR emission as well?
  for $A_v (0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0) {
    
    # Redden the spectrum for lambda < 8e4 micron
    $iwav = which $wav_hr > 8e4;
    $attn = peidust('MW', $wav_hr)/peidust('MW', 5500);
    
    $attn($iwav) .= 0;
    
    $spec2 = $starspec_hr * 10**(-0.4 * ($attn * $A_v)) + $emspec_hr * 10**(-0.4 * $attn * $nebfac * $A_v );#apply factor of 2 additional extinction for nebular lines
    
    $DL = lumdist($zz);
    $spec2 /= 4*$Pi * ($DL*1E6*$pc)**2; #Lsolar/m^2/A
    
    # Compute fluxes in filters
    
    @flux = ();
  
	  
    foreach $filt (@filters) {
      badmask($spec2->inplace,0);
      $mfilt = mag($wav_hr,$spec2,$filt); #->((0),:);  # Note unit dim is removed
      $fluxfilt = $flux0 * 10**(-0.4*$mfilt);#flux in units of Jy
      
      badmask($fluxfilt, -99999, $fluxfilt);
      $nnn = sum($fluxfilt < -999);
      
      push @flux, $fluxfilt;
      
    } # Loop over filter
    
    wbuff($modelnum* ones($logt), $logt, $logMs, $metals, $A_v, $Tau, $tburst, $logf, @flux); # Put stuff on output buffer -PEGASE
    #  wbuff($modelnum, $t, $metals, $A_v, @flux); # Put stuff on output buffer - POPSTAR
  }#loop over A_v
    $modelnum++;
  }#loop over files
  close FILES;
  
  #write to buffer
  $fitsfile = "$ENV{HOME}/Research/SEDfitting/models/stars-grid-z-".sprintf("%5.3f",$zz).".fits";
  unlink($fitsfile) if -e $fitsfile;
  wfits $buff(:,0:$counter-1), $fitsfile;
  $counter = 0;
}

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
