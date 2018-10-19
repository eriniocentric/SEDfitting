#!/usr/local/bin/perl

#following fitting procedure from Draine et al 2007 p. 881, but use lower Umin values since we do have submm data
#adapted from generateDrainegrid.pl from VNGS M51 project
#grid is redshifted to specified redshift
#filters for SDSS+UKIRT-LAS+WISE+Herschel (note now using SPIRE response curves for point sources unlike M51 project)

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
die "Usage: generateDrainegrid.p z1 z2\n" if $z2 eq "";
print "Will generate grid of Dust models from DL07 from z $z1 to $z2\n\n";


$pc = 3.086E16; # One parsec in meters
$flux0 = 3631; # AB zero-point in Jy

$massH = 1.67e-24; #mass of hydrogen in grams
$massSun = 1.98e33; #mass of Sun in g
$massHsolar = $massH / $massSun; #mass of hydrogen atom in solar masses
#Physical parameters adopted in this paper include a distance of 8.4 ± 0.6 Mpc (Feldmeier et al. 1997), a tilt angle of 20◦ , and a positional angle of the major axis of 170◦ (Tully 1974).

for ($zz=$z1; $zz<=$z2+0.00005; $zz+=0.01) {

     printf "\n\n\n################ Processing Z %10.3f ################\n\n", $zz;

     $buff = ""; # Buffer for output results


     #get dust SED from Draine and Li 2007 models
     
     for $Umin (qw[0.10 0.15 0.20 0.30 0.40 0.50 0.70 0.80 1.00 1.20
		   1.50 2.00 2.50 3.00 4.00 5.00 7.00 8.00 12.0 15.0 20.0 25.0]) {
       $Umax = '1e6';
       
       $ipah = 0;
       for $qpah (0.47,1.12,1.77,2.50,3.19,3.90,4.58) {
	 
	 $ISMdustfile = "$ENV{HOME}/Research/DustModels/DL07spec/U".$Umin."/U".$Umin."_".$Umin."_MW3.1_".sprintf("%1s0",$ipah).".txt";
	 $PDRdustfile = "$ENV{HOME}/Research/DustModels/DL07spec/U".$Umin."/U".$Umin."_".$Umax."_MW3.1_".sprintf("%1s0",$ipah).".txt";
	 
	 
	 ($wavdust,$ISMspecdust,$jnu) = rcols $ISMdustfile, {Lines=>'-1:62'};#wav in micron,flux units of ergs/s/hydrogen nucleon
	 #convert to W/A / 1E10 Mdust(Msolar)
	 $wavdust *= 1e4; #lambda in A
	 
	 ($wavdust2,$PDRspecdust,$jnu) = rcols $PDRdustfile, {Lines=>'-1:62'};#wav in micron,flux units of ergs/s/hydrogen nucleon
	 #convert to W/A / 1E10 Mdust(Msolar)
	 
	 $gamma = 0;
	 for $gamma (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1) {
	   
	   $specdust = (1-$gamma) * $ISMspecdust + $gamma * $PDRspecdust; #units of ergs/s/Hnucleon
	   $specdust /= 0.01*$massHsolar; #erg/s/Mdust
	   $specdust *= 1e-7 / $wavdust; # in W/A/Mdust (Msolar)

	   $DL = lumdist($zz);
	   
	   $specdust /= 4*$Pi * ($DL*1E6*$pc)**2; #W/m^2/A
	   
	   # Compute fluxes in filters
      
	   @flux = ();
	   
	   foreach $filt (@filters) {
	     
	     $mfilt = mag($wavdust,$specdust,$filt); #->((0),:);  # Note unit dim is removed
	     $fluxfilt = $flux0 * 10**(-0.4*$mfilt);#flux in units of Jy/Msolar
	     badmask($fluxfilt, -99999, $fluxfilt);
	     $nnn = sum($fluxfilt < -999);
	     
	     push @flux, $fluxfilt;
	     
	   } # Loop over filter
	   #print @flux;
	   wbuff($Umin, $Umax, $qpah, $ipah, $gamma, @flux); # Put stuff on output buffer
	 } #Loop over gamma ratio of PDR to ISM emission
	 $ipah++;
       }# Loop over PAH ratio qpah
     }#Loop over U radiation intensity 
     
     $fitsfile = "$ENV{HOME}/Research/SEDfitting/models/dust-grid-z-".sprintf("%5.3f",$zz).".fits";
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
