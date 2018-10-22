#find JHK limits for a given stellar mass Star Forming and Evolved galaxy out to z=2

use lib "$ENV{HOME}/Software";
use PDL;
use PDL::NiceSlice;
use KGB::Cosmology;
use KGB::PegUtils;
use KGB::SpecUtils;
use KGB::Dust;
use PDL::Graphics::PGPLOT::Window;
use PDL::Graphics::PGPLOTOptions ('set_pgplot_options');
set_pgplot_options('CharSize' => 1.0,'HardCH'=> 1.4,'HardLW'=>2.5,AspectRatio=>0.8);

$ENV{DATADIR} = "$ENV{HOME}/Software/filters";

$ABflag = 0; #use VEGA mags
$OmegaM = 0.3; $H0 = 70;
$soften = 0.03;
$Law = "SMC";
$nebfac = 2.0; # Extra factor of extinction for neb lines
$pc = 3.086E16; # One parsec in meters
$flux0 = 3631; # AB zero-point in Jy
$clight = 2.99792e8; # speed of light in  meters
$planck = 6.6260755e-34; #planck's constant in m^2 kg /s
$boltzmann = 1.3806503e-23; # m^2 kg s-2 K

# use SED for Tau = 500 Myr galaxy out to z=2.5
$win = PDL::Graphics::PGPLOT::Window->new({Device=> "$ENV{HOME}/Software/Figures/SED_ped_z=0.5,1,1.5_RCS.ps/vcps"});
#$win= PDL::Graphics::PGPLOT::Window->new({Device=>'/xs'});
$win->env(log10(0.0500),log10(10),29,33,{Axis=>'LogXY',Col=>'Black',XTitle=>'\gl (\gmm)',YTitle=>'F\d\gn\u (ergs/s/cm\u2\d)'});

$win->hold;
foreach $z (0,0.5,1,1.5) {

foreach $itime (26) { 
  
  $col = 'blue' if ($itime==6);
  $col = 'green' if ($itime==13);
  $col = 'purple' if ($itime==20);
  $col = 'red' if ($itime==26);
  
  $specfile = "$ENV{HOME}/Software/PEGASE.2/massgrid5.Z0.02_T500_tb200_logf-none.pegspec";
  
   ($t2,$Ms, $wav,$spec,$emspec) = read_peg_spec($specfile,{SPLIT_EMISSION=>1}); # Note #spec returned in W/A
  $spec *= 1E10; $emspec *= 1E10; $Ms *= 1E10; # Convert to 1E10 mass galaxy to avoid rounding
  $lin = 1;
  for $AVV (0.25) {
    
    $spec2 = $spec->dice_axis(1,$itime); 
    $emspec2 = $emspec->dice_axis(1,$itime); 
    $Ms2 = $Ms($itime);
    
    # Redden the spectrum in the same way expand_massgrid2.p does
    $attn = peidust($Law, $wav)/peidust($Law, 5500);
    
    $spec2 = $spec2 * 10**( -0.4 * $attn * $AVV ) + $emspec2 * 10**(-0.4 * $attn * $nebfac * $AVV);
    
    ($wz, $sz) = redshift_spectra($wav, $spec2, pdl($z));
    
    # Coerce 1D
    #      die "Something wrong \n" if (nelem($sz) ne $sz->getdim(0)) or (nelem($wz) ne $wz->getdim(0)); # Should never happen! 
          $wz = $wz->clump(-1)->copy; $sz = $sz->clump(-1)->copy;
    #ignore distance dimming for pedagogical purposes
    #    $DL = lumdist($z);
    
     # $sz /= 4*$Pi * ($DL*1E6*$pc)**2; 
        
#     $win->line(log10($wav(0:1212)* 1e-4),log10($spec2(0:1212)),{Color=>$col,Linestyle=>$lin});
    $win->line(log10($wz* 1e-4),log10($sz ));
    $lin+=2;
  }
}
}
#$win->text("t=300 Myr",-0.85,34.5,{Color=>blue,Charsize=>1.2});
#$win->text("t=1 Gyr",-0.85,33.9,{Color=>green,Charsize=>1.2});
#$win->text("t=3 Gyr",-0.85,32.15,{Color=>purple,Charsize=>1.2});
$win->text("t=13 Gyr",-0.85,30.3,{Color=>red,Charsize=>1.2});

if (1) {
($wv,$ffilt) = rcols('Software/filters/sys_mosaic_U.dat');
$wv /= 1e4;#convert to microns
$ffilt *= 1e35;
#$win->line(log10($wv),log10($ffilt),{Color=>'purple'});

($wv,$ffilt) = rcols('Software/filters/sys_mosaic_B.dat');
$wv /= 1e4;#convert to microns
$ffilt *= 1e35;
#$win->line(log10($wv),log10($ffilt),{Color=>'blue'});

($wv,$ffilt) = rcols('Software/filters/sys_mosaic_V.dat');
$wv /= 1e4;#convert to microns
$ffilt *= 1e33;
$win->line(log10($wv),log10($ffilt),{Color=>'green'});

($wv,$ffilt) = rcols('Software/filters/sys_mosaic_R.dat');
$wv /= 1e4;#convert to microns
$ffilt *= 1e35;
#$win->line(log10($wv),log10($ffilt),{Color=>'red'});

($wv,$ffilt) = rcols('Software/filters/sys_mosaic_I.dat');
$wv /= 1e4;#convert to microns
$ffilt *= 1e33;
$win->line(log10($wv),log10($ffilt),{Color=>'red'});

#($wv,$ffilt) = rcols('Software/filters/sys_circi_J.dat');
$wv /= 1e4;#convert to microns
$ffilt *= 1e35;
#$win->line(log10($wv),log10($ffilt),{Color=>'blue'});

#($wv,$ffilt) = rcols('Software/filters/sys_circi_H.dat');
$wv /= 1e4;#convert to microns
$ffilt *= 1e35;
#$win->line(log10($wv),log10($ffilt),{Color=>'green'});

($wv,$ffilt) = rcols('Software/filters/sys_circi_K.dat');
$wv /= 1e4;#convert to microns
$ffilt *= 1e35;
#$win->line(log10($wv),log10($ffilt),{Color=>'red'});
}

$win->close();
