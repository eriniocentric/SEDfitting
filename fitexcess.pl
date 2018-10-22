
use PDL;
use lib "$ENV{HOME}/Software";

use PDL::NiceSlice;
use PDL::Graphics::PGPLOT;
use PGPLOT;
use PDL::Graphics::PGPLOTOptions ('set_pgplot_options');
set_pgplot_options('CharSize' => 2.3,'HardCH'=> 1.5,'HardLW'=>2,'AspectRatio'=>1);
use KGB::Cosmology;
use KGB::PegUtils;
use KGB::SpecUtils;
use KGB::Dust;
use PDL::Fit::Polynomial;

($id,$spflag,$spectralclass,$z,$massKarl,$massKarl_err,$K,$Conf,$weight,$sfr2000,$sfrOII,$restUB,$gini,$assym,$fac,$facburst,$burstratio,$burstratio_avg,$burstratio_err,$bbtemp,$temp_avg,$temp_err,$massVIzK,$massVIzKerr,$zmaxVIzK,$massVIzK12,$massVIzK12err,$zmaxVIzK12,$massVIzK1234,$massVIzK1234err,$zmaxVIzK1234,$massVIzK1234_cce,$massVIzK1234_cceerr,$zmaxVIzK1234_cce,$rest3micronflux,$rest3micronflux_err,$massVIzK_nb,$massVIzK_nberr,$massVIzK12_nb,$massVIzK12_nberr,$massVIzK12_nb,$massVIzK12_nberr,$t,$t_err,$tburst,$tburst_err,$stellarburstratio,$stellarburstratio_err,$blendflag,$IRACflag) = rcols("$ENV{HOME}/Software/IRexcess/irexcess_gc.txt");

$idx1 = which($K < 20.6  & $spflag == 0 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $rest3micronflux > $rest3micronflux_err);
$idx2 = which($K < 20.6  & $spflag == 1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $rest3micronflux > $rest3micronflux_err);
$idx3 = which($K < 20.6  & $spflag == 2 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $rest3micronflux > $rest3micronflux_err);
$idx4 = which($K < 20.6  & $spflag == -1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $rest3micronflux > $rest3micronflux_err);

$idxall = which( $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $rest3micronflux > $rest3micronflux_err);


$opt = {Device=> '/xs'};#"$ENV{HOME}/Software/paper/figures/3micron_vs_sfr.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->env(-2.1,2.8,29,34,{Axis=>'LogXY',XTitle=>'SFR / M\d\(2281)\u/yr',YTitle=>'Flux Density\d\gl=3\gmm\u Greybody (W/A)',AspectRatio=>1});
$dlogrest3micronflux_err = 0.4343 * $rest3micronflux_err / $rest3micronflux;

$win->errb(log10($sfrOII($idx1)), log10($rest3micronflux($idx1)),undef,$dlogrest3micronflux_err($idx1),{Symbol=>17,Colour=>'blue',SymbolSize=>2});
$win->errb(log10($sfrOII($idx2)), log10($rest3micronflux($idx2)),undef,$dlogrest3micronflux_err($idx2),{Symbol=>17,Colour=>'green',SymbolSize=>2});
$win->errb(log10($sfrOII($idx3)), log10($rest3micronflux($idx3)),undef,$dlogrest3micronflux_err($idx3),{Symbol=>17,Colour=>'red',SymbolSize=>2});
$win->errb(log10($sfrOII($idx4)), log10($rest3micronflux($idx4)),undef,$dlogrest3micronflux_err($idx4),{Symbol=>17,Colour=>'black',SymbolSize=>2});

$win->close;

#calculate correlation coefficient for BR to SSFR without excluding any points
$ix = which( $burstratio_avg > 0 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag == 0 & ($spflag == 0 | $spflag == 2));

$sfr = $sfrOII($ix);
$ssfr = $sfrOII($ix) / (10**$massVIzK1234_cce($ix)) * 1e9;
$br  = $burstratio_avg($ix);
$br_err = $burstratio_err($ix);
$r3mu = $rest3micronflux($ix);
$r3mu_err = $rest3micronflux_err($ix);
$r3mu_wt = 1 / (0.4343 * $r3mu_err / $r3mu);

$isort = qsorti($ssfr);
$ssfr = $ssfr($isort);
$br = $br($isort);
$br_err = $br_err($isort);
$br_wt = 1 / (0.4343 * $br_err / $br);

$corrcoef_BR = sumover( ($ssfr - average($ssfr)) * ($br-average($br))) / sqrt( sumover(($ssfr - average($ssfr))**2) * sumover(($br - average($br))**2));

#fit a line to BR vs SSFR

($logbr_fit,$coefs) = fitpoly1d(log10($ssfr),log10($br),2,{Weights=> $br_wt});

#monte carlo $br and $ssfr to get errors on fit

$mctrials = 100;
$int_trials = zeroes($mctrials);
$slope_trials = zeroes($mctrials);
$ccBR_mc = zeroes($mctrials);
for ($mc = 0; $mc < $mctrials; $mc++) {

    $rnd = floor(random(nelem($ix))*(nelem($ix)));
    $rnd = sequence(nelem($ix)) if ($mc == 0);
    $br_mc         = $br->dice($rnd);
    $br_wt_mc      = $br_wt->dice($rnd);
    $ssfr_mc       = $ssfr->dice($rnd);

    ($dummy,$coefs_mc) = fitpoly1d(log10($ssfr_mc),log10($br_mc),2,{Weights=> $br_wt_mc});
    $int_trials($mc) .= $coefs_mc(0);
    $slope_trials($mc) .= $coefs_mc(1);
    $ccBR_mc($mc) .= sumover( ($ssfr_mc - average($ssfr_mc)) * ($br_mc-average($br_mc))) / sqrt( sumover(($ssfr_mc - average($ssfr_mc))**2) * sumover(($br_mc - average($br_mc))**2));
  }
dev '/xs';

points log10($ssfr),log10($br);
#points log10($ssfr),log10($br);
hold;
line log10($ssfr),$logbr_fit,{color=>red};

$logbr_err = 0.4343 * $br_err / $br;

$chi2_br = sumover(($logbr_fit-log10($br))**2 / ( $logbr_err**2) ) / (nelem($isort)-2);

#$chi2_br = sumover(($br_fit-$br)**2 / ( $br_err**2) ) / (nelem($isort)-2);

$br_model = 10**($coefs((0))) * $ssfr**($coefs((1)));

line log10($ssfr),log10($br_model),{LineStyle=>Dashed};
close_window;

#fit again for the 3 micron excess

$corrcoef_r3mu = sumover( ($sfr - average($sfr)) * ($r3mu-average($r3mu))) / sqrt( sumover(($sfr - average($sfr))**2) * sumover(($r3mu - average($r3mu))**2));

#monte carlo $br and $ssfr to get errors on fit

$mctrials = 100;
$r3muint_trials = zeroes($mctrials);
$r3muslope_trials = zeroes($mctrials);

for ($mc = 0; $mc < $mctrials; $mc++) {

    $rnd = floor(random(nelem($ix))*(nelem($ix)));
    $rnd = sequence(nelem($ix)) if ($mc == 0);
    $r3mu_mc       = $r3mu->dice($rnd);
    $wt_mc        = $r3mu_wt->dice($rnd);
    $sfr_mc       = $sfr->dice($rnd);

    ($dummy,$coefs_mc) = fitpoly1d(log10($sfr_mc),log10($r3mu_mc),2,{Weights=> $wt_mc});
    $r3muint_trials($mc) .= $coefs_mc(0);
    $r3muslope_trials($mc) .= $coefs_mc(1);
  }

dev '/xs';

points log10($sfr),log10($r3mu);

hold;

$r3mu_model = 10**($coefs_mc((0))) * $sfr**($coefs_mc((1)));

line log10($sfr),log10($r3mu_model),{LineStyle=>Dashed};
close_window;
