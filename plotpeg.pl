use PDL;
use PDL::NiceSlice;
use lib "$ENV{HOME}/Software";
use KGB::PegUtils;
use PDL::Graphics::PGPLOT;
use PDL::Graphics::PGPLOTOptions ('set_pgplot_options');
set_pgplot_options('CharSize' => 1.0,'HardCH'=> 1.1,'HardLW'=>2.5);

$tburst = qw(5000);
$tau = qw(500);
$logf = qw(-1.12);
$frac = sprintf("%4.3f",10**$logf);
($t,$Ms,$wav,$spec,$emspec) = read_peg_spec("$ENV{HOME}/Software/Massgrid5/massgrid5.Z0.02_T500_tb5000_logf-1.12.pegspec");
$spec *=1E10; $empspec *= 1E10; #Convert to 1E10 to avoid rounding

$start = which($wav >=500);
$end   = which($wav >= 100000);
$start = $start(0);
$end   = $end(0);

$logwav = log10($wav($start:$end));
system "rm $ENV{HOME}/ani*.ps";
#$win=PDL::Graphics::PGPLOT::Window->new(Device=>"$ENV{HOME}/ani.ps/cps");

for ($itime=0;$itime<=(nelem($t)-1);$itime++) {
    $age = $t(($itime));
    $win=PDL::Graphics::PGPLOT::Window->new(Device=>"$ENV{HOME}/ani$itime.ps/cps");
    $win->env(log10(500),log10(100000),max(log10($spec))-7,max(log10($spec))+0.2,{Axis=>'LogXY',XTitle=>'Wavelength(\A)',YTitle=>'Flux',Title=>"Tau = $tau; Tburst = $tburst; BurstMassFraction = $frac",AxisColour=>Black});
    $win->text("Age is $age Myrs",log10(10000),max(log10($spec))-1,{Colour=>Black});
    $spec2 = $spec->dice_axis(1,$itime); 
    $emspec2 = $emspec->dice_axis(1,$itime); 
    $flux = $spec2+$emspec2;
    $flux = $flux->clump(-1)->copy;
    $x = which($flux == 0);
    $flux($x) .= 1;
    $logflux = log10($flux($start:$end));
    $win->line($logwav,$logflux,{Colour=>red});
    $win->close();
    system "convert -rotate 90 ani$itime.ps ani$itime.pdf";
}
#$win->close();
