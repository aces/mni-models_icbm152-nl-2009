#! /usr/bin/env perl
#
# Copyright Alan C. Evans
# Professor of Neurology
# McGill University
#

# Extract a surface at 5120 polygons for a template mask.

use strict;
use warnings "all";
use File::Basename;
use File::Spec;
use File::Temp qw/ tempdir /;

use Getopt::Tabular;
use MNI::Startup;
use MNI::FileUtilities;
use MNI::DataDir;

# --- set the help & usage strings ---
my $help = <<HELP;
Required parameters:
  wm_mask.mnc  : white matter mask
  white.obj    : white matter surface (output)
HELP

my $license = <<LICENSE;
Copyright Alan C. Evans
Professor of Neurology
McGill University

LICENSE

my $usage = <<USAGE;
Usage: $ProgramName wm_mask.mnc white.obj [options]
       $ProgramName -help to list options

$license
USAGE

Getopt::Tabular::SetHelp( $help, $usage );

my $subsample = 0;

my @options = ( 
  );

GetOptions( \@options, \@ARGV ) or exit 1;
die "$usage\n" unless @ARGV == 2;

my $head_mask = shift;
my $surface = shift;

if( !( -e $head_mask ) ) {
  print "$usage\n";
  die "White matter mask must exist.\n";
}

## Note: this must be the old ellipsoid from civet-1.1.11.
my $ICBM_white_model = MNI::DataDir::dir("surface-extraction") .
                       "/white_model_320.obj";
# my $ICBM_white_model = "/export/home/claude/QUARANTINE/Linux-x86_64/share/surface-extraction/white_model_320.obj";

my $tmpdir = &tempdir( "mcubes-XXXXXX", TMPDIR => 1, CLEANUP => 1 );

&run( 'mincdefrag', $head_mask, "${tmpdir}/mask.mnc", 1, 6 );

# Run ASP to fit surface on head mask.

&run( 'cp', '-f', $ICBM_white_model, $surface );
&run_asp( $surface, "${tmpdir}/mask.mnc", 0.5, $ICBM_white_model,
          $tmpdir );

print "done\n";


# Run ASP to fit surface on head mask up to 5120 polygons only.

sub run_asp {

  my $white = shift;
  my $wm_mask = shift;
  my $isovalue = shift;
  my $white_model = shift;
  my $tmpdir = shift;

  my $self_dist2 = 0.01;
  my $n_selfs = 9;

  my $stop_threshold = 1e-3;
  my $stop_iters = 10;

  my $n_per = 5;
  my $tolerance = 1.0e-03;
  my $f_tolerance = 1.0e-06;
  my $oo_scale = 0.5;

  my @schedule = (
  #  size   sw    cw  n_itr  inc   offo   offi   si  in  out  over   sw   self
  #  ----- ----  ---- -----  ----  ----   ----   --  --  ---  ----  ----  ----
      320, 1e3,   0,    100,  50,   100,    1,   2,  50,  3,   40,  1e0,  1,
     1280, 1e3,   0,    100,  50,   100,    1,   2,  50,  3,   20,  1e0,  1,
     1280, 1e2,   0,    200,  50,   100,    1,   1,  50,  3,   20,  1e0,  5,
     5120, 1e2,   0,    200,  50,    50,    2, 1.5,  50,  3,   10,  1e0,  3,
     5120, 1e1,   0,    200,  50,    50,    5, 1.5,  25,  3,   10,  1e0,  3,
     5120, 5e0,   0,    200,  50,    50,    5, 1.5,  25,  3,   10,  1e0,  3,
  );
  my $sched_size =  13;
  my $num_steps = @schedule / $sched_size;

  &run( 'cp', '-f', $white_model, "${tmpdir}/white_model_tmp.obj" );
  $white_model = "${tmpdir}/white_model_tmp.obj";

  for( my $i = 0;  $i < @schedule;  $i += $sched_size ) {
    my ( $size, $sw, $cw, $n_iters, $iter_inc, $offo, $offi,
         $si_step, $in_dist, $out_dist, $oversample, $self_weight, 
         $self_dist ) = @schedule[$i..$i+$sched_size-1];

    &run( 'subdivide_polygons', $white_model, $white_model, $size );

    my $prev_size = `print_n_polygons $white`;
    chomp( $prev_size );
    if( $prev_size != $size ) {
      &run( 'subdivide_polygons', $white, $white, $size );
      &run( 'subdivide_polygons', $white_model, $white_model, $size );
    }

    $oversample *= $oo_scale;
    my $self2 = get_self_intersect( $self_weight, $n_selfs, $self_dist,
                                    $self_dist2 );

    my $b2 = " -boundary $offo $offi $wm_mask " .
             " $isovalue - $out_dist $in_dist 0 0 $oversample ";

    for( my $iter = 0;  $iter < $n_iters;  $iter += $iter_inc ) {
      print "echo Step ${size}: $iter / $n_iters    $sw\n";

      my $ni = $n_iters - $iter;
      if( $ni > $iter_inc )  { $ni = $iter_inc; }

      my $command = "surface_fit -mode two -surface $white $white " .
                    " -stretch $sw $white_model -.9 0 0 0 " .
                    " $b2 $self2 -step $si_step " .
                    " -fitting $ni $n_per $tolerance " .
                    " -ftol $f_tolerance " .
                    " -stop $stop_threshold $stop_iters ";
      system( $command );
    }
  }
  unlink( $white_model );
}

# from surface-extraction/deform_utils.pl

sub get_self_intersect( $$$$ ) {

    my( $self_weight, $n_selfs, $self_dist, $self_dist2 ) = @_;
    my( $self, $weight, $s, $dist );

    if( $self_weight > 0.0 ) {
        $self = "";
        $weight = $self_weight;

        for( $s = 0;  $s < $n_selfs;  ++$s ) {
            $dist = $self_dist + ($self_dist2 - $self_dist) *
                    $s / $n_selfs;
            $self = $self . " -self_intersect $weight $dist ";
            $weight *= 10.0;
        }
        $self = $self . " -self_intersect 1e8 $self_dist2 ";
    } else { 
        $self = ""; 
    }
    $self;
}


# Execute a system call.

sub run {
  print "@_\n";
  system(@_)==0 or die "Command @_ failed with status: $?";
}

