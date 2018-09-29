#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;
use File::Path qw(make_path);
use List::Util qw(min max sum);
use JSON;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin."/perl/";

use McCortexBreakpoints;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  
  print STDERR <<EOF;
Usage: ./make-circos.pl <outdir> [breakpoints.txt]

  Generate a circos plot for a breakpoint file. Can pipe in breakpoints if not
  given on the command line. Will create <outdir> and parent directories if
  they do no already exist.

  Example:

    gzip -dc breakpoints.txt.gz | ./make-circos.pl my-circos-plot -

  Output will be: my-circos-plot/circos.png

EOF

  exit(-1);
}


# Replace CHROM_UNIT with 1000 or 1000000 (1K or 1M)
my $circos_conf = <<EOF;
<<include etc/colors_fonts_patterns.conf>>
<<include ideogram.conf>>
<<include ticks.conf>>
<image>
<<include etc/image.conf>>
</image>

karyotype   = karyotype.txt
chromosomes_units           = CHROM_UNIT
chromosomes_display_default = yes

<links>
z             = 0
radius        = 0.85r
bezier_radius = 0.6r
thickness     = 2
show          = yes
color         = black
ribbon        = yes

<link>
file          = breakpoints.fw.txt
<rules>
<rule>
  condition   = not(var(intrachr))
  radius      = 0.7r
</rule>
</rules>
</link>

<link>
file          = breakpoints.rv.txt
<rules>
<rule>
  condition   = not(var(intrachr))
  radius      = 0.7r
</rule>
</rules>
</link>

</links>

track_width = 0.04
track_pad   = 0.00
track_start = 0.98

<plots>

type = heatmap
stroke_thickness = 0

# SNPs
<plot>
<<include r0r1.conf>>
color = white,red_a5,red_a4,red_a3,red_a2,red_a1,red
file  = mnps.txt
# min = 0
# max = 1000000000
</plot>

# Insertions
<plot>
<<include r0r1.conf>>
color = white,blue_a5,blue_a4,blue_a3,blue_a2,blue_a1,blue
file  = insertions.txt
# min = 0
# max = 1000000000
</plot>

# Deletions
<plot>
<<include r0r1.conf>>
color = white,black_a5,black_a4,black_a3,black_a2,black_a1,black
file  = deletions.txt
# min = 0
# max = 1000000000
</plot>

</plots>

<<include etc/housekeeping.conf>>
EOF

# save to <out>/ideogram.conf
my $ideogram_conf = <<EOF;
<ideogram>

<spacing>
default = 0.01r
break   = 2u
</spacing>

<<include ideogram.position.conf>>
<<include ideogram.label.conf>>
<<include bands.conf>>

</ideogram>
EOF

# save to <out>/ideogram.label.conf
my $ideogram_label_conf = <<EOF;
show_label       = yes
label_font       = default
label_radius     = dims(ideogram,radius) + 0.15r
label_with_tag   = yes
label_size       = 21
label_format     = eval(var(chr))
EOF

# label_parallel   = yes
# label_case       = lower
# label_format     = eval(sprintf("chr%s",var(name)))
# eval(sprintf("chr%s",var(label)))

# save to <out>/ideogram.position.conf
my $ideogram_position_conf = <<EOF;
radius           = 0.65r
thickness        = 25p
fill             = yes
fill_color       = black
stroke_thickness = 2
stroke_color     = black
EOF

# save to <out>/ideogram.conf
my $bands_conf = <<EOF;
show_bands            = yes
fill_bands            = yes
band_stroke_thickness = 2
band_stroke_color     = white
band_transparency     = 0
EOF

# Replace MULTIPLIER_VALUE,MULTIPLIER_UNIT with 1e-3,K or 1e-6,M
# save to <out>/ticks.conf
my $ticks_conf = <<EOF;
show_ticks          = yes
show_tick_labels    = yes

show_grid          = yes
grid_start         = 0.5r
grid_end           = 1.0r

<ticks>
skip_first_label     = no
skip_last_label      = no
radius               = dims(ideogram,radius_outer)
tick_separation      = 2p
min_label_distance_to_edge = 10p
label_separation = 5p
label_offset     = 5p
multiplier       = MULTIPLIER_VALUE
color            = black

<tick>
spacing        = 1u
size           = 8p
thickness      = 2p
show_label     = no
</tick>

<tick>
spacing        = 5u
size           = 12p
thickness      = 2p
show_label     = yes
label_size     = 20p
format         = %dMULTIPLIER_UNIT
</tick>

<tick>
spacing        = 10u
size           = 14p
thickness      = 2p
show_label     = yes
label_size     = 24p
format         = %dK
</tick>

</ticks>
EOF

my $r0r1_conf = <<EOF;
# set track radius values based on track counter
r1  = eval(sprintf("%fr",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))))
r0  = eval(sprintf("%fr",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))-conf(track_width)))
EOF


main(@ARGV);
exit(0);

sub open_or_die
{
  my ($path) = @_;
  my $fh;
  open($fh, ">$path") or die("Cannot open file: $path");
  return $fh;
}

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}

sub main
{
  if(@_ != 2) { print_usage(); }
  my ($outdir,$file) = @_;

  # Input file
  my $fh;

  # Output files
  my ($brkfw_fh, $brkrv_fh, $mnp_fh, $ins_fh, $del_fh);
  my ($circos_fh, $karyotype_fh);
  my ($idgrm_fh, $idgrm_label_fh, $idgrm_pos_fh);
  my ($bands_fh, $ticks_fh, $r0r1_fh);

  open($fh, $file) or die("Cannot read file $file");
  my $cb = new McCortexBreakpoints($fh);

  # 0. Create output directory
  if($outdir eq "") {
    die("Empty output directory path is a bad idea");
  }
  if($outdir eq "." || $outdir eq "..") {
    die("Current working directory not a wise outdir");
  }

  if(-e $outdir) {
    if(!(-d $outdir)) {
      die("File already exists with name of output dir: $outdir");
    }
  } else {
    make_path($outdir) or die("Cannot make output directory: $outdir");
  }

  # 1. Open output files
  $brkfw_fh = open_or_die("$outdir/breakpoints.fw.txt");
  $brkrv_fh = open_or_die("$outdir/breakpoints.rv.txt");
  $mnp_fh = open_or_die("$outdir/mnps.txt");
  $ins_fh = open_or_die("$outdir/insertions.txt");
  $del_fh = open_or_die("$outdir/deletions.txt");
  $circos_fh = open_or_die("$outdir/circos.conf");
  $karyotype_fh = open_or_die("$outdir/karyotype.txt");
  $idgrm_fh = open_or_die("$outdir/ideogram.conf");
  $idgrm_label_fh = open_or_die("$outdir/ideogram.label.conf");
  $idgrm_pos_fh = open_or_die("$outdir/ideogram.position.conf");
  $bands_fh = open_or_die("$outdir/bands.conf");
  $ticks_fh = open_or_die("$outdir/ticks.conf");
  $r0r1_fh = open_or_die("$outdir/r0r1.conf");

  # 2. Load ref coordinates from header of breakpoint file
  my %chrs = ();

  # Parse contigs from header
  my ($idx, $sum_lengths) = (1,0);

  my $hdr_txt = $cb->{'_header'};
  my $hdr_json = decode_json($hdr_txt);
  my $contigs_hdr = $hdr_json->{'commands'}[0]->{'breakpoints'}->{'contigs'};

  for my $hdr (@$contigs_hdr) {
    $chrs{$hdr->{'id'}} = {'ID' => $hdr->{'id'},
                           'length' => $hdr->{'length'},
                           'idx' => $idx};
    $idx++;
    $sum_lengths += $hdr->{'length'};
  }

  print "Genome size: $sum_lengths\n";
  my ($binsize, $chrom_units, $scale);

  if($sum_lengths >= 10**9) { $scale = 6; }
  else { $scale = 3; }

  $chrom_units = 10 ** int(log10($sum_lengths/100));
  $binsize = $chrom_units/2;
  # $binsize = 10**$scale;

  print "  histogram binsize: $binsize\n";
  print "  chrom units: $chrom_units\n";

  make_data_files($cb, $brkfw_fh, $brkrv_fh, $mnp_fh, $ins_fh, $del_fh, $binsize, \%chrs);

  close($fh);

  make_circos_file($circos_fh, $chrom_units);
  make_karyotype_file($karyotype_fh, \%chrs);
  make_ideograph_files($idgrm_fh, $idgrm_label_fh, $idgrm_pos_fh);
  make_bands_file($bands_fh);
  make_ticks_file($ticks_fh, $scale);

  print $r0r1_fh $r0r1_conf;
  close($r0r1_fh);

  print "To create $outdir/circos.png:

      cd $outdir
      circos

";
}

#
# Functions and constants
#


# chr - chr1 1 0 249250621 <colour>
# chr - chr2 2 0 243199373 chr2
# chr - chr3 3 0 198022430 chr3
# chr - chr4 4 0 191154276 chr4
sub make_karyotype_file
{
  my ($karyotype_fh,$chrhash) = @_;
  my @chrs = sort {$a->{'idx'} <=> $b->{'idx'}} values %$chrhash;
  my $num_chr = @chrs;
  for(my ($i,$j) = (0, 1); $i < $num_chr; $i++, $j++) {
    print $karyotype_fh "chr - $chrs[$i]->{'ID'} $j 0 $chrs[$i]->{'length'} white\n";
  }
  close($karyotype_fh);
}

# <out>/breakpoints.fw.txt
# <out>/breakpoints.rv.txt
#   chrFrom 200000 200000 chrTo 150000 150000
#
# <out>/mnps.txt
# <out>/insertions.txt
# <out>/deletions.txt
#   strep1 0 49999 2
#   strep1 50000 99999 1000
#   strep1 100000 149999 3000
#   strep1 150000 199999 4000
#   strep1 200000 249999 5000
#   strep1 250000 1796226 6954

sub dump_hist
{
  my ($fh,$hist,$binsize,$chrhash) = @_;
  # strep1 0 49999 2
  my ($max,$sum) = (0,0);
  my @chrs = sort {$a->{'idx'} <=> $b->{'idx'}} values %$chrhash;
  for my $chr (@chrs) {
    my $chrhist = $hist->[$chr->{'idx'}];
    for(my ($i,$pos) = (0, 0); $pos < $chr->{'length'}; $pos += $binsize, $i++) {
      my $end = min($pos+$binsize-1, $chr->{'length'});
      print $fh "$chr->{'ID'} $pos $end ".$chrhist->[$i]."\n";
    }
    $max = max($max, @$chrhist);
    $sum += sum(@$chrhist);
  }
  return ($max,$sum);
}

# Closes FILEHANDLEs passed, but not breakpoint file we are reading
sub make_data_files
{
  my ($cb, $brkfw_fh, $brkrv_fh, $mnp_fh, $ins_fh, $del_fh, $binsize, $chrhash) = @_;
  my ($seq5p, $seq3p, $pathseq, $flank5p_refs, $flank3p_refs, $cols, $callid);

  print "Creating data files...\n";

  # For more information on colours see: brewer-palettes-swatches.pdf
  #
  # light blue, dark blue, light green, dark green
  # my @cols = qw(paired-4-qual-1 paired-4-qual-2 paired-4-qual-3 paired-4-qual-4);
  #
  # light orange, dark orange, light purple, dark purple
  # my @cols = qw(spectral-10-div-7 spectral-10-div-8 spectral-10-div-9 spectral-10-div-10);
  #
  # dark blue, dark green, dark orange, dark purple
  # my @cols = qw(spectral-10-div-2 spectral-10-div-4 spectral-10-div-8 spectral-10-div-10);

  # my @cols = qw(spectral-10-div-2 spectral-10-div-4 spectral-10-div-7 spectral-10-div-9);
  #
  # my @cols = qw(plyg-10-div-2 plyg-10-div-4 plyg-10-div-7 plyg-10-div-9);
  # light/dark blue/red
  # my @cols = qw(rdbu-10-div-2 rdbu-10-div-4 rdbu-10-div-7 rdbu-10-div-9);
  # purple/green Green is inversion
  my @cols = qw(prgn-10-div-2 prgn-10-div-4 prgn-10-div-7 prgn-10-div-9);

  print "  Key: light/dark purple - Regular breakpoint\n";
  print "       light/dark green  - Inversion breakpoint\n";

  my ($num_skipped,$num_entries) = (0,0);
  my ($fw_breakpoints, $rv_breakpoints) = (0,0);
  my @mnphist = ();
  my @inshist = ();
  my @delhist = ();
  map {$mnphist[$_->{'idx'}] = [map {0} (0..$_->{'length'})]} values %$chrhash;
  map {$inshist[$_->{'idx'}] = [map {0} (0..$_->{'length'})]} values %$chrhash;
  map {$delhist[$_->{'idx'}] = [map {0} (0..$_->{'length'})]} values %$chrhash;

  for($num_entries = 0; ; $num_entries++)
  {
    ($seq5p, $seq3p, $pathseq, $flank5p_refs, $flank3p_refs, $cols, $callid) = $cb->next();
    if(!defined($seq5p)) { last; }

    my $num5p = scalar(@$flank5p_refs);
    my $num3p = scalar(@$flank3p_refs);

    if($num5p != 1 || $num3p != 1) { $num_skipped++; }
    else {
      my $flank5p = $flank5p_refs->[0];
      my $flank3p = $flank3p_refs->[0];
      my $reflen = $flank3p->{'start'} - $flank5p->{'end'} - 1;
      my $pathlen = length($pathseq);

      my ($chr5p, $chr3p) = map {$chrhash->{$_->{'chrom'}}} ($flank5p, $flank3p);
      my $start5p = $flank5p->{'start'};
      my $pos5p = $flank5p->{'end'};
      my $pos3p = $flank3p->{'start'};
      my $end3p = $flank3p->{'end'};

      # link colour to show direction
      my $lcol = "";

      my $intrachr = ($chr5p->{'idx'} == $chr3p->{'idx'}); # ignore for now
      my $first_link = ($chr5p->{'idx'} < $chr3p->{'idx'} ||
                        ($chr5p->{'idx'} == $chr3p->{'idx'} && $pos5p < $pos3p));

      my $inversion = ($flank5p->{'strand'} ne $flank3p->{'strand'});

      # Two color scheme
      if(!$inversion) { $lcol = "color=$cols[1]"; }
      else { $lcol = "color=$cols[3]"; }

      # if($first_link && !$inversion) { $lcol = "color=$cols[0]"; }
      # if(!$first_link && !$inversion) { $lcol = "color=$cols[1]"; }
      # if($first_link && $inversion) { $lcol = "color=$cols[2]"; }
      # if(!$first_link && $inversion) { $lcol = "color=$cols[3]"; }

      # if($intrachr) {
      #   if($inversion) { $lcol = "color=red,z=1"; }
      #   else { $lcol = "color=black,z=0"; }
      # }
      # else {
      #   if($first_link) { $lcol = "color=orange"; }
      #   else { $lcol = "color=green"; }
      # }

      if($flank5p->{'strand'} ne $flank3p->{'strand'}) {
        # rv breakpoint
        print $brkrv_fh "$chr5p->{'ID'} $start5p $pos5p $chr3p->{'ID'} $pos3p $end3p $lcol\n";
        $rv_breakpoints++;
      }
      elsif($flank5p->{'chrom'} ne $flank3p->{'chrom'}) {
        # fw breakpoint
        print $brkfw_fh "$chr5p->{'ID'} $start5p $pos5p $chr3p->{'ID'} $pos3p $end3p $lcol\n";
        $fw_breakpoints++;
      }
      elsif($pathlen == $reflen) {
        # mnps
        $mnphist[$chr5p->{'idx'}]->[int($pos5p/$binsize)]++;
      }
      elsif(abs($pathlen-$reflen) < 200+min($pathlen,$reflen)) {
        if($pathlen < $reflen) {
          # deletion
          $delhist[$chr5p->{'idx'}]->[int($pos5p/$binsize)]++;
        } else {
          #insertion
          $inshist[$chr5p->{'idx'}]->[int($pos5p/$binsize)]++;
        }
      }
      else
      {
        # fw breakpoint
        print $brkfw_fh "$chr5p->{'ID'} $start5p $pos5p $chr3p->{'ID'} $pos3p $end3p $lcol\n";
        $fw_breakpoints++;
      }
    }
  }

  printf("  Skipped $num_skipped / $num_entries (%.2f%%)\n",
         100*$num_skipped / $num_entries);

  my ($mnp_max, $mnp_sum) = dump_hist($mnp_fh,\@mnphist,$binsize,$chrhash);
  my ($ins_max, $ins_sum) = dump_hist($ins_fh,\@inshist,$binsize,$chrhash);
  my ($del_max, $del_sum) = dump_hist($del_fh,\@delhist,$binsize,$chrhash);

  print "    $mnp_sum MNPs (max: $mnp_max)\n";
  print "    $ins_sum INSs (max: $ins_max)\n";
  print "    $del_sum DELs (max: $del_max)\n";
  print "    $fw_breakpoints FW Breakpoints\n";
  print "    $rv_breakpoints RV Breakpoints\n";

  close($brkfw_fh);
  close($brkrv_fh);
  close($mnp_fh);
  close($ins_fh);
  close($del_fh);

  return ($mnp_max, $ins_max, $del_max, max($fw_breakpoints, $rv_breakpoints));
}

sub make_circos_file
{
  my ($circos_fh, $chrom_units) = @_;
  print "Creating circos.conf\n";
  my $circos = $circos_conf;
  $circos =~ s/CHROM_UNIT/$chrom_units/g;
  print $circos_fh $circos;
  close($circos_fh);
}

# Closes FILEHANDLEs passed
sub make_ideograph_files
{
  my ($idgrm_fh, $idgrm_label_fh, $idgrm_pos_fh) = @_;
  print "Creating ideogram.conf\n";
  print $idgrm_fh $ideogram_conf;
  print "Creating ideogram.label.conf\n";
  print $idgrm_label_fh $ideogram_label_conf;
  print "Creating ideogram.position.conf\n";
  print $idgrm_pos_fh $ideogram_position_conf;

  close($idgrm_fh);
  close($idgrm_label_fh);
  close($idgrm_pos_fh);
}

# Closes FILEHANDLE passed
sub make_bands_file
{
  my ($bands_fh) = @_;
  print $bands_fh $bands_conf;
  close($bands_fh);
}

# Closes FILEHANDLE passed
sub make_ticks_file
{
  my ($ticks_fh,$scale) = @_;
  print "Creating ticks.conf\n";

  my $ticks = $ticks_conf;
  my $unit;

  if($scale == 3) { $unit = 'K'; }
  elsif($scale == 6) { $unit = 'M'; }
  else { die('Bad scale - not K or M (10^3 or 10^6)'); }

  $ticks =~ s/MULTIPLIER_VALUE/1e-$scale/g;
  $ticks =~ s/MULTIPLIER_UNIT/$unit/g;
  print $ticks_fh $ticks;
  close($ticks_fh);
}
