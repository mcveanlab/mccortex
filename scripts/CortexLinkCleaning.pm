package CortexLinkCleaning;

use strict;
use warnings;
use Carp;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use JSON;
use List::Util qw(min max sum);

use base 'Exporter';
our @EXPORT = qw(should_keep_link find_link_cutoff json_hdr_load_contig_hist
                 calc_exp_runlens print_expect_table);

sub should_keep_link
{
  my ($expect_tbl, $max_count, $threshold, $cutoff_dist, $dist, $count) = @_;
  if(!defined($expect_tbl->[$dist])) { die("Wat."); }
  return ($dist < $cutoff_dist &&
          ($count > $max_count || $expect_tbl->[$dist]->[$count] > $threshold));
}

sub find_link_cutoff
{
  my ($expect_tbl,$kmer_size,$threshold) = @_;
  my $cutoff = $kmer_size;
  while($cutoff < @$expect_tbl && $expect_tbl->[$cutoff]->[1] < $threshold) { $cutoff++; }
  return $cutoff;
}

sub print_expect_table
{
  my ($expect_tbl, $kmer_size, $fh) = @_;
  if(!defined($fh)) { $fh = \*STDOUT; }
  for(my $i = $kmer_size; $i < @$expect_tbl; $i++) {
    print $fh $i."  ".join(' ', map {sprintf("%.3f",$_)} @{$expect_tbl->[$i]})."\n";
  }
}

sub json_hdr_load_contig_hist
{
  my ($hdr,$colour) = @_;

  my $json = (ref($hdr) ne "" ? $hdr : decode_json($hdr));
  my $length_arr = $json->{'paths'}->{'contig_hists'}->[$colour]->{'lengths'};
  my $counts_arr = $json->{'paths'}->{'contig_hists'}->[$colour]->{'counts'};

  if(!defined($length_arr)) { die("Missing lengths array [$colour]"); }
  if(!defined($counts_arr)) { die("Missing counts array [$colour]"); }

  if(@$length_arr != @$counts_arr) { die("Mismatch [$colour]: @$length_arr != @$counts_arr"); }

  my $max_len = max(@$length_arr);
  my @contigs = map {0} 0..$max_len;

  for(my $i = 0; $i < @$length_arr; $i++) {
    $contigs[$length_arr->[$i]] += $counts_arr->[$i];
  }

  return \@contigs;
}

sub factorial
{
  my ($i) = @_;
  my $ret = 1;
  for(; $i > 1; $i--) { $ret *= $i; }
  return $ret;
}

sub combinatorial
{
  my ($n,$k) = @_;
  return factorial($n)/(factorial($k)*factorial($n-$k));
}

sub combinatorial_fast
{
  my ($n,$k,$factorial_arr) = @_;
  return $factorial_arr->[$n]/($factorial_arr->[$k]*$factorial_arr->[$n-$k]);
}

sub binomial_pdf
{
  my ($n,$k,$p) = @_;
  return combinatorial($n,$k) * $p**$k * (1-$p)**($n-$k);
}

sub binomial_pdf_fast
{
  my ($n,$k,$p,$factorial_arr) = @_;
  return combinatorial_fast($n,$k,$factorial_arr) * $p**$k * (1-$p)**($n-$k);
}

sub binomial_cdf
{
  my ($n,$k,$p) = @_;
  my $binom = 0;
  for(my $i = 0; $i <= $k; $i++) { $binom += binomial_pdf($n,$k,$p); }
}

sub binomial_cdf_fast
{
  my ($n,$k,$p,$factorial_arr) = @_;
  my $binom = 0;
  for(my $i = 0; $i <= $k; $i++) { $binom += binomial_pdf_fast($n,$k,$p,$factorial_arr); }
}

sub poisson_cdf
{
  my ($lambda,$i) = @_;
  return exp(1) ** -$lambda * sum(map {($lambda ** $_) / factorial($_)} (0..$i));
}

sub poisson_cdf_fast
{
  my ($lambda,$i,$factorial_arr) = @_;
  return exp(1) ** -$lambda * sum(map {($lambda ** $_) / $factorial_arr->[$_]} (0..$i));
}

# returns 2d array of expectation
# exp[x][y] is expectation of y or fewer counts of length x
sub calc_exp_runlens
{
  my ($genome_size,$k,$maxcount,@hist) = @_;
  my @counts = (0) x scalar(@hist);
  for(my $i = $k; $i < @hist; $i++) {
    for(my $j = $i; $j < @hist; $j++) {
      $counts[$i] += ($j-$i+1)*$hist[$j];
    }
  }

  my $accum = 1;
  my @factorials = (1);
  for(my $i = 1; $i < @hist; $i++) {
    $accum *= $i;
    push(@factorials, $accum);
  }

  # print "Contig counts: @counts\n";
  my @exp = ();
  for(my $i = $k; $i < @hist; $i++) {
    $exp[$i] = [map {poisson_cdf_fast($counts[$i] / $genome_size, $_, \@factorials)} 0..$maxcount];
  }
  return (\@exp,\@counts);
}

1;
