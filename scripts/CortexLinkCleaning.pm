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
our @EXPORT = qw(should_keep_link should_keep_link_err_model
                 find_link_cutoff json_hdr_load_contig_hist
                 calc_eff_covg_hist calc_exp_runlens print_expect_table);

sub should_keep_link
{
  my ($link_model, $threshold, $cutoff_dist, $dist, $count) = @_;
  if(!defined($link_model->[$dist])) { die("Wat."); }
  return $count > 1 &&
         ($dist < $cutoff_dist &&
          ($count >= scalar(@{$link_model->[$dist]}) ||
           $link_model->[$dist]->[$count] > $threshold));
}

sub should_keep_link_err_model
{
  my ($link_model, $err_model, $threshold, $dist, $count) = @_;
  return $count > 1 &&
         ($count > scalar(@{$link_model->[$dist]}) ||
          $err_model->[$dist]->[$count] <= $threshold * $link_model->[$dist]->[$count]);
}

sub find_link_cutoff
{
  my ($link_model,$kmer_size,$threshold) = @_;
  my $cutoff = $kmer_size;
  while($cutoff < @$link_model && $link_model->[$cutoff]->[1] < $threshold) { $cutoff++; }
  return $cutoff;
}

sub print_expect_table
{
  my ($expect_tbl, $fh) = @_;
  if(!defined($fh)) { $fh = \*STDOUT; }
  for(my $i = 0; $i < @$expect_tbl; $i++) {
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

sub calc_factorial_arr
{
  my ($len) = @_;

  my $accum = 1;
  my @factorials = (1);

  for(my $i = 1; $i < $len; $i++) {
    $accum *= $i;
    push(@factorials, $accum);
  }

  return @factorials;
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
  if($lambda == 0) { return 0; }
  return (exp(1) ** -$lambda) * sum(map {($lambda ** $_) / factorial($_)} (0..$i));
}

sub poisson_cdf_fast
{
  my ($lambda,$i,$factorial_arr) = @_;
  if($lambda == 0) { return 0; }
  return (exp(1) ** -$lambda) * sum(map {($lambda ** $_) / $factorial_arr->[$_]} (0..$i));
}

sub calc_eff_covg_hist
{
  my (@hist) = @_;
  my @eff_covg = (0) x scalar(@hist);
  for(my $i = 0; $i < @hist; $i++) {
    for(my $j = $i; $j < @hist; $j++) {
      $eff_covg[$i] += ($j-$i+1)*$hist[$j];
    }
  }
  return @eff_covg;
}

# returns 2d array of expectation
# model->[x]->[y] is expectation of y or fewer counts of length x
#   x is of range 0..length($eff_covg_arr)
#   y is of range 0..$max_count
sub calc_exp_runlens
{
  my ($genome_size, $max_count, $eff_covg_arr, $err_rate) = @_;

  my $covg_len = scalar(@$eff_covg_arr);
  my @factorials = calc_factorial_arr($max_count+1);

  if(!defined($err_rate)) { $err_rate = 0.01; }

  my @link_model = ();
  my @err_model = ();

  for(my $i = 0; $i < $covg_len; $i++) {
    my $contig_arrival = $eff_covg_arr->[$i] / $genome_size;
    $link_model[$i] = [map {    poisson_cdf_fast($contig_arrival,           $_, \@factorials)} 0..$max_count];
    $err_model[$i]  = [map {abs($_)} map {1 - poisson_cdf_fast($contig_arrival*$err_rate, $_, \@factorials)} 0..$max_count];
  }

  return (\@link_model,\@err_model);
}

1;
