#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -euo pipefail
set +o posix

cmd=$0

function usage {
  >&2 echo "usage $cmd [--simplify|--dot] <kmer> <file1> [...]"
  >&2 echo "  prints pdf to stdout, so please remember to redirect"
  >&2 echo "  e.g. $cmd 5 <(echo ACAACACGT) <(echo CCACACAA) > out.pdf"
  exit -1
}

script_args=
mkpdf=1

while [[ $# -gt 2 ]]
do
  if [[ ($1 == "--simplify") ]]
  then
    script_args=$1
    shift
  elif [[ $1 == "--dot" ]]
  then
    mkpdf=0
    shift
  else
    usage
  fi
done

if [[ $# -ne 2 || !( $1 =~ ^[0-9]+$ ) ]]
then
  usage
fi

kmer=$1
shift

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
MCCORTEX="$DIR/bin/mccortex"
CTX2GRAPHVIZ="$DIR/scripts/perl/mccortex-graph-to-graphviz.pl"
if [[ !(-e $MCCORTEX) || !(-x $MCCORTEX) ]]
then
  echo "Did you compile McCortex? I cannot run `$MCCORTEX`"
  exit -1
fi

files=$(printf " --seq %s" $@; printf "\n")

if [[ $mkpdf == 1 ]]; then
  $MCCORTEX $kmer build -q -k $kmer --sample seq2pdf $files - | \
    $CTX2GRAPHVIZ -k $kmer $script_args - | \
    dot -Tpdf
else
  $MCCORTEX $kmer build -q -k $kmer --sample seq2pdf $files - | \
    $CTX2GRAPHVIZ -k $kmer $script_args -
fi
