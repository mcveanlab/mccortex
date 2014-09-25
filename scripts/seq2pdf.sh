#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -euo pipefail
set +o posix

function usage {
  echo "usage ./seq2pdf.sh [--simplify|--dot] <kmer> <file1> [...]"
  echo "  prints pdf to stdout, so please remember to redirect"
  echo "  e.g. ./seq2pdf.sh 5 <(echo ACAACACGT) <(echo CCACACAA) > out.pdf"
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

maxk=$[ ( ($kmer + 31) / 32 ) * 32 - 1 ]

if [[ $[ $kmer & 1 ] -eq 0 || $kmer -lt 3 ]]
then
  echo kmer is not odd and greater than 2
  exit -1
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CTX="$DIR/../bin/ctx$maxk"
CTX2GRAPHVIZ="$DIR/cortex_to_graphviz.pl"
if [[ !(-e $CTX) || !(-x $CTX) ]]
then
  echo "Did you compile for MAXK=$maxk? I cannot run $CTX"
  exit -1
fi

files=$(printf " --seq %s" $@; printf "\n")

if [[ $mkpdf == 1 ]]; then
  $CTX build -q -k $kmer --sample seq2pdf $files - | \
    $CTX2GRAPHVIZ -k $kmer $script_args - | \
    dot -Tpdf
else
  $CTX build -q -k $kmer --sample seq2pdf $files - | \
    $CTX2GRAPHVIZ -k $kmer $script_args -
fi
