Data and code to reproduce experiments

In order to run any of these experiments, you must run the follow commands to
compile McCortex and download the required data:

    cd ..
    for k in 31 63 95 127; do make MAXK=$k; done
    cd results/data
    ./download.sh

You must also fetch and download the PhiX data yourself from Illumina's
basespace. Details are in results/data/PhiX/about.txt
