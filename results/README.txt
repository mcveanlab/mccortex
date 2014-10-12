Data and code to reproduce experiments

In order to run any of these experiments, you must run the follow in the root
directory of McCortex:

    cd libs
    make common core
    cd ..
    for k in 31 63 95 127; do make MAXK=$k; done

Then you need to download the data required:

    cd results/data/
    ./download.sh

You must also fetch and download the PhiX data yourself from Illumina's
basespace. Details are in results/data/PhiX/about.txt
