Isaac Turner
2014 Sep 24

In order to re-run the experiment, from the mccortex directory run:

    # Fetch libraries needed
    cd libs && make core common && cd..
    # Compile McCortex
    make MAXK=31
    # Generate the reference from chr22
    cd results/data/chr22/uniq_flanks && make && cd ../../../..
    # Run the experiment
    cd results/exp_abc/1MbpHg19
    make

That's it!

On my macbook run time is 40 mins
