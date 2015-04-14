Scripts to compare jelly-hash to mccortex internal hash table.

JellyHash is a low memory hash table:
https://github.com/noporpoise/jelly-hash.git

Results are included for linux and mac. 

./benchmark-tables.sh
./times.sh results20150409thurs.linux.log > results20150409thurs.linux.csv
./times.sh results20150409thurs.mac.log > results20150409thurs.mac.csv
./stats.R results20150409thurs.linux.csv > results20150409thurs.linux.txt
./stats.R results20150409thurs.mac.csv > results20150409thurs.mac.txt
grep mean results20150409thurs.linux.txt | awk '{if(!s){print $4"\t"$6}; s=!s; }'
grep mean results20150409thurs.linux.txt | awk '{if( s){print $4"\t"$6}; s=!s; }'
grep mean results20150409thurs.mac.txt   | awk '       {print $4"\t"$6}'


results20150409thurs.linux.log
- run on WTCHG banyan
- Processor: 8 x Intel(R) Xeon(R) 8 core CPU running at 2.70GHz
- Memory: 2TB of RAM
- OS: CentOS release 6.5; linux kernel version 2.6.32
- Compiler:
    $ cc -v
    Using built-in specs.
    Target: x86_64-redhat-linux
    Configured with: ../configure --prefix=/usr --mandir=/usr/share/man --infodir=/usr/share/info --with-bugurl=http://bugzilla.redhat.com/bugzilla --enable-bootstrap --enable-shared --enable-threads=posix --enable-checking=release --with-system-zlib --enable-__cxa_atexit --disable-libunwind-exceptions --enable-gnu-unique-object --enable-languages=c,c++,objc,obj-c++,java,fortran,ada --enable-java-awt=gtk --disable-dssi --with-java-home=/usr/lib/jvm/java-1.5.0-gcj-1.5.0.0/jre --enable-libgcj-multifile --enable-java-maintainer-mode --with-ecj-jar=/usr/share/java/eclipse-ecj.jar --disable-libjava-multilib --with-ppl --with-cloog --with-tune=generic --with-arch_32=i686 --build=x86_64-redhat-linux
    Thread model: posix
    gcc version 4.4.7 20120313 (Red Hat 4.4.7-4) (GCC)
- mccortex commit 95061472a50d8ec882963bafd45ff1a1b3070cd3
- jelly-hash commit 1454199956d93ed626201b594042ea378bfc7136
- command: ~/mccortex/results/hash_table_benchmark$ ./benchmark-tables.sh >& results20150409thurs.linux.log

results20150409thurs.mac.log
- run on Isaac's laptop (macbook pro)
- MacBook Pro (13-inch, Mid 2010)
- Processor: 2.66 GHz Intel Core 2 Duo
- Memory: 4 GB 1067 MHz DDR3
- OS: Mac OS X (10.10.3)
- Compiler:
    $ cc -v
    Apple LLVM version 6.1.0 (clang-602.0.49) (based on LLVM 3.6.0svn)
    Target: x86_64-apple-darwin14.3.0
    Thread model: posix
- mccortex commit 95061472a50d8ec882963bafd45ff1a1b3070cd3
- jelly-hash commit 1454199956d93ed626201b594042ea378bfc7136
- command: ~/mccortex/results/hash_table_benchmark$ ./benchmark-tables.sh --lowmem >& results20150409thurs.mac.log
- Note only includes 80M entry hash tables
