## Footprint Detection ##
> by Shane Neph


Overview
=========
The core program that implements the DNaseI [footprinting description] from: An expansive human regulatory lexicon encoded in transcription factor footprints, Nature 489, 83–90 (2012).

Given a set of sequencing tag counts (integers) at each base of any region, this program creates an unthresholded list of candidate footprints.


Build
=====
// requires g++ version 4.7 or newer  
make -C src/  
This will create subdirectory bin/ and place _fp2012_ in it.


Program Options:
================
```
fp2012
	[--flankmin <bases>  = 6]
	[--flankmax <bases>  = 12]
	[--centermin <bases> = 6]
	[--centermax <bases> = 100]
	[--maxthold <value>  = 10]
	File-Full-O-Integers
```

The footprint occupancy score (FOS) of a candidate footprint is defined as:  
```FOS = (C+1)/L + (C+1)/R,```
where C is the average number of tags over the central/core region of a potential footprint, and L (R) is the average tag level found in the left (right) flanking region.  A lower FOS is a more significant score.  Any case where L or R is less than or equal to C is ignored, and, consequently, no division by zero can occur.

--flankmin  
--flankmax set the min/max number of flanking bases over which to find the max mean value for L and R  

--centermin  
--centermax set the min/max number of bases over which to find the min mean value for C  

--maxthold should be ignored so that the program uses the default value of 10


Input
=====
This program accepts a file full of integers that represent the number of cleavage events assigned to each base.  That is, each base receives an integer that shows the number of uniquely-mapping sequencing tags with 5' ends mapping to that position.  For those bases with none, they receive a zero.  Since cleavage occurs in between two nucleotides, choose the 5'-most of these to represent each cleavage event.  It is not strictly necessary to use only uniquely-mapping tags though it is common in practice.

You may use a dash (-) to denote that input comes from stdin.

Regions of interest can be broken up into subsequences by file.  We recommend that you partition your data by chromosome and run the _fp2012_ program on each of these in parallel.  Note that the [bedextract] and [unstarch] programs are designed to stream data by chromosome very efficiently.  Use the applicable tool to stream the data, and then select the column of interest using the built-in _cut_ command.


Results
=======
The output of this program consists of unthresholded candidate footprints.  The output is deterministic so if you run everything again with the same inputs, you will receive the same output.  The file format has 8 columns.

1. The leftmost position of the left-flanking region
2. The start of the central/core potential footprint (1 bp beyond the end of left-flanking region)
3. The start of the right-flanking region (1 bp beyond the end of the core region)
4. The end of the right-flanking region (1 bp beyond the end of the right-flanking region)
5. The FOS
6. The mean tag level of the left-flanking region (L)
7. The mean tag level of the core region (C)
8. The mean tag level of the right-flanking region (R)

All candidate footprints' core regions are disjoint and they do not abutt.  The set of chosen candidate footprints are optimized relative to each other in cases where they would overlap or abutt, and each candidate is optimized over the input parameter settings.  You will next want to threshold results further, including the removal of candidate footprints with too many unmapped bases in their core regions.  Another common filter is to ignore results outside of DNaseI FDR 1% [hotspot regions].

Note that the program does not know about chromosomes.  Further, it reads in and interprets the first integer as belonging to absolute position zero.  If you feed it something that does not start at base zero, you need to adjust the output using the first input base's offset.

Even if you paste on chromosome information to the beginning of the output, be careful as the output is not in [sorted BED order].

Typically, one thresholds the potential footprints based upon some metric that utilizes the FOS, rearranges columns 2&3 to 1&2, pastes on appropriate chromosome information as the first field, and then [sorts] to obtain the final result.  Using this procedure, you end up with 0-based [start,end) footprint calls where columns 2 and 3 hold the core footprint start and end positions.


Portability
===========
This program runs fine on Linux systems, Mac OS X, and BSD systems.  It is written in standard C++, so it should compile and run natively on Windows systems though a different build manager would need to replace the simple makefile.  The makefile hardcodes g++ as the compiler.  Change CC in the makefile (for example, to clang++) along with build flags as you see fit.


Tips
====
One method of sticking zeroes in for bases that have no per-base number of cuts [using bedops], where column 5 holds the base counts, looks like:
```
  bedops -c -L <bases-with-tag-counts> \
    | awk '{ for (i=$2; i<$3; ++i) { print $1,i,i+1,".",0; } }' \
    | bedops -u - <bases-with-tag-counts> \
    | cut -f5
```

Note that ```<bases-with-tag-counts>``` must be a [properly sorted] BED file, and the output of this command sequence can be piped directly into the _fp2012_ program.  Here, all bases beyond the last found in ```<bases-with-tag-counts>``` will have no integer representation.  This could affect results in only the slightest way (the very last footprint call if several conditions are all met).  You can add another file to the _bedops -c_ call to put zeroes all of the way to the end of a chromosome for completeness if that is a concern.  Ask a question on [our forum] if advice is needed.


Performance and scalability
===========================
We regularly run this program on deeply-sequenced data using a compute cluster.  We break the genome up by chromosome and submit each to a cluster of modest machines.  We typically set flanking search parameters to 3-10 and the center/core footprint search sizes at 6-40.  With this method, you can expect full results in less than one hour with a genome roughly the size of that for human.  One can implement tricks to restrict inputs to less than a whole chromosome (for example, restrict to 1% FDR DNaseI hotspots) in order to speed up computations considerably.  The tradeoff is a significantly larger amount of bookkeeping to create inputs and to glue the final results together properly.  This seems deceptively simple to do, and it's easy to overlook pitfalls.  We recommend running _fp2012_ with data from an entire chromosome.

The program can use a bit of main memory and we recommend 2G or more RAM.  Surprisingly, feeding the program all zeroes gives the worst case memory performance.  That is something that we plan to address in the future.  Note that we have never had any memory issues in practice when using real data sets with 30 million or more uniquely-mapping sequencing tags.


[footprinting description]: http://www.nature.com/nature/journal/v489/n7414/extref/nature11212-s1.pdf
[sorted BED order]: https://bedops.readthedocs.org/en/latest/content/reference/file-management/sorting/sort-bed.html
[sorts]: https://bedops.readthedocs.org/en/latest/content/reference/file-management/sorting/sort-bed.html
[properly sorted]: https://bedops.readthedocs.org/en/latest/content/reference/file-management/sorting/sort-bed.html
[using bedops]: https://bedops.readthedocs.org/en/latest/content/reference/set-operations/bedops.html
[bedextract]: https://bedops.readthedocs.org/en/latest/content/reference/set-operations/bedextract.html
[unstarch]: https://bedops.readthedocs.org/en/latest/content/reference/file-management/compression/unstarch.html
[our forum]: http://bedops.uwencode.org/forum/
[hotspot regions]: https://github.com/StamLab/hotspot
