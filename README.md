## Footprint Detection ##
> Shane Neph


Overview
=========
The core program that implements the DNaseI [footprinting description] from: An expansive human regulatory lexicon encoded in transcription factor footprints, Nature 489, 83–90 (2012).

Given a set of sequencing tag counts (integers) at each base of any region, this program creates an unthresholded list of candidate footprints.


Build
=====
// requires g++ version 4.7 or newer  
make -C src/


Program Options:
================
```
detect-cache
	[--flankmin <bases>  = 6]
	[--flankmax <bases>  = 12]
	[--centermin <bases> = 6]
	[--centermax <bases> = 100]
	[--maxthold <value>  = 10]
	File-Full-O-Integers
```

The footprint occupancy statistic (FOS) of a candidate footprint is defined as:
FOS = (C+1)/L + (C+1)/R, where C is the average number of tags over the central/core region of a footprint, while L (R) is the average tag level found in the left (right) flanking region.

--flankmin/--flankmax set the minimum/maximum number of bases to consider when calculating L and R around a potential footprint.
--centermin/--centermax set the minimum/maximum size (# bp) of all candidate footprints.
--maxthold should be ignored and use the default value of 10


Input
=====
This program accepts a file full of integers that represent the number of sequencing tags at each base over the region of interest, including zeroes where no tags are found at a particular base.

Regions of interest can be easily broken up into sequences of integers and run in parallel on a cluster or multiple-core system.  Breaking up your data by chromosome may be the most natural method to run things in parallel.


Results
=======
The output of this program consists of unthresholded candidate footprints.  The file format has 9 columns.

Column #:
---------
1. The leftmost position of the left-flanking region
2. The start of the central/core potential footprint (1 bp beyond the end of left-flanking region)
3. The start of the right-flanking region (1 bp beyond the end of the core region)
4. The end of the right-flanking region (1 bp beyond the end of the right-flanking region)
5. The FOS
6. The mean of the left-flanking region (L)
7. The mean of the core region (C)
8. The mean of the right-flanking region (R)

Note that the program does not know about chromosomes.  Further, it reads in and interprets the first integer as belonging to absolute position 0.  So, if feed it something that does not start at base 0, you need to adjust the output using the input base offset.

Even if you paste on chromosome information to the beginning of the output, be careful as the output is not in [sorted BED order].

Typically, one would threshold the potential footprints based upon some metric that utilizes the FOS, paste on appropriate chromosome information as the first field, rearrange columns to convert the output to a BED-like format, and then sort to obtain the final result.

[footprinting description]: http://www.nature.com/nature/journal/v489/n7414/extref/nature11212-s1.pdf
[sorted BED order: https://bedops.readthedocs.org/en/latest/content/reference/file-management/sorting/sort-bed.html]
