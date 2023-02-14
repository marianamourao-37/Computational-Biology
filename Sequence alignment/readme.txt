There are two functions: smithwaterman and path_find. The latter is an auxiliary to the first, which is the main program and the one you should run.

Provide three inputs:
s1 and s2 - The sequences you wish to align. Write them in char format, ex: 'WPIWPC'.
d - The linear gap penalty.

You will obtain two outputs:
score - The score of the optimal local alignment(s).
list_opt - A matrix in which each column is an optimal local alignment. The first row corresponds to the first sequence while the second row corresponds to the second sequence.


EXAMPLE:

If you introduce the following in your command window:

	[a,b] = smithwaterman('WPIWPC','IIWPI',4)

You will get an output:

	a = 30

	b = 2x2 string array
	"WPI"	"IWP"
	"WPI"	"IWP"
 