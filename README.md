## General Instructions to Install and Run OrthoMCL (v1.4)

### OrthoMCL and MCL Background

OrthoMCL is an algorithm developed by L. Li, C. Stoeckert, and D. Ross specifically for identifying ortholog clusters for eukaryotes. Grouping orthologous genes allows us to infer similar functions across-species and make claims regarding conservation. OrthoMCL aims to identitfy 'recent' paralogs and orthologs via an all-vs-all BLASTP (i.e. BLASTP each proteome within and against each other). These homologs are then converted into a network-like graph where genes are nodes and edges are weighted by normalized BLAST results. MCL (see next section) is then applied to this graph. MCL offers the benefit of "evaluat[ing] the global pattern of sequence similarities among provisionally grouped sequences during clustering". OrthoMCL offered advantages over tools at the time such as INPARANOID (only does pairwise species comparison) and COG (which can be fooled by multidomain proteins as it considers local similarity over global similarity). 

<div style="text-align:center">
	<div>
		<img src="https://genome.cshlp.org/content/13/9/2178/F1.large.jpg" alt="omcl-flowchart" width="250"/>
	</div>
	<i>OrthoMCL pipeline, excerpt from Li et al., (2003)</i>
</div>

<div style="text-align:center">
	<div>
		<img src="https://genome.cshlp.org/content/13/9/2178/F2.large.jpg" alt="omcl-flowchart" width="250"/>
	</div>
	<i>Creating a similairty matrix from BLAST results, excerpt from Li et al., (2003)</i>
</div>


Markov Clustering (MCL) is an unsupervised clustering algorithm based on the simple idea that random walks in a network will typically stay within a cluster. Namely it simulates network flow by using Markov matrices and altering inflation/expansion operators until a steady state matrix is reached (converegence). Where expansion "corresponds to computing random walks of higher length, which means random walks with many steps" and inflation "changes the probabilities associated with the collection of random walks departing from one particular node (corresponding with a matrix column) by favouring more probable walks over less probable walks" (i.e. increasing the probability of intra-cluster walks). Put more simply, we alter between expansion (which encourages intercluster flow) and inflation (discourages intercluster flow) until equilibrium is reached. It is important to note that one can adjust the inflation operator and change cluster tightness as mentioned in the OrthoMCL paper.

- [MCL Documentation](https://micans.org/mcl/)
- [Useful MCL vis] (https://micans.org/mcl/index.html?sec_mcl)
- [OrthoMCL original paper] (https://genome.cshlp.org/content/13/9/2178.full)
 
### Installing OrthoMCL

As this protocol is based on R. Patel's 2012 Expressolog paper which used an older version of OrthoMCL and installing OrthoMCL v2 requires more dependencies such as MySQL, I decided to use OrthoMCL v1.4. For tutorials explaining how to install and use v2, please see [here](https://www.biostars.org/p/199083/) and [here](http://darencard.net/blog/2018-01-12-orthomcl-tutorial/). OrthoMCL simply needs your list of proteomes in fasta format (ending .fa) to work! Here are the list of dependencies I had to install onto the cluster for v1.4:

- OrthoMCL (attached, but found [here](https://orthomcl.org/common/downloads/software/unsupported/v1.4/))
- MCL (unpacked from above)
- BioPerl [Bio::SearchIO](https://metacpan.org/pod/Bio::SearchIO) and [Storable](https://metacpan.org/pod/Storable) Perl modules
- Legacy Blast (attached, but also found [here](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/)) which has blastall

<b> NB: Some of the instructions for installing the above are specific to the cluster I am using so your mileage may vary! Namely, the cluster I am using does not allow users to have root access so perl modules must be installed locally. </b>


1. Make sure you have perl installed and it is version 5+. ```perl --version```
2. Download cpanm, and local::lib to a local dir such as ~. Add appropriate environmental variables. [Credits](https://stackoverflow.com/questions/2980297/how-can-i-use-cpan-as-a-non-root-user):

    ```
    wget -O- http://cpanmin.us | perl - -l ~/perl5 App::cpanminus local::lib
	eval `perl -I ~/perl5/lib/perl5 -Mlocal::lib`
	echo 'eval `perl -I ~/perl5/lib/perl5 -Mlocal::lib`' >> ~/.bash-rc
	echo 'export MANPATH=$HOME/perl5/man:$MANPATH' >> ~/.bash_rc
    ```

3. Install Storable and Bio::SearchIO. ```cpanm Storable``` ```cpanm Bio::SearchIO```
4. Unpack MCL and compile. [Credits](https://github.com/sanger-pathogens/companion/blob/master/test/travis.setup-orthomcl.sh). Load gcc module if necessary. After unpacking:

    ```
	cd mcl-*
	./configure
	make -j3
	make install
    ```
5. Unpack BLAST.     
6. Unpack orthomcl and edit ```orthomcl_module.pm``` with settings below (these instructions are in the README):
	- $PATH\_TO_ORTHOMCL: the orthomcl directory itself (i.e. this directory)
	- $BLASTALL: location of blastall you unpacked earlier
	- $BLAST_NOCPU: number of CPUs, incredibly useful for clusters
	- $FORMATDB: location of formatdb you unpacked earlier
	- $MCL: location of compiled mcl, try shmcl/mcl
	- $ORTHOMCL\_DATA_DIR: location of your fasta files (they should all be in one dir)
	
	Also adjust these variables to our expressolog paper parameters:
	
	- $BLAST\_PVALUE\_CUTOFF_DEFAULT         = 1e-10;
	- $PERCENT\_IDENTITY\_CUTOFF_DEFAULT     = 60;
	- $PERCENT\_MATCH\_CUTOFF_DEFAULT        = 60;
	- $MCL\_INFLATION_DEFAULT               = 2.2;
7. Write your submission script to your cluster, remembering to adjust the # of CPUs and walltime. My run of 3 proteomes (~50k genes each) took ~11h on a 2 node/80 CPU cluster. You may have to write the perl script to tell the cluster to look into your local directory for perl modules. 

	```
	perl -I ~/perl5/lib/perl5 orthomcl.pl --mode 1 --fa_files "genome1.fa,genome2.fa,genome3.fa"
	```

8. Your output should come in a folder with the MON_DD, with a file named ```all_orthomcl.out```. It should look something like this:

	```
	ORTHOMCL4913(4 genes,3 taxa):    AT1G03910.1(TAIR10_pep_20101214_updated) AT1G03910.2(TAIR10_pep_20101214_updated) Os03t0800100-01(IRGSP-1.0_protein_2019-06-26) Solyc10g085410.2.1(ITAG3.2_proteins)
	```

