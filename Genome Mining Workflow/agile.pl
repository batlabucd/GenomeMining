#!/usr/bin/perl
use strict;
use warnings;
my $datestring = localtime();
print "Time start: ".$datestring."\n";

my %var=&controls();
if($var{'check'} eq "false"){
    die "Error in control file, one or more variable in incorrect format or missing\n";
}
print "Parameter file read successfully\n";
my $filename=$var{'runname'}."Contigs";
my $orig=$var{'genomefile'};
my $sister_taxa=$var{'sistertaxa'};
my $specname=$var{'runname'};
my $finbin=$filename."_files";
my $overhang=$var{'overhang'};my $shearlen=$var{'shearlength'};
my $genome="Sheared_".$shearlen.$filename."\.txt";
my $augmodel=$var{'aug_model'};

############################################################## Code to shear genome file ###########################################
if(!-e $genome){
    print "Shearing genome into smaller contigs\n";
    open(SHEAR, "$orig");open(SHEAROUT, ">>$genome");
    my $concount=0;
    my $toshear="";
    while(<SHEAR>){
	$_=~s/\s//g;
	if($_=~m/\>/){
	    if($concount==0 || $toshear!~m/[A-Z]+/){
		$concount++;
	    }
	    elsif(length($toshear)<$shearlen){
		print SHEAROUT ">Contig".$concount."\n".$toshear."\n";
		$concount++;$toshear="";
	    }
	    else{
		for(my $i=1;$i<length($toshear);$i+=($shearlen-$overhang)){
		    my $offset=$i-1;
		    my $subcontig = substr($toshear, $offset,$shearlen);
		    print SHEAROUT ">Contig".$concount."\n".$subcontig."\n";
		    $concount++;
		}
		$toshear="";
	    }
	}
	else{
	    $toshear.= $_;
	}
    }
    close SHEAR;close SHEAROUT;
}
else{
    print $genome." exists already, genome was not sheared\n";
}

########################################## make blast formatted genome #########################################
if(! -e $specname."\.nsq"){
    print "Formatting blast file\n";
    my $makedb=$var{'makeblastdbpath'};
    `$makedb/makeblastdb -in $genome -dbtype=nucl -out $specname`;#uncomment for blastrun
}
else{
    print $specname."\.nsq Blast databse already exists, skipping makeblastdb step!\n";
}

################################# Loop through all query genes and blast them against target species ##############################
open(IN, "$sister_taxa");
my @blastarray=<IN>;
@blastarray=split(/\>/,join('',@blastarray));
for(my $i=1;$i<scalar(@blastarray);$i++){
    my $filename=$filename."_".$i."_Results.txt";
    open(OUT, ">seq.fa") or die "Problem with initial blast. Most likely an error in the query file\n";
    print OUT ">".$blastarray[$i];
    if(! -e $filename){#if the blast result does not already exist
	my $blastpath=$var{'blastpath'};
	`$blastpath/blastn  -db $specname -query seq.fa -out $filename`;#run blast
    }
    else{
	print $filename." already exists, not running blast!\n";
    }
    
    close OUT;
}
close IN;


################################# Loop through all blast output files, look for paralogs and mask non-functional ones ##############################
print "blast completed\nCycling through contigs for maker\n";
my $num="";
my $handle=my $ohandle="";
my @array=(<$filename*>);
my @targetgenes;
my $cname="";
my $prot="";
my $genetag;
my $gn;my $ty;my $ok;
my $contigs;
my $chimers;
foreach my $blasthit(@array){#loop through each blast file
    if($blasthit=~m/\_([0-9]+)\_/){
	my $num=$1;
	$handle=$filename."_".$num;
        $ohandle=$handle."_out.txt";
    }
    print $blasthit."\n";
    $chimers=$specname."_possible_chimers";
    my @cons;
    my $targetgenesfiles=$var{'target_genes'};
    open(TARS,"$targetgenesfiles");
    @targetgenes=<TARS>;#read in target genes to find gene name for potentially dodgy sequence;
    close TARS;
    chomp @targetgenes;
    open(IN, "$blasthit") or die "Cannot find $blasthit\n\n";#open in file
    my @lines=<IN>;#read file in 
    close IN;
    $ty=0;#variable for counting problematic hits
    my $qlen=0;#length of query
    my $qline="";
    @lines=split(/\>/,join('',@lines));#code to add contig name to blast hits
    for(my $k=1;$k<scalar(@lines);$k++){
	if($lines[$k]=~m/(Contig[0-9]+)/){
	    $cname=$1;
	    $lines[$k]=~s/Score/Score***$cname/;
	}
    }
    my $joined=join('',@lines);
    @lines=split(/Score/,$joined);#split file on each score (ie hit)  
    if($lines[0]=~m/Query\=([\S\s\n]+)Length\=/){#get query gene name                                                                                 
	$prot= $1."\n";
    }
    while($prot=~s/\(([A-Za-z0-9]+)//){#get query gene name
	$gn=$1;
	if($gn~~@targetgenes){
	    $genetag=$gn;
	}
    }
    if($lines[0]=~m/Length\=([0-9]+)/){#get length of query 
	$qlen=$1;
    }
    my $queryseq=""; my $tmp=""; my $revcomp="";
    for(my $i=1;$i<scalar(@lines);$i++){#loop through all hits
	my $query=$lines[$i];
	my $starti=0; my $endi=0; my $refstarti=0;my $refendi=0;#4 variables to describe beginning and end of query and hit regions
	$ok=0;
	my $tmpquery=$query;
	$queryseq="";
	my $w1=0;my $w2=0;
	if($lines[$i]=~m/\*(Contig[0-9]+)/){#get contig name for first hit 
	    $contigs.="\n".$1;
	}
	if($query=~m/Query[\s]+([0-9]+)/){#get start position fo query 
	    $starti=$1;
	}
	if($query=~m/Sbjct[\s]+([0-9]+)/){#get start position of hit
	    $refstarti=$1;
	}
	while($tmpquery=~s/(Query[\s]+[0-9]+[\s]+([A-Z\-a-z]+)[\s]+[0-9]+\n)//){#take the actual gene sequence of query region
	    $qline=$1;
	}
	while($tmpquery=~s/(Sbjct[\s]+[0-9]+[\s]+([A-Z\-a-z]+)[\s]+([0-9]+)\n)//){#take gene sequence of hit region 
	    $queryseq.=$2;
	    $refendi=$3;
	}
	if($qline=~m/([0-9]+)\n/){
	    $endi=$1;# get end position
	}
	if($endi<$starti){#if theyre reverse compliment
	    $contigs.="_rev ";#take note     
	    $tmp=$starti;#swap start and end position numbers
	    $starti=$endi;
	    $endi=$tmp;
	    $revcomp = reverse($queryseq);#also reverse compliment the sequence to forward direction (to look for stop codons) 
	    $queryseq =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
	}
	else{
	    $contigs.=" ";
	}
	if($refendi<$refstarti){
	    $tmp=$refstarti;
	    $refstarti=$refendi;
	    $refendi=$tmp;
	}
	#For each hit, cycle through all other hits, and look if another hit covers the same region of the query that the first hit covers, essentially indicating that they are possible paralogs
	$w1=$endi-$starti;#get width of region covered by hit
	my $previousjcontig="";#scalar for previous contig
	my $revcom;my $prevrev; my $sline;my $totwidth;my $diff; my $range;
	for(my $j=$i+1;$j<scalar(@lines);$j++){#loop through all the hits one ahead of i
	    my $otherq=$lines[$j];#take second hit for overlap comparison  
	    if($otherq=~m/\*(Contig[0-9]+)/){#get contig name
		$previousjcontig=$1;
	    }
	    my $startj=0;my $endj=0;my $refstartj=0;my $refendj=0;#variables for start and end of sequences in the i+1 hit 
	    my $b=$otherq;my $hitseq;
	    $hitseq="";#all the same as in the i loop
	    if($otherq=~m/Query[\s]+([0-9]+)/){
		$startj=$1;
	    }
	    if($otherq=~m/Sbjct[\s]+([0-9]+)/){
		$refstartj=$1;
	    }
	    while($b=~s/(Query[\s]+[0-9]+[\s]+([A-Z\-a-z]+)[\s]+[0-9]+\n)//){
		$sline=$1;
	    }
	    while($b=~s/(Sbjct[\s]+[0-9]+[\s]+([A-Z\-a-z]+)[\s]+([0-9]+)\n)//){
		$hitseq.=$2;
		$refendj=$3;
	    }
	    if($sline=~m/([0-9]+)\n/){
		$endj=$1;
	    }
	    if($endj<$startj){
		$prevrev="_rev ";
		$tmp=$startj;
		$startj=$endj;
		$endj=$tmp;
		$revcom= reverse($hitseq);
		$hitseq=~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
	    }
	    else{
		$prevrev=" ";
	    }
	    if($refendj<$refstartj){
		$tmp=$refstartj;
		$refstartj=$refendj;
		$refendj=$tmp;
	    }
	    $w2=$endj-$startj;#get width of j sequence
	    my $min=0;my $max=0;
	    if($endi eq $qlen){#if the last position happens to match the length of the query sequence, then you have a match up to the end of sequence
		if($queryseq=~m/TAA$/ || $queryseq=~m/TGA/ || $queryseq=~m/TAG/){#remove the possible stop codon, to avoid false posiitve pseudogenes
		    $queryseq=~s/[A-Z][A-Z][A-Z]$//;
		}
	    }
	    if($endj eq $qlen){#same for the j sequence
		if($hitseq=~m/TAA$/ || $hitseq=~m/TGA/ || $hitseq=~m/TAG/){
		    $hitseq=~s/[A-Z][A-Z][A-Z]$//;
		}
	    }
	    if($starti<$startj && $starti<$endj && $starti<$endi){#get the max and min start and end position for the 2 hits
		if($starti!=$startj){
		    $min=$starti;
		}
		else{
		    $min=$starti;
		}
	    }
	    elsif($startj<$starti && $startj<$endj && $startj<$endi){
		if($starti!=$startj){
		    $min=$startj;
		}
		
		else{
		    $min=$startj;##mod this text
		}
	    }
	    if($endi<$endj){
		$max=$endj;
	    }
	    elsif($endj<$endi){
		$max=$endi;
	    }
	    elsif($endj==$endi){
		$max=$endi;
	    }
	    if($min==0){
	    }
	    else{
		$totwidth=$w1+$w2;#get total width of 2 regions, i and j  
		$range=$max-$min;#get the range
		$diff=0;
		my $qstatus="";my $hstatus="";
		my $blank;my $blank2;
		$diff=sqrt(($totwidth-$range)*($totwidth-$range));#make the difference between them positive in instances of when theyre not 
		my $overthresh=$var{'paralogminoverlap'};
		if($range<$totwidth && $diff>$overthresh){#if overlapping, but more so than a threshold
		    $queryseq=~s/[\s\-]+//g;
		     $hitseq=~s/[\s\-]+//g;
		    if($starti<$startj && $endi>$endj || $startj<$starti && $endj>$endi){#if one sequence is entirely contained in another 
			if(length($queryseq)>length($hitseq)){
			    if($queryseq=~m/[A-Z]+$hitseq[A-Z]*/ || $queryseq=~m/[A-Z]*$hitseq[A-Z]+/){
		     #if your single query regions is hitting multiple areas, but those multiple areas have identical sequences, then just let them pass, there is nothing that can be done:
			    }
			    else{
				$qstatus=&checkforstop($queryseq);#send query sequence to look for stop codons in frame 
				if($qstatus eq "TRUE"){# if no frame in query hit region has a functional CDS based on presence of stop codon 
				    $blank="<=".$refstarti."...".$refendi."=> ";#make note of region for masking
				    unless($blank~~@cons){
					$contigs.="<=".$refstarti."...".$refendi."=> ";
					push @cons, $blank;
				    }
				}
				$hstatus=&checkforstop($hitseq);#do same for the i+1 seq
				if($hstatus eq "TRUE"){
				    $blank2="<=".$refstartj."...".$refendj."=> ";
				    unless($blank2~~@cons){
					$tmp=$previousjcontig.$prevrev;
					unless($tmp~~@cons){
					    $contigs.="\n ".$previousjcontig.$prevrev." ";
					    push @cons, $tmp;
					}
					$contigs.="<=".$refstartj."...".$refendj."=> ";
					push @cons, $blank2;
				    }
				}
				if($qstatus ne "FALSE" && $hstatus ne "FALSE"){
				    $ty++;
				}
			    }
			}
			elsif(length($hitseq)>length($queryseq)){
			    if($hitseq=~m/[A-Z]+$queryseq[A-Z]*/ || $hitseq=~m/[A-Z]*$queryseq[A-Z]+/){
			    }
			    else{
				$qstatus=&checkforstop($queryseq);
				if($qstatus eq "TRUE"){
				    $blank="<=".$refstarti."...".$refendi."=> ";
				    unless($blank~~@cons){
					$contigs.="<=".$refstarti."...".$refendi."=> ";
					push @cons, $blank;
				    }
				}
				$hstatus=&checkforstop($hitseq);
				if($hstatus eq "TRUE"){
				    $blank2="<=".$refstartj."...".$refendj."=> ";
				    unless($blank2~~@cons){
					$tmp=$previousjcontig.$prevrev;
					unless($tmp~~@cons){
					    $contigs.="\n ".$previousjcontig.$prevrev." ";
					    push @cons, $tmp;
					}
					$contigs.="<=".$refstartj."...".$refendj."=> ";
					push @cons, $blank2;
				    }
				}
				if($qstatus ne "FALSE" && $hstatus ne "FALSE"){
				    $ty++;
				}
			    }
			}
			elsif(length($hitseq) eq length($queryseq)){
			    if($hitseq eq $queryseq){
			    }
			    else{
				$qstatus=&checkforstop($queryseq);
				if($qstatus eq "TRUE"){
				    $blank="<=".$refstarti."...".$refendi."=> ";
				    unless($blank~~@cons){
					$contigs.="<=".$refstarti."...".$refendi."=> ";
					push @cons, $blank;
				    }
				}
				$hstatus=&checkforstop($hitseq);
				if($hstatus eq "TRUE"){
				    $blank2="<=".$refstartj."...".$refendj."=> ";
				    unless($blank2~~@cons){
					$tmp=$previousjcontig.$prevrev;
					unless($tmp~~@cons){
					    $contigs.="\n ".$previousjcontig.$prevrev." ";
					    push @cons, $tmp;
					}
					$contigs.="<=".$refstartj."...".$refendj."=> ";
					push @cons, $blank2;
				    }
				}
				if($qstatus ne "FALSE" && $hstatus ne "FALSE"){
				    $ty++;
				}
			    }
                        }
                    }#end of if one contained in the other
		    else{# sequences overlap, but one not contained in other
			$qstatus=&checkforstop($queryseq);
			if($qstatus eq "TRUE"){
			    $blank="<=".$refstarti."...".$refendi."=> ";
			    unless($blank~~@cons){
				$contigs.="<=".$refstarti."...".$refendi."=> ";
				push @cons, $blank;
			    }
			}
			$hstatus=&checkforstop($hitseq);
			if($hstatus eq "TRUE"){
			    $blank2="<=".$refstartj."...".$refendj."=> ";
			    unless($blank2~~@cons){
				$tmp=$previousjcontig.$prevrev;
				unless($tmp~~@cons){
				    $contigs.="\n ".$previousjcontig.$prevrev." ";
				    push @cons, $tmp;
				}
				$contigs.="<=".$refstartj."...".$refendj."=> ";
				push @cons, $blank2;
			    }
			}
			if($qstatus ne "FALSE" && $hstatus ne "FALSE"){
			    $ty++;
			}
		    }
		}#end of if overlap and >treshold   
	    }
	}
    }
    my $overall=0;
    if($ty>0){
	$overall++;
	open(CHIM, ">>$chimers");
	print CHIM $genetag."\n";
	close CHIM;
    }
    $targetgenesfiles=$var{'target_genes'};
    open(R,"$targetgenesfiles");#open list of genes
    my @goi=<R>;#take genes of interest into array
    chomp @goi;
    close R;
    my @genelist=[];my $gcount++;
    open(IN, "$blasthit");#open in results file 
    my @hits=<IN>;#read file in        
    close IN;
    @hits=split(/\>/,join('',@hits));#split result file based on contig
    my $allcons=$hits[0];#look at start of blast results
    @lines=split(/\n/,$allcons);#split it on each line
    my $gname;my $name;my $score;my $e;my $revtrue;my $hitseq;
    foreach my $line(@lines){
	my $tl=$line;
	while($tl=~s/\(([A-Za-z0-9]+)\)//){#if the gene matches a gene of interest
	    $gname=$1;
	    if($gname~~@goi){
		#chomp $gname;
		$handle=$handle."_".$gname."\_out.txt";
	    }
	}
	if($line=~m/(Contig[0-9]+)[\s]+([0-9\.]+)[\s]+([\-0-9\.e]+)/){#get contig name
	    $name=$1;$score=$2;$e=$3;#variables for contig name, score and e value
	    $revtrue=0;#check for reverse compliment
	    $hitseq="";
	    if($e<"1e-5" && $score>=40){#min score and evalue, maybe mess with this
		for(my $j=1;$j<scalar(@hits);$j++){#loop through all blast hits in file
		    if($hits[$j]=~m/[\s]*$name[\s\n]+/){#find the matching contig
			$hitseq=$hits[$j];#tmp var for contig
			$j=scalar(@hits)+1;#stop looking fo rcontig
			my $f=0;my $r=0;
			if($hitseq=~m/Sbjct[\s]+([0-9]+)[\s]+[A-Za-z]+[\s]+([0-9]+)/){#check start and end of first line ofhit, see if its reverse comp
			    $f=$1;$r=$2;
			}
			if($f>$r){
			    $revtrue=1;#if end is bigger than start, it is rev
			}
		    }
		}
	    }
	    my $threshevalue=$var{'eval'};my $threshscore=$var{'score'};
	    if($e<$threshevalue && $score>=$threshscore){
		if($name~~@genelist || $name."\_rev"~~@genelist){#if contig has already been counted
		}
		else{#otherwise add it to list
		    if($revtrue==1){
			$genelist[$gcount]=$name."\_rev";
			$gcount++;
		    }
		    else{
			$genelist[$gcount]=$name;
			$gcount++;
		    }
		}
	    }
	}
    }
    open(OUT, ">$handle");
    foreach my $v(@genelist){
	if($v=~m/ARRAY/){}
	else{
	    print OUT $v."\n";
	}
    }
    close OUT;
#Get gene sequences for the contigs!
    my @mask=[];
    if($contigs=~m/\</){
	@mask=split(/\n/,$contigs);
    }
    else{
	@mask=[];
    }
    open(LIST, "$handle");#open list of contigs to get
    my @contigsarray=<LIST>;#read into array
    close LIST;
    open(IN, "$genome");#open genome file
    my $i=0;
    my $final="ConcatSeq_".$handle;#output file
    open(OUT, ">$final");#open output final file
    print OUT ">Concatenated sequence\n";#Print out fastaheader for supercontig
    my $tarseq="";my $a="";my $b="";
    while(<IN>){#read genom file
	if($_=~m/\>/){#if the line has fasta header
	    if($i==2 && $tarseq=~m/[\S]+/){
	    }
	    if($tarseq=~m/[\S]+/ && $i ne "0"){
		$tarseq=~s/\s//g;
		my $mcon="";my $startrem; my $stoprem;my $newstring;my $cover;
		foreach my $masked(@mask){
		    if($masked=~m/(Contig[0-9\_rev]+)/){
			$mcon=$1."\n";
			if($mcon eq $a || $mcon eq $b){
			    while($masked=~s/\<\=([0-9]+)\.\.\.([0-9]+)\=\>//){
				$startrem=$1;$stoprem=$2;
				$cover=$stoprem-$startrem;
				$newstring="";
				for(my $n=0;$n<$cover;$n++){
				    $newstring.="N";#Mask regions representing non-functional paralog
				}
				if(length($tarseq) > $startrem && length($tarseq) > $stoprem && length($tarseq)>0){
				    substr($tarseq,$startrem,$cover)=$newstring;
				}
			    }
			}
		    }
		}
		print OUT $tarseq;
		$tarseq="";
	    }
	    $i=0;#switch is 0
	    $a="";$b="";$tarseq="";
	    if($_=~m/(Contig[0-9]+)[\s\n]/){#if it has contig name
		$a=$1."\n";$b=$1."\_rev\n";#check if forward or reverse in list file
		if($a~~@contigsarray){
		    $i=1;#if its on the list, switch is on
		}
		elsif($b~~@contigsarray){
		    $i=2;#if its on the list, switch is on for rev comp
		}
	    }
	}
	else{#if not a header, then is dna seqeunce
	    if($i==1 || $i==2){#if switch is on
		$tarseq.=$_;
	    }
	}
    }
    close IN;
    `rm $handle`;
    close OUT;

################################# Run Maker ##############################
    open(IN, "$sister_taxa");
    my @query=(<IN>);#take in genes
    @query=split(/\>/,join('',@query));
    my @genes=(<ConcatSeq*>);#list concatenated blast hits
    my $gene;my $tardir;my$logfile;my $ofile;my $genann; my $seqann;my $qhmm;my $qhmm2;my $uniann;my $unidna;my $aug_gene;my $codon;my $protein="";
    foreach my $file(@genes){#for each concat file
	$gene="";
	if($file=~m/\_([A-Za-z0-9]+)\_out/){#get name of gene that was query from filename
	    $gene=$1;
	    foreach my $y(@query){#foreach sister taxa gene
		$protein="";
		if($y=~m/\($gene\)/){#if gene matches one in name
		    $name=$specname."_".$gene;
		    open(OUT, ">tempgene");#open out temp file
		    print OUT ">".$y;#print out gene to tempfile 
		    close OUT;
		    `rm tempgene.prot`;
		    open(OUT, ">tempgene.prot");
		    my $head="";my $curseq=""; 
		    if($y=~m/(.*)\n([A-Za-z\s\n]+)/){
			$head=">".$1."\n";$curseq=$2;
		    }
		    $curseq =~ s/\s//g;
		    for(my $i=0;$i<(length($curseq)-2);$i+=3){
			$codon=substr($curseq,$i,3);
			$protein.=&codons($codon);
		    }
		    print OUT "$head$protein\n";
		    close OUT;
		    open(OPT, "<maker_opts.ctl");#open in maker opts file
		    open(OPOUT, ">pre.ctl");#make new opt file 
		    while(<OPT>){
			if($_=~m/(^genome=)/){
			    print OPOUT $1.$file."\n";
			}
			elsif($_=~m/(^est\=)/){
			    print OPOUT $1."tempgene\n";
			}
			elsif($_=~m/(^protein\=)/){
			    print OPOUT $1."tempgene\.prot\n";
			}
			elsif($_=~m/(^augustus\_species\=)/){
			    print OPOUT $1.$augmodel."\n";
			}
			else{
			    print OPOUT $_;
			}
		    }
		    close OPOUT;
		    close OPT;
		}#finished making maker opt file
	    }
	    my $makerpath=$var{'makerpath'};
	    $tardir=$name."\.maker\.output";
	    if(! -d $name."\.maker\.output"){
		`$makerpath/maker  -base $name -fix_nucleotides -qq  maker_exe.ctl pre.ctl maker_bopts.ctl`;
	    ###################################################### Train SNAP HMM #################################################
		$tardir=$name."\.maker\.output";
		$logfile=$name."_master_datastore_index.log";
		$ofile=$name."_gff";
		`$makerpath/gff3_merge -d $tardir/$logfile -o $tardir/$ofile`;
		`$makerpath/maker2zff -n $tardir/$ofile`;
		`mv genome.dna $tardir/`;
		`mv genome.ann $tardir/`;
		$genann=$tardir."/genome.ann";
		$seqann=$tardir."/genome.dna";
		$qhmm=$tardir."/".$gene."\.hmm";
		`$makerpath/../exe/snap/fathom -categorize 1000 $genann $seqann`;
		`$makerpath/../exe/snap/fathom -export 1000 -plus uni.ann uni.dna`;
		`$makerpath/../exe/snap/forge export.ann export.dna`;
		`$makerpath/../exe/snap/hmm-assembler.pl $gene . >$qhmm`;
		`mv *count $tardir/`;`mv *model $tardir/`;`mv *duration $tardir/`;`mv transitions $tardir/`;`mv phaseprefs $tardir/`;`mv export* $tardir/`;`mv *.ann $tardir/`;`mv *.dna $tardir/`;
		open(READ, "pre.ctl");
		open(OUT4, ">newctl");
		while(<READ>){
		    if($_=~m/(^snaphmm\=)/){
			print OUT4 $1.$qhmm."\n";
		    }
		    elsif($_=~m/(^est2genome\=)/){
			print OUT4 $1."0\n";
		    }
		    elsif($_=~m/(^protein2genome\=)/){
			print OUT4 $1."0\n";
		    }
		    else{
			print OUT4 $_;
		    }
		}
		close READ; close OUT4;
		`mv newctl pre.ctl`;
		`$makerpath/maker -base $name -fix_nucleotides -qq  maker_exe.ctl pre.ctl maker_bopts.ctl`;
		##################################################### End of first HMM snap training ###################################
		###################################################### Train second round of SNAP HMM ################################################# 
		`$makerpath/gff3_merge -d $tardir/$logfile -o $tardir/$ofile`;
		`$makerpath/maker2zff -n $tardir/$ofile`;
		`mv genome.dna $tardir/`;
		`mv genome.ann $tardir/`;
		$genann=$tardir."/genome.ann";
		$seqann=$tardir."/genome.dna";
		$qhmm2=$tardir."/".$gene."2\.hmm";
		`$makerpath/../exe/snap/fathom -categorize 1000 $genann $seqann`;
		`$makerpath/../exe/snap/fathom -export 1000 -plus uni.ann uni.dna`;
		`$makerpath/../exe/snap/forge export.ann export.dna`;
		`$makerpath/../exe/snap/hmm-assembler.pl $gene . >$qhmm2`;
		`mv *count $tardir/`;`mv *model $tardir/`;`mv *duration $tardir/`;`mv transitions $tardir/`;`mv phaseprefs $tardir/`;`mv export* $tardir/`;`mv *.ann $tardir/`;`mv *.dna $tardir/`;
		$uniann=$tardir."/uni.ann";
		$unidna=$tardir."/uni.dna";
		$aug_gene=$tardir."/aug_genbank";
		open(READ, "pre.ctl");
		open(OUT4, ">newctl");
		while(<READ>){
		    if($_=~m/(^snaphmm\=)/){
			print OUT4 $1.$qhmm2."\n";
		    }
		    else{
			print OUT4 $_;
		    }
		}
		close READ; close OUT4;
		`mv newctl pre.ctl`;
		`$makerpath/maker -base $name -fix_nucleotides -qq  maker_exe.ctl pre.ctl maker_bopts.ctl`;
	    }
	    else{
		print $tardir." exists already, so Maker did not run\n";
	    }
	
	    `rm $file`;
    
	    my $mpiblast=$tardir."/mpi_blastdb/";
	    my $voiddir=$tardir."/".$specname."\_".$gene."\_datastore/9F/1D/Concatenated/theVoid.Concatenated/";
	    if(-d $voiddir || -e $voiddir){
		`rm -r $voiddir`;#remove voi directory to save space                                                                                                                
	    }
	    if(-d $mpiblast || -e $mpiblast){
		`rm -r $mpiblast`;#remove to save space                                                                                                                        
	    }
	}
    }
}



 ###################################################### Run EMBOSS pairwise aligner #################################################
#`perl run_emboss.pl $sister $specname`;
my $targetgenesfiles=$var{'target_genes'};
open(IN, "$targetgenesfiles");
my @tar=<IN>;
close IN;
chomp @tar;
my $refseq=$sister_taxa;
my $qspec=$specname;
my $targene;
open(IN, "$refseq");
my @seqs=<IN>;
@seqs=split(/\>/,join('',@seqs));
my $length;my $afile="";my $bfile="";my $ncount;
my $finseq;
foreach my $seq(@seqs){
    $finseq="";
    $seq=">".$seq;
    my $a=$seq;
    my $h;my $tar2="xxx";
    while($a=~s/\(([A-Za-z0-9]+)\)//){
        $targene=$1;
        $tar2=$targene;
        $length=0;
    }
    if($tar2~~@tar){
        open(OUT, ">seqa");
        print OUT $seq;
        close OUT;
        if($seq=~m/\n([A-Z\n]+)/){
            $h=$1;
            $h=~s/\s//g;
            $length=length($h);
        }
        $bfile=$qspec."_".$targene."\.maker\.output/".$qspec."\_".$targene."\_datastore/9F/1D/Concatenated/Concatenated.maker.transcripts.fasta";
       my $backupfile=$qspec."_".$targene."\.maker\.output/".$qspec."\_".$targene."\_datastore/9F/1D/Concatenated/originalout.txt";
	
	if(-e $bfile){
	    open(U, "$bfile");
	    $ncount=0;
	    while(<U>){
		if($_=~m/\>/){
		    $ncount++;
		}
	    }
	    if($ncount==1){
	    }
	    else{
		print "multiple hits found\n";
		open(IN, "$bfile");
		open(OUT, ">tmpfile1");
		
		my @tmparray=<IN>;
		@tmparray=split(/\>/,join('',@tmparray));
		my $best="";
		my $aed=1;my $caed;
		foreach my $y(@tmparray){#If multiple outputs, take the one with best score, as only should be one gene present (other may be due to contigs having genes you dont want)
		    if($y=~m/AED\:[\s]*([0-9\.]+)/){
			$caed=$1;
			if($caed<$aed){
			    $best=$y;
			    $aed=$caed;
			}
		    }
		}
		print  OUT ">".$best;
		close IN;
		close OUT;
		`cp $bfile $backupfile`;
		`mv tmpfile1 $bfile`;
	    }
	}
	my $mpiblast=$qspec."_".$targene."\.maker\.output/mpi_blastdb/";
	my $voiddir=$qspec."_".$targene."\.maker\.output/".$qspec."\_".$targene."\_datastore/9F/1D/Concatenated/theVoid.Concatenated/";
	if(-d $voiddir || -e $voiddir){
	`rm -r $voiddir`;#remove void directory to save space
	}
	if(-d $mpiblast || -e $mpiblast){
	`rm -r $mpiblast`;#remove to save space
	}
	my $outfile=$qspec."_".$targene."_embossmatcher\.aln";
	my $mlen=0;
	my @checkarray;
	if(-e $bfile){
	    open(MAKER, "$bfile");
	    my @maker=<MAKER>;
	    my $maker_seq="";
	    $maker_seq=join('',@maker);
	    my $newseq="";
	    if($maker_seq=~m/\n([A-Z\n]+)/){
		$newseq=$1;
	    }
	    $newseq=~s/\s//g;
	    $mlen=length($newseq);
	    open(CHECKDONE, "$bfile");#open it into array
	    @checkarray=<CHECKDONE>;
	    close CHECKDONE;
	}
	if(-e $bfile && scalar(@checkarray)>1 && ! -e $outfile){#as long as the file exists, and that it has a full sequence
	    my $matcherpath=$var{'matcherpath'};
	    `$matcherpath/matcher -asequence seqa -bsequence $bfile  -outfile $outfile`;#run the aligner delet this
	}
	open(K, ">>$outfile");
        print K "\nRefseq length is ".$length."\n";
        print K "\nMakerseq length is ".$mlen."\n\n";
        close K;
        print $outfile."\n";
        open(ALN, "$outfile");
        my $cutseq="";
	my $refalnseq="";
        my $gaps=0;
        my $identity=0;
        my $alnlen=0;
        my $matching=0;
        my $gapperc="";
        my $final=$qspec."_Final_Genes.fa";
        my $headerforseq="XXXX";my $percentage_of_maker_kept;my $ref_coverage;
        open(FINAL, ">>$final");
        while(<ALN>){
            if($_=~m/maker\-[\s]+([A-Z\-]+)/){
                $cutseq.=$1;
            }
            elsif($_=~m/august[\s]+([A-Z\-]+)/){
                $cutseq.=$1;
            }
            elsif($_=~m/snap[\S]+[\s]+([A-Z\-]+)/){
                $cutseq.=$1;
            }
            elsif($_=~m/Identity\:[\s]+([0-9]+)\/([0-9]+)[\s]+\(([0-9\.]+)\%/){
                $identity=$3;
		$alnlen=$2;
                $matching=$1;
            }
            elsif($_=~m/Gaps\:[\s]+([0-9]+)\/[0-9]+[\s]+\(([0-9\.]+)/){
                $gaps=$1;
                $gapperc=$2;
            }
            elsif($_=~m/1\:[\s]+([\S][\S][\S])/){
                $headerforseq=$1;
            }
            elsif($_=~m/^$headerforseq[\S]+[\s]+([A-Z\-]+)/ || $_=~m/^[\s]+$headerforseq[\S]+[\s]+([A-Z\-]+)/ ){
                $refalnseq.=$1;
	    }
        }
	$refalnseq=~s/[\s\-]+//g;
        
        $cutseq=~s/\-//g;
	$percentage_of_maker_kept=0;
	$ref_coverage=0;
        unless($length==0 || $mlen==0){
            $percentage_of_maker_kept=sprintf("%.2f",(length($cutseq)/$mlen)*100);#metric
	    $ref_coverage=sprintf("%.2f",(length($refalnseq)/$length)*100);#changed form length(h) to just $length
        }
	my $refthresh=$var{'coverage'};
        if($ref_coverage>$refthresh){
	   my $check=0;my $aa1=my$aa2=my $aa3=0;my $output_header;
	   my $testseq=$cutseq;
	   my $f2=substr($cutseq, 1, length($cutseq)-1);
	   my $f3=substr($cutseq, 2, length($cutseq)-1);
	   open(OUT1, ">>t1");
	   print OUT1 ">Seq1\n".$testseq."\n".">Seq2\n".$f2."\n".">Seq3\n".$f3."\n";
	   close OUT1;
	   my @aa;my $codon;my $protein;
	   for(my $p=0;$p<(length($testseq)-2);$p+=3){
	       $codon=substr($testseq, $p,3);
	       $protein.=&codons($codon);
	   }
	   $aa[0]=">Seq1\n".$protein."\n";
	   $protein="";
	   for(my $p=0;$p<(length($f2)-2);$p+=3){
	       $codon=substr($f2, $p,3);
	       $protein.=&codons($codon);
           }
           $aa[1]=">Seq2\n".$protein."\n";
	   $protein="";
	   for(my $p=0;$p<(length($f3)-2);$p+=3){
	       $codon=substr($f3, $p,3);
               $protein.=&codons($codon);
           }
           $aa[2]=">Seq3\n".$protein."\n";
	   $aa1=0;$aa2=0;$aa3=0;
	   $output_header=">".$qspec."_".$targene."_IdentwithRefspec_".$identity."_CoverageofRef_".$ref_coverage."_percentage_of_original_maker_gene_".$percentage_of_maker_kept."\n";
	   while($aa[0]=~s/\*//){$aa1++};while($aa[1]=~s/\*//){$aa2++};while($aa[2]=~s/\*//){$aa3++};
	   $finseq="";
	   if($aa1<$aa2 && $aa1<$aa3){
	       $finseq=$output_header.$testseq."\n";
	   }
	   elsif($aa2<$aa1 && $aa2<$aa3){
	       $finseq=$output_header.$f2."\n";
	   }
	   elsif($aa3<$aa1 && $aa3<$aa2){
	       $finseq=$output_header.$f3."\n";
	   }
	   my $fbpath=$var{'framebotpath'};
	   `java -jar $fbpath/FrameBot.jar index seqa protindex`;
	   open(NUC, ">nuc");
	   print NUC $finseq;
	   close NUC;
	   `java -jar $fbpath/FrameBot.jar framebot -o frame.txt -l 20  protindex nuc`;
	   open(IN3, "frame.txt_corr_nucl.fasta");
	   my $framefixed="";
	   while(<IN3>){
	       $framefixed.=$_;
	   }
	   close IN3;
	   if($framefixed=~m/[A-Za-z\S]+/){
	       print FINAL $framefixed;
	       $framefixed="";#modified this
	   }
	   else{
	       print FINAL $finseq;
	       $finseq="";###modified this
	   }
	}
    }
    close FINAL;
    if(-e 'seqa'){
    `rm seqa`;
    }
}
`mkdir $finbin`;
`mv $filename* $finbin/`;
my $final_maker_out=$specname."_Maker_Run";
`mkdir $final_maker_out`;
`mv $specname* $final_maker_out/`;
$datestring = localtime();
print "Time finished: ".$datestring."\n";



sub checkforstop{
    my $qseq=$_[0];
    my $f1=$qseq;
    my $ret="UNKNOWN";
    my $f2=substr($qseq,1,length($qseq));
    my $f3=substr($qseq,2,length($qseq));
    my $f1count=0;my $f2count=0;my $f3count=0;
    my $codon;
    for(my $i=0;$i<length($f1);$i++){
	$codon=substr($f1,$i,3);
	if($codon eq  "TGA" || $codon eq "TAA" || $codon eq "TGA"){
	    $f1count++;
	}
	$i=$i+2;
    }
    for(my $i=0;$i<length($f2);$i++){
        $codon=substr($f2,$i,3);
        if($codon eq  "TGA" || $codon eq "TAA" || $codon eq "TGA"){
            $f2count++;
        }
        $i=$i+2;
    }
    
    for(my $i=0;$i<length($f3);$i++){
        $codon=substr($f3,$i,3);
        if($codon eq  "TGA" || $codon eq "TAA" || $codon eq "TGA"){
            $f3count++;
        }
        $i=$i+2;
    }
    if($f1count >0 && $f2count>0 && $f3count>0){
	$ret="TRUE";
    }
    elsif($f1count==0 || $f2count==0 || $f3count==0){
	$ret="FALSE";
    }
    return $ret;
}

sub codons{
    my $codon=$_[0];
    $codon = uc $codon;
    my (%nucs)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G','NNN'=>'X','GCN'=>'A','RAY'=>'B','TGY'=>'C','GAY'=>'D','GAR'=>'E','TTY'=>'F','GGN'=>'G','CAY'=>'H','ATH'=>'I','AAR'=>'K','TTR'=>'L','CTN'=>'L','YTR'=>'L','AAY'=>'N','CCN'=>'P','CAR'=>'Q','CGN'=>'R','AGR'=>'R','MGR'=>'R','TCN'=>'S','AGY'=>'S','ACN'=>'T','GUN'=>'V','TAY'=>'Y','SAR'=>'Z','TAR'=>'*','TRA'=>'*');
    if(exists $nucs{$codon}){
	return $nucs{$codon};
    }
    else{
	return "X";
    }
}



sub controls{
    open(IN, "ctrlfile") or die "parameter file not found\n";
    my %hash;
    $hash{'check'}='true';
    while(<IN>){
	if($_=~m/[\S]+/){
	    if($_=~m/Overhang\:\"([0-9]+)\"\;/){
		$hash{'overhang'}=$1;
	    }
	    elsif($_=~m/ShearLength\:\"([0-9]+)\"\;/){
		$hash{'shearlength'}=$1;
	    }
	    elsif($_=~m/GenomeFile\:\"([\S]+)\"\;/){
		$hash{'genomefile'}=$1;
	    }
	    elsif($_=~m/SisterTaxa\:\"([\S]+)\"\;/){
		$hash{'sistertaxa'}=$1;
	    }
	    elsif($_=~m/TestName\:\"([\S]+)\"\;/){
		$hash{'runname'}=$1;
	    }
	    elsif($_=~m/OverlapDiff\:\"([0-9]+)\"\;/){
		$hash{'paralogminoverlap'}=$1;
	    }
	    elsif($_=~m/Evalue\:\"([\S]+)\"\;/){
		$hash{'eval'}=$1;
	    }
	    elsif($_=~m/Score\:\"([\S]+)\"\;/){
		$hash{'score'}=$1;
	    }
	    elsif($_=~m/Ref\_Coverage\:\"([0-9]+)\"\;/){
		$hash{'coverage'}=$1;
	    }
	    elsif($_=~m/Target_genes\:\"([\S]+)\"\;/){
		$hash{'target_genes'}=$1;
	    }
	    elsif($_=~m/Augustus_model\:\"([\S]+)\"\;/){
		$hash{'aug_model'}=$1;
	    }
	    elsif($_=~m/Maker_Path\:\"([\S]+)\"\;/){
		if(-d $1){
		    $hash{'makerpath'}=$1;
		}
		else{
		    die "Error, path not found: \n".$1."\n";
		}
	    }
	    elsif($_=~m/Blastn_Path\:\"([\S]+)\"\;/){
		if(-d $1){
		    $hash{'blastpath'}=$1;
		}
		else{
		    die "Error, path not found: \n".$1."\n";
		}
	    }
	    elsif($_=~m/Makeblastdb_Path\:\"([\S]+)\"\;/){
		if(-d $1){
		    $hash{'makeblastdbpath'}=$1;
		}
		else{
		    die "Error, path not found: \n".$1."\n";
		}
	    }
	    elsif($_=~m/Matcher_Path\:\"([\S]+)\"\;/){
		if(-d $1){
		    $hash{'matcherpath'}=$1;
		}
		else{
		    die "Error, path not found: \n".$1."\n";
		}
	    }
	    elsif($_=~m/Framebot_Path\:\"([\S]+)\"\;/){
		if(-d $1){
		    $hash{'framebotpath'}=$1;
		}
		else{
		    die "Error, path not found: \n".$1."\n";
		}
	    }
	    else{
		print "line is: ".$_;
		$hash{'check'}='false';
	    }
	}
	else{die "parameter file exists but is empty!\n";
	}
    }
    return %hash;
}
    
    
