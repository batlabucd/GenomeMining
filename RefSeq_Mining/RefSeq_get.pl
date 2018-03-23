use LWP::Simple;
$name=$ARGV[0];
if($name!~m/[A-Za-z]+/ || $name=~m/[\s\n]+/){
    die "please use alphanumeric or directory-friendly characters in your folder name. No spaces allowed!\n";
}
else{
    `mkdir $name`;
}
$logfile=$name."_Log.txt";
open(LOG, ">>$logfile");
#######################################################################
#Code to loop through all genomes and pull out all matching gene names 
#######################################################################
print "Now looping through all genes and genomes....\n\n";
open(G, "genes_to_get");#open in list of genes
@genes=<G>;#read into memory
open(IN, "genomespec");#open in list of genomes
while(<IN>){#loop through genomes
    if($_=~m/([A-Za-z0-9\s]+)\:(.*)/ && $_!~m/\#/){
        $spec=$1;
        $file=$2;
    }
    print $spec."\n";#print species
    $spec=~s/ /\_/g;
    $ofile="Target_".$spec."_Genes.txt";#make outfile
    $fail="Failed_".$ofile.".fa";#make fail file
    open(FAIL, ">>$fail");
    print $ofile."\n";
    open(OUT, ">>$ofile");
    open(GENOME, "$file");
    @array=(<GENOME>);
    @array=split(/\>/,join('',@array));
    foreach $gene(@genes){
        $genecount=0;
        chomp $gene;
        foreach $seq(@array){
            $seq2=uc $seq;
            if($seq2=~m/\($gene\)/){
                print OUT ">".$seq;
                $genecount++;
            }
        }
        if($genecount==0){
            print FAIL $gene."\n";
        }
    }
    close OUT;
    close FAIL;
}
#######################################################################
#Code to read gi nums and get the CDS from genbank entries
#######################################################################
print "Acquiring CDS sequences for each species....\n\n";
@list=(<Target*txt>);#list of species gene files
foreach $file(@list){#loop through each file
    $speciesdata="";#variable to contain all data
    $ofile="Final_".$file;
    print LOG $ofile."\n";
    open(OUT, ">>$ofile");
    @gi=();#declar array to hold gi numnbers
    if($file=~m/Failed/){
    }
    else{
	print $file."\n";#print filename              
        open(IN, "$file");#open in gene file
        $specname="";#variable to hold species name
	if($file=~m/Target\_([\S]+\_[\S]+)\.txt/){#get species name form gene file
	    $specname=$1;
	}
	while(<IN>){#loop through file line by line
	    if($_=~m/\>/){#if its fasta header
		push @gi, $_;#push header into gi array
	    }
	}
	$i=0;#set i as 0
	while($i<scalar(@gi)){#so long as i is smaller than array length
	    @slice=@gi[$i..$i+99];#get slice of array thats 100 entries long
	    $idstring="";#declare id string, blank it
	    foreach $x(@slice){#for each header in the 100 long slice
		$tag="";#declare tag blank
		if($x=~m/\>([\S]+)/){#get gi number, assuming its the first load of text after the greater sign
		    $tag=$1;#tag is gi header
		}
		if($tag=~m/[\S]+/){#as long as tag exists, add to id string. Tage may not exitss and i may exceednumber of actual entries leading to ,,,,
		    $idstring.=$tag."\,";
		}
	    }
	    $i+=100;#add 100 to i for next slice 100 long
	    $idstring=~s/\,$//;#remove final comma from id string
	    $url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='.$idstring.'&rettype=gb&retmode=text';#URL for getting seqs for ids off ncbi nucleotide
	    print LOG $url."\n";
	    $genbankdata="";
	    $genbankdata=get($url);#merge genbank data into one giant string
	    $failflag=0;
	    if($genbankdata!~m/SOURCE/ || $genbankdata!~m/[\S]+/){
		print "failed to get url. Trying again, will try a max of 3 times\n";
		while($failflag<=3){
		    $genbankdata=get($url);
		    sleep 4;
		    if($genbankdata=~m/[\S]+/){
			$failflag=4;
			print "Failed, but we got it!\n";
		    }
		    else{
			print "Failed again on attempt ".$failflag." out of 3\n\n";
			$failflag++;
		    }
		}
	    }
	    $speciesdata.=$genbankdata;
	    sleep 4;
	}
	print LOG $file."\n\n".$speciesdata."\n";
	@toparse=split(/\/\//,$speciesdata);#split string on each genbank entry
	for($j=0;$j<scalar(@toparse);$j++){#loop through gb entry
	    $dataentry=$toparse[$j];#store data
	    $genename=$seq=$product=$cds="";$start=$end=0;#declare shit load of vars
	    if($dataentry=~m/gene\=\"([\S]+)\"\n/){#get genename
		$genename=$1;
		if($dataentry=~m/CDS[\s]+[\<]*([0-9]+)\.\.([0-9]+)[\>]*/){#get coordinates for start, end
		    $start=$1;$end=$2;
		    if($dataentry=~m/ORIGIN([\S\n\s]+)/){#get dna seq
			$seq=uc $1;
			$seq=~s/\s//g;
			$seq=~s/[0-9]*//g;#remove all the shite
			if($start>0){#so long as the strat coordinate isnt zero
			    print OUT ">".$specname."_(".$genename.")_location=0..0\n";#make species header in form of >[Species name]_[(genename)]
			    $s1=$start-1,$e1=$end+1;#set start as one less for substring function
			    $len=$e1-$start;#length is end cooridniate+1 - start coordinate for substring
			    $cds=substr($seq,$s1,$len);#cds = substring
			    print OUT $cds."\n";
			}
		    }
		}
	    }
	}
    }
    close OUT;
    close IN;
}
close LOG;
#######################################################################
#Code to get the longest canonical transcript for each file
#######################################################################
print "Figuring out longest canonical transcript for each gene file...\n\n";
@target=(<Final_*Genes.txt>);#get list of all final gene files
foreach $file(@target){#foreach one
    open(IN, "genes_to_get");#open original list of genes to get
    @genes=<IN>;#read into array
    close IN;
    $new=$file;#variable
    $new=~s/\s/\_/g;#
    foreach $gene(@genes){#loop through each gene
        $i=0;#counter
        chomp $gene;#remove return char
        @trans=[];#make array
        $tr=0;#
        open(IN, "$file");#open in final gene file
        @array=<IN>;#read data into memory
        @array=split(/\>/,join('',@array));#split on gene
        foreach $x(@array){#foreach gene
            $x=~s/gene\=//;
            $x=">".$x;
            $x=~s/\[/\(/g;
            $x=~s/\]/\)/g;
            if($x=~m/\($gene\)/ || $x=~m/gene\=$gene\]/){#if sequence matches the gene
                $trans[$tr]=$x;#store in memory
                $tr++;
                $i++;
            }
            else{
                $t=$x;#check if perhaps it is the gene, but theres an issue with upper/lower case (carry over from an old script, may not be needed)
                while($t=~s/\(([A-Za-z0-9]+)\)//){
                    $r=$1;
                    $ugene=uc $r;
                    if($ugene eq $gene){
                        $x=~s/\($r\)/\($gene\)/;
                        $trans[$tr]=$x;
                        $tr++;
                        $i++;
                    }
                }
            }
        }
        $long=0;
        $longest="";
        foreach $y(@trans){#loop through each sequence matching the target gene
            if($y=~m/location\=[0-9]+\.\.[0-9]+([\S\s\n]+)/){#get actual gene sequence
                $seq=$1;
                $seq=~s/\s//g;#remove whitspace
                $len=length($seq);#get seq length
                if($len>$long){#compare lengths 
                    $long=$len;#if seq is longer
                    $longest=$y;#make it the curent longest
                }
            }
            else{
            }
        }
        $species="";
        if($file=~m/Final\_([0-9A-Za-z\_]+)\_Genes.txt/){
            $species=$1;
        }
        $ofile=$species."_longest_transcript\.fa";
        open(OUT, ">>$ofile");
        print OUT $longest;
        close OUT;
    }
}
#######################################################################
#Code to generate individual gene files from the longest transcripts
#######################################################################
@transcripts=(<*longest_transcript.fa>);#array for longest transcript files
open(FILE, "genes_to_get");#open list of genes
print "If only one gene found across all species, no output will be generated\n\n";
while(<FILE>){#loop through genes
    $targene=$_;#get gene name
    chomp $targene;
    $gcount=0;
    $data="";
    print "Finding all occurances of ".$targene."\n";
    $out="Final_gene_".$targene."\.txt";#make output
    foreach $trfile(@transcripts){#fore each longest transcript file
	open(IN2, "$trfile");#open file
	@fasta=(<IN2>);#read seqs into array
	@fasta=split(/\>/,join('',@fasta));#split array into genes
	foreach $z(@fasta){#loop through all genes
	    if($z=~m/\($targene\)/){#if gene seq matches target gene
		$data.=">".$z;#take it
		$gcount++;
	    }
	}
    }
    if($gcount>1){#if more than 1 instance of the gene was found, print them out
	open(OUTPUT, ">>$out");
	print OUTPUT $data;
	close OUTPUT;
	$data="";
    }
    else{
	print "either no or only 1 instance of ".$targene." was found, so no final file has been created for it\n\n";
	$data="";
    }
}
close FILE;
#######################################################################
print "Run now complete, moving files to folder: ".$name."/ \n";
`mv *Genes.txt $name/`;
`mv *longest_transcript.fa $name/`;
`mv Final*txt $name/`;
`mv Failed* $name/`;
exit;
