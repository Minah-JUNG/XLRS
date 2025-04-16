#!/usr/bin/perl
if($#ARGV < 1){
print "This script is optimized for bi-allelic filtered VCF file.\n";
exit;
} else {
	$file = $ARGV[0];
	$detail = $ARGV[1];

	open(FILE,"<$file");
	open(OUTPUT,">$detail\_filtered.hapmap");
	open(STATISTICS,">$detail\_SNP_statistics.txt");
	@array = <FILE>;
	$limit = @array;

print STATISTICS "SNP\tGene\tMajor_allele\tMinor_allele\tMAF\tMissing_allele\tNo._of_homozygous_line\tNo._of_heterozygous_line\tNo._of_missing_line\tNo._of_total_line\tHomozygous_rate\tHeterozygous_rate\tMissing_rate\tTotal_proportion\tDominant_alternative\n";
	for($i=0; $i < $limit; $i++){
	$line = $array[$i]; chomp $line;
	@element = split(/\t/,$line); $cols = @element - 1;
	$head = $element[0];
	
	@subset = @element[9..$cols];
	$list_of_acc = join("\t",@subset);
	if($head eq "#CHROM"){
			
	print OUTPUT "rs\talleles\tchrom\tpos\tstrand\tassembly\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t$list_of_acc\n";  $start = $i+1; last;} else { next;}
	}

	for($k=$start; $k < $limit; $k++){
	$line = $array[$k]; chomp $line;
	@element = split(/\t/,$line); $cols = @element - 1;
	$rs = "$element[0]" . "_" . "$element[1]"; $chrom = $element[0];

	$chrom =~ s/HsG0//g;
	$chrom =~ s/HsG//g;

	 $pos = $element[1]; $strand = "."; $assembly = "$detail"; $RefAllele = "$element[3]";

	$alternative = $element[4];
	($af, $as, $at) = split(/,/,$alternative,3);

	 %ALLELE = (0 => $RefAllele,
                     1 => $af,
                     2 => $as,
                     3 => $at);

	@subset = @element[9..$cols]; $lim_s = @subset;

		$homo = 0; $hetero = 0; $miss = 0;

		for($s=0; $s < $lim_s; $s++){
		$allele_info = $subset[$s];

		($all1, $mix) = split(/\//,$allele_info,2);
		$all2 = substr($mix,0,1);
		$print1 = $ALLELE{$all1};
		$print2 = $ALLELE{$all2};

		if($all1 eq ".") {

		$new[$s] = "NN";

		$miss = $miss + 1;

		} else {

		$new[$s] = "$print1" . "$print2";

			if($print1 ne $print2){
			
			$hetero = $hetero + 1;			

			} else {

			$homo = $homo + 1;

			}


		}

		}

	$sequence = join("",@new);
	@SEQ = split(//,$sequence);
	@sorted_new = sort @SEQ;

	$sorted_sequence = join("",@sorted_new);
	
	$limit1 = @SEQ;
	
	$SNPS = join("\t",@new);

	$A=0; $T=0; $G=0; $C=0; $N=0;
	for($j=0; $j < $limit1; $j++){
        $character = $sorted_new[$j];

        if($character eq "A"){
        $A = $A + 1;
        } elsif($character eq "T"){
        $T = $T + 1;
        } elsif($character eq "G"){
        $G = $G + 1;
        } elsif($character eq "C") {
        $C = $C + 1;
        } else {
	$N = $N + 1;
	}

        }

	@count = ($A, $T, $G, $C);

	$missing_allele = $N;
	%matching = (
		     $A => "A",
                     $T => "T",
                     $G => "G",
                     $C => "C");
        @sort_count = sort { $b <=> $a } @count;
        $major = $matching{$sort_count[0]};
        $minor = $matching{$sort_count[1]};


	if($minor eq $RefAllele){

	$tag = "O";

	} else {

	$tag = "X";

	}

	if($sort_count[0] == $sort_count[1]){
	
	$allele = "$RefAllele" . "/" . "$alternative";

	} else {

        $allele =  "$major" . "/" . "$minor";

	}


	$MAF = $sort_count[1] / ($sort_count[0] + $sort_count[1] + 0.000000000000000000000000000000000000000000000000001) * 100;

	$missing_rate = $N / ($A + $T + $G + $C + $N) * 100;	

	$MISSING = $miss / ($homo + $hetero + $miss) * 100;

	$HETEROZYGOUS = $hetero / ($homo + $hetero + $miss) * 100;

	$HOMOZYGOUS = $homo / ($homo + $hetero + $miss) * 100;

	$Num_acc = $homo + $hetero + $miss;

	$TOTAL = $MISSING + $HETEROZYGOUS + $HOMOZYGOUS; 

	print OUTPUT "$rs\t$allele\t$chrom\t$pos\t$strand\t$assembly\tNA\tNA\tNA\tNA\tQC+\t$SNPS\n";

	print STATISTICS "$rs\t$detail\t$sort_count[0]\t$sort_count[1]\t$MAF\t$missing_allele\t$homo\t$hetero\t$miss\t$Num_acc\t$HOMOZYGOUS\t$HETEROZYGOUS\t$MISSING\t$TOTAL\t$tag\n";
		
}


}
close FILE;
close OUTPUT;
close STATISTICS;
