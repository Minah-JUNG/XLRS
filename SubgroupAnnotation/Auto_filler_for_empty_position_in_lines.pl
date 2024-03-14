#!/usr/bin/perl
if($#ARGV < 0){
exit;
} else {

        $file = $ARGV[0];
#       $max = $ARGV[1];

        open(FILE,"<$file");
        open(OUTPUT,">Filled_$file");

        @array = <FILE>;
        $limit = @array;

	@col_number = ();

	for($k=0; $k < $limit; $k++){

	$line_read = $array[$k]; chomp $line_read;

	@arr = split(/\t/,$line_read);

	$ncol = @arr; 

	$col_number[$k] = $ncol;

	}

	@sorted = sort {$b <=> $a} @col_number;
        
        $max = $sorted[0];

        for($i=0; $i < $limit; $i++){

        $line = $array[$i]; chomp $line;
        @element = split(/\t/,$line);

        @push = ();

                for($k=0; $k < $max; $k++){

                if($element[$k] eq ""){

                $push[$k] = "-";

                } else {

                $push[$k] = $element[$k];

                }

                }

        $string = join "\t", @push;

        print OUTPUT "$string\n";

        }

}

close FILE;
close OUTPUT;
