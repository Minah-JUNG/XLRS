#!/usr/bin/perl
if($#ARGV < 0){
print "perl script input_file\n";
} else {

        $file = $ARGV[0];
        open(FILE,"<$file");
        open(OUTPUT,">Top_Selected_$file");

        @array = <FILE>;
        $limit = @array;

        $key0 = "";

        for($i=0; $i < $limit; $i++){

        $line = $array[$i]; chomp $line;

        ($key, $contents) = split(/\t/,$line,2);

        if($key ne $key0){

        $key0 = $key;

        print OUTPUT "$line\n";

        } else {

        }

        }

}
close OUTPUT;
close FILE;
exit;
