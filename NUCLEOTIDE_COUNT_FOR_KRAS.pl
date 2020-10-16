#!/usr/bin/perl
# 



open (INPUT3, "$ARGV[0]_amplicon.txt") or die "\nAmp file could not be found\n\n";
while (<INPUT3>) {
	chomp ($_);
	($amp,$amplicon,$start,$end,$total_length) = split(/\t/,$_);
	$amplicon =~ tr/acgt/ACGT/;
	@AMP = split('',$amplicon);
	$start_term = "$AMP[0]"."$AMP[1]"."$AMP[2]"."$AMP[3]"."$AMP[4]"."$AMP[5]";
	$end_term = "$AMP[-6]"."$AMP[-5]"."$AMP[-4]"."$AMP[-3]"."$AMP[-2]"."$AMP[-1]";
	$results = "$ARGV[0]"."_SLX_14586_C5K9C_read_counts_hash.txt";
	open(OUTPUT, ">>$results");

#	print OUTPUT "\n";
#$index =1;
#	while ($index <= 150){
#		print OUTPUT "$index\t";
#		$index = $index +1;
#	}
#	print OUTPUT "\n";
	
	open (INPUT, "filenames.txt") or die "\nFile names file could not be found\n\n";
	while (<INPUT>) {
		chomp ($_);
		$input_file = $_;
		if ($input_file =~ /.fq/){
			open (INPUT2, "$input_file") or die "\n$_ could not be found\n\n";
			($lib,$fld,$hft) = split(/\./,$input_file);
			
			@TOTAL_BASE =();
			@A_BASE =();
			@C_BASE =();
			@G_BASE =();
			@T_BASE =();
		
			$index = 0;
			while ($index < scalar@AMP){
				push (@TOTAL_BASE, "0");
				push (@A_BASE, "0");
				push (@C_BASE, "0");
				push (@G_BASE, "0");
				push (@T_BASE, "0");
				$index = $index +1;
			}
		
			while (<INPUT2>) {
				chomp ($_);
			
			
			
				
				if (($_ =~ /^$start_term/) && ($_ =~ /$end_term$/) && (length($_) == $total_length)){
					$index = $start +1;
					$correct_amplicons = $correct_amplicons +1;
#					print "$amplicon\n";
					while ($index <= $end +1){
				
						$term = "$AMP[$index-6]"."$AMP[$index-5]"."$AMP[$index-4]"."$AMP[$index-3]"."$AMP[$index-2]"."$AMP[$index-1]".
							"([ACGT])"."$AMP[$index+1]"."$AMP[$index+2]"."$AMP[$index+3]"."$AMP[$index+4]"."$AMP[$index+5]"."$AMP[$index+6]";
							
						
						if ($_ =~ /$term/){
							$TOTAL_BASE[$index] = $TOTAL_BASE[$index] +1;
							if ($1 eq "A"){
						 		$A_BASE[$index] = $A_BASE[$index] +1;
							}
							if ($1 eq "C"){
							 	$C_BASE[$index] = $C_BASE[$index] +1;
							}
							if ($1 eq "G"){
							 	$G_BASE[$index] = $G_BASE[$index] +1;
							}
							if ($1 eq "T"){
							 	$T_BASE[$index] = $T_BASE[$index] +1;
							}
						}
						$index = $index +1;
					}
				}	
			}
		
			close (INPUT2);
#			print OUTPUT "$input_file\n";
#			foreach $val(@AMP){
#				print OUTPUT "$val\t";
#			}
			print OUTPUT "$fld"."_A\t";
			foreach $val(@A_BASE){
				print OUTPUT "$val\t";
			}
			print OUTPUT "\n";
			print OUTPUT "$fld"."_C\t";
			foreach $val(@C_BASE){
				print OUTPUT "$val\t";
			}
			print OUTPUT "\n";
			print OUTPUT "$fld"."_G\t";
			foreach $val(@G_BASE){
				print OUTPUT "$val\t";
			}
			print OUTPUT "\n";
			print OUTPUT "$fld"."_T\t";
			foreach $val(@T_BASE){
				print OUTPUT "$val\t";
			}
			print OUTPUT "\n";
			print OUTPUT "$fld"."_Sum\t";
			foreach $val(@TOTAL_BASE){
				print OUTPUT "$val\t";
			}
			print OUTPUT "\n";
			$index = $start +1;

		}
		print "Done $amp $fld...\n";
	}
	close (INPUT);
	close (OUTPUT);
}
print "\nCorrect amplicons = $correct_amplicons\n";
close (INPUT3);
close (OUTPUT2);
exit;
