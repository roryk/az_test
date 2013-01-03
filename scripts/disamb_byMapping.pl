#!/usr/bin/perl 

# To compare to SAM files from explant sequencing
# keep singletons

my $sam1 = shift;   #human sam file
my $sam2 = shift;   #mouse sam file

my %pair;  #paried reads

##############
#process human mapping
##############

print STDERR "processing $sam1 human mapping: ", `date`;

open(SAM1, $sam1) or die "can not open $sam1\n";

my %sam1;
keys(%sam1) = 7000000;

my $n = 0;

while(defined($line=<SAM1>)) {
    if ($line !~ /^@/) {
	my @a = split(/\t/, $line);
	$n++;

	$pair{$a[0]} = $a[0];

	my $dir = $a[1] & 0x40 ? "1" : "2";
	my $read = "$a[0]_$dir";  #readID and direction (1 or 2)

	my $nm = $1 if ( /NM:i:(\d+)/ );
	my $nh = $1 if ( /NH:i:(\d+)/ );
	my $gp = /XO:i:(\d+)/ ? $1 : 0;
	my $score = $nm + $nh + $gp;

	if ($score > 20) {
	    print "human score: $score\n";
	}

	my $mate = $a[1] & 0x8 ? 0 : 1; # Mate is unmapped

	if (!$sam1{$read}) {
	    $sam1{$read} = {S => $score, Q => $a[4], M => $mate};
	}

	# overwrite if new mapping is better
	elsif ($sam1{$read} && $sam1{$read}->{S} > $score ) {
	    $sam1{$read} = {S => $score, Q => $a[4], M => $mate};
	}

	print STDERR "$n mapping processed in human mapping for $sam1\n" if ( $n%10000000 == 0 );
    }
}

my $num = keys(%sam1);
print STDERR "$n total mapping processed in human mapping, $num of unique reads mapped\n";
close(SAM1);

##############
#process mouse mapping
##############

print STDERR "processing $sam2 mouse mapping: ", `date`;;

open(SAM2, $sam2) or die "can not open $sam2\n";


my %sam2;
keys(%sam2) = 7000000;

my $n = 0;

while (defined($line=<SAM2>)) {
    if ($line !~ /^@/) {
	my @a = split(/\t/, $line);
	$n++;

	$pair{$a[0]} = $a[0];

	my $dir = $a[1] & 0x40 ? "1" : "2";
	my $read = "$a[0]_$dir";  #readID and direction (1 or 2)

	my $nm = $1 if ( /NM:i:(\d+)/ );
	my $nh = $1 if ( /NH:i:(\d+)/ );
	my $gp = /XO:i:(\d+)/ ? $1 : 0;
	my $score = $nm + $nh + $gp;

	if ($score > 20) {
	    print "mouse score: $score\n";
	}

	my $mate = $a[1] & 0x8 ? 0 : 1; # Mate is unmapped

	if (!$sam2{$read}) {
	    $sam2{$read} = {S => $score, Q => $a[4], M => $mate};
	}

	# overwrite if new mapping is better
	elsif ($sam2{$read} && $sam2{$read}->{S} > $score ) {
	    $sam2{$read} = {S => $score, Q => $a[4], M => $mate};
	}

	print STDERR "$n mapping processed in mouse mapping for $sam2\n" if ( $n%10000000 == 0 );
    }
}

my $num = keys(%sam2);
print STDERR "$n total mapping processed in mouse mapping, $num of unique reads mapped\n";
close(SAM2);

########################
#differentiate
########################

open (HUMAN, ">human.disambreads.txt");
open (MOUSE, ">mouse.disambreads.txt");
open (AMB, "ambigouou.disambreads.txt");

open (OUT, ">map_class_pe4.txt");

foreach my $key (keys %pair) {
    my $t1 = "$key"."_1";
    my $t2 = "$key"."_2";  #mate

    my ($h1, $h2) = ($sam1{$t1}->{S} ? $sam1{$t1}->{S} : 1000, $sam1{$t2} ? $sam1{$t2}->{S} : 1000);  #human mapping score, set it to 1000 if unmapped
    my ($m1, $m2) = ($sam2{$t1}->{S} ? $sam2{$t1}->{S} : 1000, $sam2{$t2} ? $sam2{$t2}->{S} : 1000);  #mouse mapping score, set it to 1000 if unmapped

    if (($h1==$m1)&&($h2==$m2)) {
	$identity{$key} = "amb";
	print OUT "$key amb $h1==$h2==$m1==$m2\n";
	print AMB "$key\n";
	$amb++;
    }

    elsif (($h1==$m2)&&($h2==$m1)) {
	$identity{$key} = "amb";
	print OUT "$key amb $h1==$h2==$m1==$m2\n";
	print AMB "$key\n";
    }

    else {
	my @score = sort {$a->[0] <=> $b->[0] } ([$h1,"human"], [$h2, "human"], [$m1, "mouse"], [$m2, "mouse"]);
	
	if ($score[0]->[1] eq $score[1]->[1]) {  #top 2 are from same species
	    $identity{$key} = $score[0]->[1];

	    if ($score[0]->[1] eq "human") {
		print HUMAN "$key\n";
	    }
	    elsif ($score[0]->[1] eq "mouse") {
		print MOUSE "$key\n";
	    }

	    print OUT "$key $score[0]->[1] by first 2 reads $h1==$h2==$m1==$m2\n";
	}

	else {                #top 2 are from different species 
	    if ($score[0]->[0] < $score[1]->[0]) {   #disambiguation by first 2 reads
		$identity{$key} = $score[0]->[1];

		if ($score[0]->[1] eq "human") {
		    print HUMAN "$key\n";
		}
		elsif ($score[0]->[1] eq "mouse") {
		    print MOUSE "$key\n";
		}

		print OUT "$key $score[0]->[1] by unmatched first 2 reads $h1==$h2==$m1==$m2\n";
	    }

	    else {                                   #first 2 reads with same scores, disambiguation by second 2 reads
		$identity{$key} = $score[2]->[1];

		if ($score[2]->[1] eq "human") {
		    print HUMAN "$key\n";
		}
		elsif ($score[2]->[1] eq "mouse") {
		    print MOUSE "$key\n";
		}

		print OUT "$key $score[2]->[1] by unmatched 2nd 2 reads $h1==$h2==$m1==$m2\n";
	    }
	}
    }
}

