#!/usr/bin/perl 

# To compare to SAM files from explant sequencing
# keep singletons
#sam file name convention:
#map to human genome:  sample_name.human.sam
#map to mouse genome:  sample_name.mouse.sam

my $sam1 = shift;   #input sam file, from mapping to human genome
my $sam2 = shift;   #input sam file, from mapping to mouse genome
my $out_dir = shift;   #output directory for 4 sam files

system("mkdir $out_dir");

@sam1 = split(/\//, $sam1);
$sam1[-1] =~ /(.*)\.human.sam/;
my $name = $1;    #sample name

#output sam files:
open (HUMAN, ">$out_dir/$name.disambiguousHuman.sam");
open (MOUSE, ">$out_dir/$name.disambiguousMouse.sam");
open (AMB1, ">$out_dir/$name.ambiguousHuman.sam");
open (AMB2, ">$out_dir/$name.ambiguousMouse.sam");

my %pair;  #paried reads

##############
#process human mapping
##############

print STDERR "processing $sam1 human mapping: ", `date`;

open(SAM1, $sam1) or die "can not open $sam1\n";

my %sam1;   #readName=>score and the line
keys(%sam1) = 7000000;

my $n = 0;

while(defined($line=<SAM1>)) {
    if ($line =~ /^@/) {   #print header
	print HUMAN "$line";
	print AMB1 "$line";
    }

    else {
	my @a = split(/\t/, $line);
	$n++;

	$pair{$a[0]} = $a[0];

	my $dir = $a[1] & 0x40 ? "1" : "2";
	my $read = "$a[0]_$dir";  #readID and direction (1 or 2)

	my $nm = $1 if ( /NM:i:(\d+)/ );
	my $nh = $1 if ( /NH:i:(\d+)/ );
	my $gp = /XO:i:(\d+)/ ? $1 : 0;
	my $score = $nm + $nh + $gp;

	my $mate = $a[1] & 0x8 ? 0 : 1; # Mate is unmapped

	if (!$sam1{$read}) {
	    $sam1{$read} = {S => $score, Q => $a[4], M => $mate, A => $line};
	}

	# overwrite if new mapping is better, all alignments are kept
	elsif ($sam1{$read} && $sam1{$read}->{S} > $score ) {
	    $sam1{$read} = {S => $score, Q => $a[4], M => $mate, A => "$sam1{$read}->A$line"};
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

my %sam2;     #readName=>score and the line
keys(%sam2) = 7000000;

my $n = 0;

while (defined($line=<SAM2>)) {
    if ($line =~ /^@/) {   #print header
	print MOUSE "$line";
	print AMB2 "$line";
    }
    else {
	my @a = split(/\t/, $line);
	$n++;

	$pair{$a[0]} = $a[0];

	my $dir = $a[1] & 0x40 ? "1" : "2";
	my $read = "$a[0]_$dir";  #readID and direction (1 or 2)

	my $nm = $1 if ( /NM:i:(\d+)/ );
	my $nh = $1 if ( /NH:i:(\d+)/ );
	my $gp = /XO:i:(\d+)/ ? $1 : 0;
	my $score = $nm + $nh + $gp;

	my $mate = $a[1] & 0x8 ? 0 : 1; # Mate is unmapped

	if (!$sam2{$read}) {
	    $sam2{$read} = {S => $score, Q => $a[4], M => $mate, A => $line};
	}

	# overwrite if new mapping is better, all alignments are kept
	elsif ($sam2{$read} && $sam2{$read}->{S} > $score ) {
	    $sam2{$read} = {S => $score, Q => $a[4], M => $mate, A => "$sam2{read}->{A}$line"};
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

my %human;      #reads assigned human
my %mouse;      #reads assigned mouse
my %amb;        #ambiguous reads

foreach my $key (keys %pair) {
    my $t1 = "$key"."_1";
    my $t2 = "$key"."_2";  #mate

    print "$t1==$t2\n";

    my ($h1, $h2) = ($sam1{$t1}->{S} ? $sam1{$t1}->{S} : 1000, $sam1{$t2} ? $sam1{$t2}->{S} : 1000);  #human mapping score, set it to 1000 if unmapped
    my ($m1, $m2) = ($sam2{$t1}->{S} ? $sam2{$t1}->{S} : 1000, $sam2{$t2} ? $sam2{$t2}->{S} : 1000);  #mouse mapping score, set it to 1000 if unmapped

    if (($h1==$m1)&&($h2==$m2)) {
	$amb{$key} = $key;
	print "$key=amb\m\n";
    }

    elsif (($h1==$m2)&&($h2==$m1)) {
	$amb{$key} = $key;
	print "$key=amb\m\n";
    }

    else {
	my @score = sort {$a->[0] <=> $b->[0] } ([$h1,"human"], [$h2, "human"], [$m1, "mouse"], [$m2, "mouse"]);
	
	if ($score[0]->[1] eq $score[1]->[1]) {       #top 2 are from same species

	    if ($score[0]->[1] eq "human") {
		$human{$key} = $key;
	    }
	    elsif ($score[0]->[1] eq "mouse") {
		$mouse{$key} = $key;
	    }
	}

	else {                                       #top 2 are from different species 
	    if ($score[0]->[0] < $score[1]->[0]) {   #disambiguation by first 2 reads
		if ($score[0]->[1] eq "human") {
		    $human{$key} = $key;
		}
		elsif ($score[0]->[1] eq "mouse") {
		    $mouse{$key} = $key;
		}
	    }

	    else {                                   #first 2 reads with same scores, disambiguation by second 2 reads
		if ($score[2]->[1] eq "human") {
		    $human{$key} = $key;
		}
		elsif ($score[2]->[1] eq "mouse") {
		    $mouse{$key} = $key;
		}
	    }
	}
    }
}

#####################
#uptput 4 sam files
#####################

foreach my $key (keys %pair) {
    my $t1 = "$key"."_1";
    my $t2 = "$key"."_2";  #mate

    if ($human{$key}) {
	if ($sam1{$t1}) {
	    print HUMAN "$sam1{$t1}->{A}";
	}
	if ($sam1{$t2}) {
	    print HUMAN "$sam1{$t2}->{A}";
	}
    }
    if ($mouse{$key}) {
	if ($sam2{$t1}) {
	    print MOUSE "$sam2{$t1}->{A}";
	}
	if ($sam2{$t2}) {
	    print MOUSE "$sam1{$t2}->{A}";
	}
    }
    if ($amb{$key}) {
	if ($sam1{$t1}) {
	    print AMB1 "$sam1{$t1}->{A}";
	}
	if ($sam1{$t2}) {
	    print AMB1 "$sam1{$t2}->{A}";
	}
	if ($sam2{$t1}) {
	    print AMB2 "$sam2{$t1}->{A}";
	}
	if ($sam2{$t2}) {
	    print AMB2 "$sam1{$t2}->{A}";
	}
    }
}
    
	
	
