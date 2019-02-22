#! /usr/bin/perl -w

use strict;

my $file;
my $align_switches = '-a';
my $n_threads = 0;
for (my $i=0; $i < scalar @ARGV; ++$i) {
    if ($ARGV[$i] eq '-h' || $ARGV[$i] eq '--help') {
	&help();
    }
    elsif ($ARGV[$i] eq '-A' || $ARGV[$i] eq '--aligned') { $align_switches = '-A'; }
    elsif ($ARGV[$i] eq '-P' || $ARGV[$i] eq '--pairwise') { $align_switches = '-P'; }
    elsif ($ARGV[$i] eq '-T' || $ARGV[$i] eq '--threads') {
	if ($i+1 < scalar @ARGV && $ARGV[$i+1] =~ /^[0-9]/) {
	    $n_threads = $ARGV[++$i];
	}
	else { die "-T/--threads require a number as next argument.\n"; }
    }
    elsif ($i == scalar @ARGV-1) { $file = $ARGV[$i]; }
    else { die "Unrecognized argument: $ARGV[$i].\n"; }
}

if (!(-e $file)) { die "Could not find file named $file.\n"; exit; }

if ($align_switches eq '-a') {
    print STDERR "Aligning sequences\n";
    my $outname = $file;
    while (-e $outname) { $outname .= ".aln"; }
    if ($n_threads > 0) { $align_switches = "-T $n_threads"; }
    system "pairalign -a -n -v $align_switches $file > $outname";
    $file = $outname;
    $align_switches .= " --format pairfst -A";
    print STDERR "Aligned sequences\n";
}
elsif ($align_switches eq '-A') {
    if ($n_threads > 0) { $align_switches = "-T $n_threads"; }
    $align_switches .= " --format fasta -A";
}
elsif ($align_switches eq '-P') {
    if ($n_threads > 0) { $align_switches = "-T $n_threads"; }
    $align_switches .= ' --format pairfst -A';
}


print STDERR "Making guide tree to make groups\n";
my $njtree = `pairalign -n -m -j $align_switches $file | treeator -n`;
#print $njtree;
my $n_taxa = `echo '$njtree' | treebender -n`;
chomp $n_taxa;
if ($n_taxa && $n_taxa > 2) {
    print STDERR "Made guide tree with $n_taxa tips.\n";
    &look_for_alignment_groups($njtree,$align_switches,$file);
}
else { print "Less than three taxa, no need for multiple alignment.\n"; }

sub look_for_alignment_groups {
    my $tree = shift;
    my $aligned_switches = shift;
    my $alignment_file = shift;
    my @taxa = `echo '$tree' | treebender -t '\n'`;
    my $pairwise = 0;
    if ($aligned_switches =~ /pairfst/) { $pairwise = 1; }
    chomp @taxa;
    my $pair;
    my $name1;
    my $name2;
    open my $ALIGNMENT, '>', "tempXXXalignment.fst" or die "Could not open tempXXXalignment.fst: $!.\n" ;
    open ALIGN, '<', $alignment_file or die "Could not open file: $alignment_file!\n";
    #print STDERR "Processing taxa in:\n", $tree;
    print STDERR "Parsing relavent sequences from alignment.\n";
    while (<ALIGN>) {
	chomp;
	if (/^>(.*)/) {
	    if (defined($name1) && (defined($name2) || !$pairwise)) {
		if (&in_array(\$name1,\@taxa) && (!$pairwise || &in_array(\$name2,\@taxa))) {
		    print $ALIGNMENT $pair;
		}
		undef $name1;
		undef $name2;
		$pair='';
	    }
	    if (defined $name1) { $name2 = $1; }
	    else {$name1 = $1; }
	    $pair .= '>' . $1;
	    $pair .= "|all\n";
	}
	else { $pair .= $_ . "\n"; }
    }
    if (defined($name1) && (defined($name2) || !$pairwise)) {
    	if (&in_array(\$name1,\@taxa) && (!$pairwise || &in_array(\$name2,\@taxa))) {
	    print $ALIGNMENT $pair;
	}
    }
    close $ALIGNMENT;
    print STDERR "Checking if group of ", scalar @taxa, " taxa is alignable.\n";
    my @alignment_groups = `pairalign $aligned_switches --group alignment_groups tempXXXalignment.fst`;
    my $alignable = 'no';
    foreach (@alignment_groups) {
	if (/Alignment groups: [^_]+_A/) { $alignable = 'yes'; last; }
    }
    if ($alignable eq 'yes') {
	foreach (@taxa) { print $_, ' '; }
	print "\n";
    }
    else {
	print STDERR "Not alignable, spliting tree an checking sub groups.\n";
	my @trees = `echo '$tree' | treebender --null_short_branches 0 | treebender --split m --split_stop t:2`;
	#print @trees;
	if (scalar @trees > 1) {
	    foreach (@trees) {
		my $n_taxa;
		if (/^\(/) {
		    $n_taxa = `echo '$_' | treebender -n`;
		}
		else { $n_taxa = 1; }
		if ($n_taxa && $n_taxa > 1) {
		    &look_for_alignment_groups($_,$aligned_switches,$alignment_file);
		}
		else {
		    #print STDERR "Fewer than 3 taxa in resulting group.\n";
		    #print "#";
		    #if ($n_taxa > 1) {
		#	system "echo '$_' | treebender -t";
		#	#print "\n";
		#    }
		#    else {
			s/;$//;
			print "$_";
		 #   }
		}
	    }
	}
	else {
	    print STDERR "Failed to split tree.\n";
	    print '# ';
	    system "echo '$tree' | treebender -t";
	}
    }
}

sub in_array{
    my $value_ref = shift;
    my $array_ref = shift;
    foreach (@$array_ref) { if ($_ eq $$value_ref) { return 1; } }
    return 0;
}
sub help {
    print "miniptera is a script to find the largest alignabel group acording to MAD score\nusing the phylommand package. Phylommand needs to be in the users PATH. The\ninput should be in fasta format and given in a file as last argument\n\n";
    print "miniptera.pl [-A/-h] file.fst\n\n";
    print "--aligned / -A    will tell miniptera that the sequences are already aligned.\n";
    print "--help / -h       will print this help.\n";
    print "--pairwise / -P   will tell miniptera that the sequences are in pairwise fasta\n";
    print "                  and aligned.\n";
    print "--threads / -T    will give the number of threads to use. Require that pairalign\n";
    print "                  was compiled to use threads.\n";
    exit;
}
