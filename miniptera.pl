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

if ($align_switches ne '-A') {
    print STDERR "Aligning sequences\n";
    my $outname = $file;
    while (-e $outname) { $outname .= ".aln"; }
    if ($n_threads > 0) { $align_switches = " -T $n_threads"; }
    system "pairalign -a -n -v $align_switches $file > $outname";
    $file = $outname;
    $align_switches .= "-A --format pairfst";
    print STDERR "Aligned sequences\n";
}
else {
    if ($n_threads > 0) { $align_switches = " -T $n_threads"; }
    $align_switches .= " --format fasta";
}


print STDERR "Making guide tree for making groups\n";
my $njtree = `pairalign -n -m -j -v $align_switches $file`; # | treebender --split m --split_stop t:2`;
my $n_taxa = `echo $njtree | treebender -n`;
print STDERR "Made guide tree with $n_taxa tips.\n";
if ($n_taxa > 2) {
    &look_for_alignment_groups($njtree,$align_switches,$file);
}
else { print "Less than three taxa, no need for multiple alignment"; }

sub look_for_alignment_groups {
    my $tree = shift;
    my $aligned_switches = shift;
    my $alignment_file = shift;
    my @taxa = `echo $tree | treebender -t '\n'`;
    chomp @taxa;
    my $alignment='';
    open ALIGN, '<', $alignment_file or die "Could not open file: $alignment_file!\n";
    my $print = 'no';
    my $pair_in_group = 'first';
    print STDERR "Parsing relavent sequences from alignment.\n";
    while (<ALIGN>) {
	chomp;
	if (/^>(.*)/) {
	    my $name = $1;
	    my $test = &in_array(\$name,\@taxa);
	    if ($test && ($pair_in_group eq 'first' || $pair_in_group eq 'yes')) { $print = 'yes'; }
	    else { $print = 'no'; }
	    if ($aligned_switches =~ /pairfst/ && $pair_in_group eq 'first' && $test) { $pair_in_group = 'yes'; }
	    elsif ($aligned_switches =~ /pairfst/ && $pair_in_group eq 'first' && !$test) { $pair_in_group = 'no'; }
	    else { $pair_in_group = 'first'; }
	    if ($print eq 'yes') { $alignment .= $_ . "|all\n"; }
	}
	if ($print eq 'yes') { $alignment .= $_ . "\n"; }
    }
    print STDERR "Checking if group is alignable.\n";
    my @alignment_groups = `pairalign $aligned_switches --group alignment_groups`;
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
	my @trees = `echo $tree | treebender --split m --split_stop t:2`;
	foreach (@trees) {
	    my $n_taxa = `echo $_ | treebender -n`;
	    if ($n_taxa && $n_taxa > 3) {
		&look_for_alignment_groups($_,$aligned_switches,$alignment_file);
	    }
	    else {
		print "Loose taxa: ";
		system "echo $_ | treebender -t";
		print "\n";
	    }
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
    print "miniptera is a script to find the largest alignabe group acording to MAD score\nusing the phylommand package. Phylommand needs to be in the users PATH. The\ninput should be in fasta format and given in a file as last argument\n\n";
    print "miniptera.pl [-A/-h] file.fst\n\n";
    print "--aligned / -A    will tell miniptera that the sequences are already aligned.\n";
    print "--help / -h       will print this help.\n";
    exit;
}
