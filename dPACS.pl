#!/usr/bin/perl

use Docopt;
use List::Util qw( min max );
use IPC::Open2;
my $protofile = '/home/galaxy/galaxy/tools/test/proto.511';    #path to REBASE proto file
my $TMPDIR    = "/tmp/dpc.$$";
$out = `mkdir /tmp/dpc.$$`;

=head1 COPYRIGHT 

Copyright © 2018, Syngenta. All rights reserved.

=head1 LICENSE INFORMATION

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

THIS PROGRAM IS MADE AVAILABLE FOR DISTRIBUTION WITHOUT ANY FORM OF WARRANTY TO THE 
EXTENT PERMITTED BY APPLICABLE LAW.  THE COPYRIGHT HOLDER PROVIDES THE PROGRAM \"AS IS\" 
WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT  
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM LIES
WITH THE USER.  SHOULD THE PROGRAM PROVE DEFECTIVE IN ANY WAY, THE USER ASSUMES THE
COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION. THE COPYRIGHT HOLDER IS NOT 
RESPONSIBLE FOR ANY AMENDMENT, MODIFICATION OR OTHER ENHANCEMENT MADE TO THE PROGRAM 
BY ANY USER WHO REDISTRIBUTES THE PROGRAM SO AMENDED, MODIFIED OR ENHANCED.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL THE 
COPYRIGHT HOLDER BE LIABLE TO ANY USER FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL,
INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE
PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE
OR LOSSES SUSTAINED BY THE USER OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO 
OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER HAS BEEN ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGES.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.


=head1 NAME

dPACS.pl
 
=head1 USAGE

    dPACS.pl --wildseq=<oligo> --mutseq=<oligo> --mutations=<mutation number> --scan=<left|right|straddle|all> [--cut=<all|wild|mut>] [--zero] [--revcomp]
    dPACS.pl -h | --help | --perldoc

    <oligo> = nucleotide sequence
    <mutation number> = number of mutations in restriction site  
    --scan = left|right|straddle|all (restriction site mutations to left or right of snp site, straddling snp site, or anywhere).
    --cut = all|wild|mut (filter for restriction sites which will cut all, only the wildtype or only the mutant, sequence, Default = all)
    --zero = for all enzymes that cut, include --mutations=0
    --revcomp = reverse complement all sequences

    output is to the screen

=head1 DEPENDENCIES
    
    REBASE proto file from http://rebase.neb.com/rebase/link_proto ($protofile)
    EMBOSS configured with version of REBASE used above
    Access to TMP space ($TMPDIR)
 
=cut

$opts = docopt(%args);
foreach my $arg (@ARGV)
    {
        if ( $arg =~ /--wildseq=(\S+)/ )
            {
                $wildseq = $1;
                $wildseq =~ s/["' ]//g;
            }
        if ( $arg =~ /--mutseq=(\S+)/ )
            {
                $mutseq = $1;
                $mutseq =~ s/["' ]//g;
            }
        if ( $arg =~ /--scan=(left|right|straddle|all)/ )
            {
                $scan = $1;
            }
        if ( $arg =~ /--cut=(all|wild|mut)/ )
            {
                $cut = $1;
            }
        if ( $arg =~ /--mutations=(\d+)/ )
            {
                if ( $1 eq '0' )
                    {
                        $mutswant = '.';
                    }
                else
                    {
                        $mutswant = $1;
                    }
            }
        if ( $arg eq "--zero" )
            {
                $zero = '1';
            }
        if ( $arg eq "--perldoc" )
            {
                my $pr = `perldoc $0`;
                print $pr;
                exit;
            }
        if ( $arg eq "--revcomp" )
            {
                $revcomp = '1';
            }

    }

# Create a set of pattern files for fuzznuc to use. Take the rebase PROTO file, reformat it for fuzznuc, and split it by site length.
open( PROTO, $protofile );
while (<PROTO>)
    {
        last if ( $_ =~ /TYPE I ENZYMES|TYPE III ENZYMES/ );
        unless ( $_ =~ /REBASE|Rich|TYPE|=-=|------/ )
            {
                my @inline = split( /\s+/, $_ );
                if ( $#inline == 1 )
                    {
                        $patt = uncut( $inline[1] );
                        $len  = length($patt);
                        ${$len}{ $inline[0] } = $patt;
                        push @list, $len unless $seen{$len}++;
                        $protosite{ $inline[0] } = $inline[1];
                    }
                if ( $#inline == 2 )
                    {
                        $patt = uncut( $inline[1] );
                        $len  = length($patt);
                        ${$len}{ $inline[0] } = $patt;
                        push @list, $len unless $seen{$len}++;
                        $protosite{ $inline[0] } = $inline[1] . " " . $inline[2];
                    }
                if ( $#inline == 3 )
                    {

                        $patt = uncut( $inline[2] );
                        $len  = length($patt);
                        ${$len}{ $inline[0] } = $patt;
                        push @list, $len unless $seen{$len}++;
                        $protosite{ $inline[0] } = $inline[1] . " " . $inline[2] . " " . $inline[3];
                    }

            }
    }
close(PROTO);

foreach my $len (@list)
    {
        open( OUTFILE, ">$TMPDIR/$len.patt" ) or die "can't write\n";

        #print "writing $TMPDIR/$len.patt\n";
        foreach my $key ( keys %$len )
            {
                print OUTFILE ">$key\n$$len{$key}\n";
            }
        close(OUTFILE);
    }

# The user might want to revcomp the sequences - let them do so.

$wildseq = revcomp($wildseq) if $revcomp;
$mutseq  = revcomp($mutseq)  if $revcomp;

unless ( length($wildseq) == length($mutseq) )
    {
        die "Input sequences should be the same length.\nIf you have gaps, please indicate them with - characters.\n";
    }
else
    {
        if ( $wildseq =~ /-/ )
            {
                my $tmpseq = $wildseq;
                $wildseq = $mutseq;
                $mutseq  = $tmpseq;
            }

    }

if ( $mutseq =~ /-/ )
    {
        $gapseq = $mutseq;
        $mutseq =~ s/-//g;
    }

# Next we need to know the coordinates of the difference between the sequences...
# Iterate over the length of $wildseq, and compare it to $mutseq, or $gapseq, on
# a character-by-character basis. Store diffs in $ROI and @roi.

my @wildseq = split( //, $wildseq );
my @mutseq  = split( //, $mutseq );
my @gapseq  = split( //, $gapseq ) if $gapseq;

if ($gapseq)
    {
        for ( $cnt = 0; $cnt < $#wildseq; $cnt++ )
            {

                unless ( $wildseq[$cnt] eq $gapseq[$cnt] )
                    {
                        my $add = $cnt + 1;
                        push @roi, $add;
                    }
            }
    }
else
    {
        for ( my $cnt = 0; $cnt < $#wildseq; $cnt++ )
            {
                unless ( $wildseq[$cnt] eq $mutseq[$cnt] )
                    {
                        my $add = $cnt + 1;
                        push @roi, $add;
                    }
            }
    }

$ROI = join( ',', @roi );

#print "Region of interest is $ROI\n";

# Run fuzznuc against the two sequences (wildseq and mutseq *not* gapseq) - restsites might be created
# by a deletion.
# There is a problem when we start making oligo-based mutations to a mutant sequence which has a gap
# in it, however. Since the mutation process is sequence-specific, it makes sense to leave the gaps in
# during the in silico mutation process, since that is more likely to reflect the process in vitro.
# We want to scan for possible restriction sites on gapless sequence though, since we need to find
# unmutated res sites. Solution is to do this, then apply a subroutine later to detect restriction pattern
# differences between mutatted and unmutated sequences.
# Unfortunately, can't pipe the sequences to fuzznuc, need files.

open( SEQFILE, ">$TMPDIR/wildseq.fasta" );
print SEQFILE ">wildseq\n$wildseq\n";
close(SEQFILE);
open( SEQFILE, ">$TMPDIR/mutseq.fasta" );
print SEQFILE ">mutseq\n$mutseq\n";
close(SEQFILE);

foreach $len (@list)
    {

        # Set missmatch to one less than the pattern length.
        $miss     = $len - 1;
        $pattfile = "\@$TMPDIR/$len" . ".patt";
        @wildfuzz = `fuzznuc -pattern $pattfile -seq $TMPDIR/wildseq.fasta -pmismatch $miss -auto -stdout\n`;
        @mutfuzz  = `fuzznuc -pattern $pattfile -seq $TMPDIR/mutseq.fasta -pmismatch $miss -auto -stdout\n`;
        foreach my $line (@wildfuzz)
            {
                $line =~ s/^\s+//g;
                $line =~ s/\s+|:/\t/g;
                $line =~ s/\t$//g;
                unless ( $line =~ /#|Start/ || $line =~ /^$/ )
                    {
                        push @wildcuts, $line;
                    }
            }
        foreach my $line (@mutfuzz)
            {
                $line =~ s/^\s+//g;
                $line =~ s/\s+|:/\t/g;
                $line =~ s/\t$//g;
                unless ( $line =~ /#|Start/ || $line =~ /^$/ )
                    {
                        push @mutcuts, $line;
                    }
            }

    }

# The coordinates in @mutcuts are not going to match those in @wildcuts, in the case of an indel.
# So: produce a new @gapcuts with the coordinates adjusted.

if ($gapseq)
    {
        $gapstart = min @roi;
        $gaplen   = $#roi + 1;

        # We already tested (in the gap subroutine) whether the gap coordinates are contiguous.
        foreach $line (@mutcuts)
            {
                my @line = split( /\s+/, $line );

                if ( $line[0] > $gapstart )
                    {
                        $line[0] += $gaplen;
                    }
                if ( $line[1] > $gapstart )
                    {
                        $line[1] += $gaplen;
                    }
                push @gapcuts, "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]";
                print CHK "MUTT\t$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\n\n";
            }
    }

#close(CHK);

# Extract the deltas from the arrays. From this point on, we don't need @mutcuts if we've got @gapcuts.

if (@gapcuts)
    {
        @mutcuts = @gapcuts;
        $mutseq  = $gapseq;
    }

# Extract the lines from @wildcuts and @mutcuts which are not in the other - these are our potentially differentiating patterns.
undef %w_seen;
undef %m_seen;
undef @wild_only;
undef @mut_only;

foreach $line (@wildcuts)
    {
        $w_seen{$line} = 1;
    }
foreach $line (@mutcuts)
    {
        $m_seen{$line} = 1;
    }
foreach $line (@wildcuts)
    {
        unless ( $m_seen{$line} )
            {
                push @wild_only, $line;
            }
    }
foreach $line (@mutcuts)
    {
        unless ( $w_seen{$line} )
            {
                push @mut_only, $line;
            }
    }



# Combine the outputs, but only use interesting (right mut number) lines.
if ( $cut eq 'wild' )
    {
        push @alldiffs, @wild_only;
    }
elsif ( $cut eq 'mut' )
    {
        push @alldiffs, @mut_only;
    }
elsif ( $cut eq 'all' )
    {
        push @alldiffs, @mut_only, @wild_only;
    }
else
    {
        push @alldiffs, @mut_only, @wild_only;
    }

foreach $line (@alldiffs)
    {
        my @line = split( /\s+/, $line );
        if ( $line[5] eq $mutswant )
            {
                push @combined, $line;
            }

    }



# We'll also need the whole thing, to scan for interfering sites ($zero)
push @everycut, @mutcuts, @wildcuts;

foreach $line (@combined)
    {
        my @line = split( /\s+/, $line );
        if ( $mutswant == $line[5] )
            {
                if ( $scan eq 'all' && gap_ambig($line) eq 'True' )
                    {

                        push @repout, $line;
                        push @hitenz, $line[3] unless $seenhit{ $line[3] }++;
                    }
                elsif ( $scan eq 'left' && min(@roi) > max( getoffsets($line) ) && gap_ambig($line) eq 'True' )
                    {

                        push @repout, $line;
                        push @hitenz, $line[3] unless $seenhit{ $line[3] }++;

                    }
                elsif ( $scan eq 'right' && max(@roi) < min( getoffsets($line) ) && gap_ambig($line) eq 'True' )
                    {

                        push @repout, $line;
                        push @hitenz, $line[3] unless $seenhit{ $line[3] }++;
                    }
                elsif ( $scan eq 'straddle' && min( getoffsets($line) ) < min(@roi) && max( getoffsets($line) ) > max(@roi) && gap_ambig($line) eq 'True' )
                    {

                        push @repout, $line;
                        push @hitenz, $line[3] unless $seenhit{ $line[3] }++;

                    }

            }
    }



my %seenzero;
if ($zero)
    {
        foreach $line (@everycut)
            {
                my @line = split( /\s+/, $line );
                foreach $enzyme (@hitenz)
                    {

                        if ( $line[5] eq '.' && $line[3] eq $enzyme )
                            {
                                push @repout, $line unless $seenzero{$line}++;
                            }
                    }
            }

    }

foreach $enzyme ( sort(@hitenz) )
    {
        print "-" x length($wildseq);
        print "\n\n";
        print "$enzyme\t$protosite{$enzyme}\n";
        foreach $line (@repout)
            {
                my @line = split( /\s+/, $line );
                if ( $line[3] eq $enzyme )
                    {
                        print "\n";

                        print report($line);
                        print "\n";
                    }
            }
    }

# All of the temporary files go in $TMPDIR. Clean it out when we are finished.
clean();

sub getoffsets
    {
        # returns a list of mutation points, coordinates are on $wildseq or $mutseq
        my @mutpoints;
        my @line     = split( /\s+/, $_[0] );
        my @ambig    = split( //,    $line[4] );
        my @nonambig = split( //,    $line[6] );
        for ( my $cnt = 0; $cnt < length( $line[4] ); $cnt++ )
            {
                unless ( snuc( $ambig[$cnt], $nonambig[$cnt] ) )
                    {
                        my $addme = $cnt + $line[0];
                        push @mutpoints, $addme;
                    }
            }
        return @mutpoints;
    }

sub gap_ambig
    {

        #What if we have a condition like this?:
        #BetI    W^CCGGW
        #              v vvv
        #              WCCGGW
        #GGTCGATTGATCTTCCAATATGCTAGTTTCAACAACTCTCGTTCTTTACACTT
        #GGTCGATTGATCTTCCAAT-TGCTAGTTTCAACAACTCTCGTTCTTTACACTT
        #                   ^
        #12345678901234567890123456789012345678901234567890123
        #         1         2         3         4         5
        #fuzznuc says it will distinguish, but actually, it won't, because the deletion just moves a cuttable
        #nuc into the pattern. must check for this, only way is to use something like EMBOSS restrict
        my @mutlist;
        my $cutline    = $_[0];
        my $testwild   = $wildseq;
        my $testmut    = $mutseq;
        my @line       = split( /\s+/, $_[0] );
        my @res_patt   = split( //, $line[4] );
        my $res_patt_s = $line[4];
        my @seq_patt   = split( //, $line[6] );
        my $seq_patt_s = $line[6];

        for ( my $cnt = 0; $cnt < $#res_patt + 1; $cnt++ )
            {
                unless ( snuc( $seq_patt[$cnt], $res_patt[$cnt] ) )
                    {
                        push @mutlist, $cnt;
                    }
            }
        my $mutlist = join( ',', @mutlist );

        # Iterate over the mutations. In each case, use substr to actually make the mutation to the
        # wild type and mutant sequences, so we can compare the results of restrict(EMBOSS) to see if there is a difference.
        # It doesn't make much sense to try to mutate a gap character, and in practice the oligos wouldn't be designed to do that.
        # Do we remove the gap characters *before* or *after* making the substitution in the mutant sequence? We are only going to
        # call this subroutine if it has gaps, but since the mutation process is *sequence* specific rather than *site* (i.e. offset)
        # specific, I'm making the substitution first.
        foreach $pos (@mutlist)
            {

                # get offset to mutating position in actual sequence, $pos is the posn in the pattern
                $adjpos = $pos + $line[0] - 1;

                # find the nucleotide at that position
                my $neednuc = substr( $res_patt_s, $pos, 1 );

                # mutate wild-type sequence (won't have any gaps)
                substr( $testwild, $adjpos, 1 ) = $neednuc;

                # check if the posn in the mutant sequence is a gap, if not, mutate
                my $change_mut_nuc = substr( $testmut, $adjpos, 1 );
                unless ( $change_mut_nuc eq '-' )
                    {
                        substr( $testmut, $adjpos, 1 ) = $neednuc;
                    }
            }

        # lose the gaps in the mutant sequence
        $testmut =~ s/-//g;

        open( WILD, ">$TMPDIR/wild" );
        open( MUT,  ">$TMPDIR/mut" );
        print WILD ">wild\n$testwild\n";
        print MUT ">mut\n$testmut\n";
        close(WILD);
        close(MUT);
        my $wild_cuts_rest_line = `restrict $TMPDIR/wild -warning 0 -error 0 -enzymes $line[3] -stdout -auto | grep -c +`;
        my $mut_cuts_rest_line  = `restrict $TMPDIR/mut -warning 0 -error 0 -enzymes $line[3] -stdout -auto | grep -c +`;
        my $restrict_out        = `restrict $TMPDIR/wild -warning 0 -error 0 -enzymes $line[3] -stdout -auto`;

        #print GDBG $restrict_out;
        my $restrict_out = `restrict $TMPDIR/mut -warning 0 -error 0 -enzymes $line[3] -stdout -auto`;

        #print GDBG $restrict_out;
        chomp $wild_cuts_rest_line;
        chomp $mut_cuts_rest_line;

        if ( "$wild_cuts_rest_line" eq "$mut_cuts_rest_line" )
            {
                return "False";
            }
        else
            {
                return "True";
            }

    }

sub report
    {
        my $report;
        my @input    = split( /\t+/, $_[0] );
        my @ambig    = split( //,    $input[4] );
        my @nonambig = split( //,    $input[6] );
        for ( 2 .. $input[0] ) { $report .= ' '; }
        for ( my $cnt = 0; $cnt < length( $input[4] ); $cnt++ )
            {
                if ( snuc( $ambig[$cnt], $nonambig[$cnt] ) )
                    {
                        $report .= ' ';
                    }
                else
                    {
                        $report .= 'v';
                    }
            }
        $report .= "\n";
        for ( 2 .. $input[0] ) { $report .= ' '; }
        $report .= "$input[4]\n";
        $report .= "$wildseq\n";
        $report .= "$mutseq\n";
        $end = min(@roi) - 3;
        for ( 0 .. ( $end + 1 ) ) { $report .= " "; }
        for ( 0 .. $#roi ) { $report .= "^"; }
        $report .= "\n";
        $report .= scale($wildseq);
        return $report;

    }

sub debug
    {
        my $dreport;
        my @input    = split( /\t+/, $_[0] );
        my @ambig    = split( //,    $input[4] );
        my @nonambig = split( //,    $input[6] );
        $dreport .= "Debug\n";
        $dreport .= "$_[0]\n";

        my @offset = getoffsets( $_[0] );
        my $max    = max(@offset);
        my $min    = min(@offset);
        $dreport .= "Offsets = max $max min $min\n";
        $dreport .= $offset;
        $dreport .= "\nROI =\t" . $ROI . "\n";

        return $dreport;
    }

sub scale
    {
        #This just creates a nifty scale, the same length as $_[0]
        my $ret;
        $length = length( $_[0] );
        $tens   = 0;

        #Units
        for ( $cnt = 1; $cnt < $length + 1; $cnt++ )
            {
                $ret .= ( $cnt % 10 );

            }
        $ret .= "\n";

        #Tens
        for ( $cnt = 1; $cnt < $length + 1; $cnt++ )
            {
                if ( ( $cnt % 10 ) == 0 )
                    {
                        $tens++;
                        $ret .= "$tens";
                        if ( $tens == 9 ) { $tens = 0 }
                    }
                else
                    {
                        $ret .= " ";
                    }
            }
        $ret .= "\n\n";
        return $ret;
    }

sub snuc
    {
        # Input is two nucleotide codes, one possibly an ambiguity code. Returns true if they match
        my @A = ( 'A', 'R', 'W', 'M', 'D', 'H', 'V', 'N' );
        my @G = ( 'G', 'R', 'S', 'K', 'B', 'D', 'V', 'N' );
        my @C = ( 'C', 'Y', 'S', 'M', 'B', 'H', 'V', 'N' );
        my @T = ( 'T', 'Y', 'W', 'K', 'B', 'D', 'H', 'N' );
        my $pri;
        my $sec;
        if ( uc( $_[0] ) eq 'A' | uc( $_[0] ) eq 'C' | uc( $_[0] ) eq 'G' | uc( $_[0] ) eq 'T' )
            {

                $pri = uc( $_[0] );
                $sec = uc( $_[1] );
            }
        else
            {
                $pri = uc( $_[1] );
                $sec = uc( $_[0] );
            }

        if ( $pri eq 'A' )
            {

                foreach $e (@A)
                    {
                        if ( $sec eq $e )
                            {
                                return '1';
                            }
                    }
            }
        elsif ( $pri eq 'C' )
            {
                foreach $e (@C)
                    {
                        if ( $sec eq $e )
                            {
                                return '1';
                            }
                    }
            }
        elsif ( $pri eq 'G' )
            {
                foreach $e (@G)
                    {
                        if ( $sec eq $e )
                            {
                                return '1';
                            }
                    }
            }
        elsif ( $pri eq 'T' )
            {
                foreach $e (@T)
                    {
                        if ( $sec eq $e )
                            {
                                return '1';
                            }
                    }
            }

    }

sub revcomp
    {

        my $revcomp = reverse( $_[0] );
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
    }

sub uncut
    {

        # This just removes the cruft from the pattern field in PROTO
        my $site = $_[0];
        $site =~ s/\^//g;
        $site =~ s/,//g;
        return $site;
    }

sub clean
    {

        my $out = `rm -f /tmp/*patt /tmp/*fasta /tmp/longout /tmp/shortout > /dev/null 2>&1`;
        my $out = `rm -rf $TMPDIR`;

    }


