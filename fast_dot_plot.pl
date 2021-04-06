#!/usr/bin/perl -w
# NOTE: YOU MUST CHANGE THE LINE ABOVE TO POINT TO
# THE FULL PATH OF THE PERL EXECUTABLE ON YOUR SYSTEM.

# System requirements: requires perl and perl modules
# DB_File and GD.pm, available from CPAN.
 
# For usage information, run this program with the flag
# -h.

#To produce plot with printed word size, step size, and sequence length(s), uncomment lines 378, 384 and 387

#Code provided by H. Skaletsky

require v5.6.0;

use GD;
use DB_File;
use strict;
use vars qw/$VERSION/;

BEGIN {

	    $VERSION = '1.0';

	    if (!defined $GD::VERSION || $GD::VERSION ne 1.32) {
			print STDERR "\n$0:\n";
			print STDERR "WARNING -- GD VERSION 1.32 REQUIRED.\n";		
			if (defined $GD::VERSION){
				print STDERR "You have GD version $GD::VERSION.\n";
				print STDERR "If $0\n";
				print STDERR "does not work correctly you will\n";
				print STDERR "have to install GD version 1.32.\n\n";
		    }
		}

	    if (!defined $DB_File::VERSION || $DB_File::VERSION ne 1.814) {
			print STDERR "\n$0:\n";
			print STDERR "WARNING -- DB_File VERSION 1.814 REQUIRED.\n";
			if (defined $DB_File::VERSION){
				print STDERR "You have DB_File version $DB_File::VERSION.\n";
				print STDERR "If $0\n";
				print STDERR "does not work correctly you will\n";
				print STDERR "have to install DB_File version 1.814.\n\n";
		    }
		}
}

use Getopt::Long;

sub usage();
sub main();
sub print_square($$$$$$$$$$);
sub print_triangle($$$$$$$$$);
sub print_png ($$);
sub revcomp($);


main();

sub usage() {
    print qq/
    
USAGE: $0 -i <in file> [optional args]

OPTIONAL ARGUMENTS:

        -j <2nd file> 
        -o <out file>
        -d <dot file>
        -w <word size> 
        -s <step size>
        -p <pixels>
        -t <title>
          
        -h print this help message


Create a PNG ("portable network graphics") file
that displays a triangular dot plot of the input
sequence against itself, or a square plot against 
a second sequence with option -j.





<in file>   is a fasta format file from which the
            sequence is taken.

            This is the only required parameter.
			


<2nd file>  is a fasta format file from which a
            second sequence is taken.
	
            default:
            null -- triangular dot plot created           
			


<out file>  is the PNG file created.

            default:
            <in file(s)>_<word_size>_<step_length>_<pixels>.png
			


<dot file>  A database file containing the 
            positions of words of <word size>
            encountered in your sequence.
            This database is not human-readable,
            and can get quite large, so make 
            sure that you have plenty of disk 
            space available. 
            
            Databases size scales linearly 
            with sequence length. With default 
            parameters,databases are are roughly 
            150k per base-pair of sequence.
            
            Complex sequences use more disk
            space than repetitive ones.
			
            default:
            <infile(s)>_<word_size>_<step_length>.dot
			


<word size> is the word size for a match.  A dot
            is printed if there is a perfect match
            of length <word size>.
            
            default:
            <word size> = 100
			


<step size> is the number of bases to move the
            word for each dot.
            
            default:
            <step size> = 1
			


<pixels>    is the width of the plot in pixels and
            defines the resolution of the image.

            default:
            <pixels> = 600
			


<title>     is a title to place in the output.

            default:
            null -- no title appears



If for some reason this program is disrupted during
a run, restart with the same parameters and it will
automatically resume where it left off.

Completed runs can be re-plotted from the existing 
dot file using a different resolution by changing -p

-h causes this message to printed.

(Version $VERSION)

/;
}

sub main() {

	#Define the arguments
    my ($seqfile, $second, $outfile, $dotfile, $word, $step, $pixels, $title, $help);
	
	#Perform Checking -- remind the user that TMTOWTDI.
    if (!GetOptions(
    	'infile=s'  => \$seqfile,
    	'jfile=s'   => \$second,
		'wordlen=i' => \$word,
		'step=i'    => \$step,
		'outfile=s' => \$outfile,
		'title=s'   => \$title,
		'pixels=i'   => \$pixels,
		'dotfile=s' => \$dotfile,
		'help'      => \$help)) {
		usage;
		exit -1;
    }
    
    if ($help) {
		usage;
		exit;
    }
    
	if (!defined $seqfile){
		usage; exit -1;
    }
    
	if (!defined $second){
		print "Creating triangular plot from: $seqfile\n";
		print "\t use -j to make a square plot\n";
	}
	
    if (!defined $word){
    	$word = 100;
    	print "Using word size 100\n";
		print "\tuse -w to set word size\n";
    }elsif ($word <= 0){
		print STDERR "$0 -w <word size> must be >= 0\n";
		exit -1;
    }
    
	if (!defined $step){
		$step = 1;
    	print "Using step size 1\n";
		print "\tuse -s to set step size\n";
	}elsif ($step <= 0) {
		print STDERR "$0 -s <step length> must be >= 0\n";
		exit -1;
    }	
    
    if (!defined $pixels){
    	$pixels = 600;
    	print "Using image size 600 pixels\n";
		print "\tuse -p to set pixels\n";
	}elsif ($pixels <= 0) {
		print STDERR "$0 -s <pixels> must be >= 0\n";
		exit -1;
    }	
	
	if (!defined $dotfile) {
		$seqfile =~ m/\/?([^\/\.]+)\.?\w*$/;
		$dotfile = $1;
		if ($second){
			$second =~ m/\/?([^\/\.]+)\.?\w*$/;
			$dotfile .= "_" . $1;
		}
		$dotfile .= "_" . $word . "_" . $step;
		#$dotfile =~ tr/[\.\/]/_/;
		$dotfile .= ".dot";
		print "Storing matches in: $dotfile\n";
		print "\tuse -d to name dot file\n";
    }
		
	if (!defined $outfile){
		$seqfile =~ m/\/?([^\/\.]+)\.?\w*$/;
		$outfile = $1;
		if ($second){
			$second =~ m/\/?([^\/\.]+)\.?\w*$/;
			$outfile .= "_" . $1;
		}
		$outfile .= "_" . $word . "_" . $step . "_" . $pixels;
		#$outfile =~ tr/[\.\/]/_/;
		$outfile .= ".png";
		print "Printing picture in: $outfile\n";
		print "\tuse -o to name out file\n";
    }

	    
    $title = '' unless defined $title;

	#Open the database file -- this makes us fast!
	my %DOT = ();
	dbmopen(%DOT, "$dotfile", 0666) || die "Cannot open $dotfile: $!";

	#Define some variables we need
	my ($xlen, $ylen) = (0, 0);
	my ($i, $j, $k, $key) = (0, 0, 0, "");
	my $seq = "";
	
	#Start reading the first sequence, avoid a file slurp
	
    open(IN, $seqfile) || die "Cannot open $seqfile: $!";
    while (<IN>){
    	if (m/^>/){
    		next;
    	}
    	chomp;
    	$k += length $_;
    	$seq .= uc($_);
    	while ($seq && (($i + $word) < $k) && ($j + $word) < $k){
    		#wait until we hit the place where we left off!
    		if ($DOT{xpos} && $i < $DOT{xpos}){
    		}elsif ($i < $j){
    			$DOT{xpos} = $i;
    		}elsif ($i == $j){
    			$DOT{xpos} = $i;
    			$key = substr($seq, 0, $word);
    			#treat Ns as mismatches
    			unless ($key =~ m/N/){
    				if ($DOT{$key}){
    					$DOT{$key} .= "x:$i\t";
    				}elsif($DOT{revcomp($key)}){
    					$DOT{revcomp($key)} .= "x:$i\t";
    				}else{
    					$DOT{$key} .= "x:$i\t";
    				}
    			}
    			$j += $step;
    		}
    		substr($seq, 0, 1) = "";
    		$i++;
    	}
	    	
	}
	close IN;
	$xlen = $k;
	
	#Is there another sequence?

	    
	if ($second){

		#Looks like there is.
		#Flush out these variables:

		($i, $j, $k, $key) = (0, 0, 0, "");
		$seq = "";
		
		#Read the second sequence, avoid a file slurp
		
    	open(IN, $second) || die "Cannot open $second: $!";	
	    while (<IN>){
	    	if (m/^>/){
	    		next;
	    	}
	    	chomp;
	    	$k += length $_;
	    	$seq .= uc($_);
	    	while ($seq && (($i + $word) < $k) && ($j + $word) < $k){
				#wait until we hit the place where we left off!
	    		if ($DOT{ypos} && $i < $DOT{ypos}){
	    		}elsif ($i < $j){
	    			$DOT{ypos} = $i;
	    		}elsif ($i == $j){
	    			$DOT{ypos} = $i;
	    			$key = substr($seq, 0, $word);
	    			#treat Ns as mismatches
	    			unless ($key =~ m/N/){
	    				if ($DOT{$key}){
	    					$DOT{$key} .= "y:$i\t";
	    				}elsif($DOT{revcomp($key)}){
	    					$DOT{revcomp($key)} .= "y:$i\t";
	    				}else{
	    					$DOT{$key} .= "y:$i\t";
	    				}
	    			}
	    			$j += $step;
	    		}
	    		substr($seq, 0, 1) = "";
	    		$i++;
	    	}
	    	
	    }
	    close IN;
	    $ylen = $k;
	}
	# Close that DB
	
	dbmclose %DOT;
	
    # Create and print the output.
    
    my ($img, $white, $black, $gray, $margin, $x0, $y0, $x1, $y1);
	$margin = 50;
	$x0 = $y0 = $margin;
	$x1 = $y1 = $pixels - $margin;
    $img = new GD::Image($pixels, $pixels);
    $white = $img->colorAllocate(255,255,255);
	$black = $img->colorAllocate(0,0,0);
	$gray = $img->colorAllocate(187,187,187);
    $img->interlaced('true');
    #$img->string(gdLargeFont, $x0, $y0 - 35, "$title  -w=$word -s=$step", $black);
	
	
    # Square or Triangle?
    
    if ($second){
    	#$img->string(gdLargeFont, $x0, $y0 - 15, "$xlen bp x $ylen bp", $black);
    	print_square($img, $x0, $x1, $y0, $y1, $black, $gray, $xlen, $ylen, $dotfile);
    }else{
		#$img->string(gdLargeFont, $x0, $y0 - 15, "$xlen bp", $black);
		print_triangle($img, $x0, $x1, $y0, $y1, $black, $gray, $xlen, $dotfile);
    }
    
    # Print this thing out
    
    print_png($img, $outfile);
    
    print "Plot Complete!\n";
} 
 

sub print_square ($$$$$$$$$$){

	my ($img, $x0, $x1, $y0, $y1, $black, $gray, $xlen, $ylen, $dotfile) = @_;
	
	my $del = ($x1 - $x0) / $xlen;
	if ($del > (($y1 - $y0) / $ylen)){
		$del = (($y1 - $y0) / $ylen);
	}	
	my %hash = ();
	my ($xd, $yd, $key, $value);
	my %TEMPFILE = ();
	dbmopen(%TEMPFILE, "$dotfile", 0666) || die "Can't open database $dotfile: $!";
	#avoid slurping again, pull values out one by one!
	while (($key,$value) = each %TEMPFILE){
		if ($key =~ m/[ACGT]+/){
			while ($value =~ m/([xy]):(\d+)\t(.*)/){
				($hash{$1}{$2}, $value) = (1, $3);	
			}
		}
		foreach $xd (keys(%{$hash{x}})){
			foreach $yd (keys(%{$hash{y}})){
				$img->setPixel(($x0 + ($del * $xd)), ($y0 + ($del * $yd)), $black);
			}
		}
		%hash = ();
	}
	dbmclose %TEMPFILE;
	$img->rectangle($x0, $y0, ($x0 + ($del * $xlen)), ($y0 + ($del * $ylen)), $gray);
}

 
sub print_triangle ($$$$$$$$$){

	my ($img, $x0, $x1, $y0, $y1, $black, $gray, $xlen, $dotfile) = @_;
	
	my $del = ($x1 - $x0) / $xlen;
	my %hash = ();
	my ($xd, $yd, $key, $value);
	my %TEMPFILE = ();
	dbmopen(%TEMPFILE, "$dotfile", 0666) || die "Can't open database $dotfile: $!";
	#avoid slurping again, pull values out one by one!
	while (($key,$value) = each %TEMPFILE){
		if ($key =~ m/[ACGT]+/){
			while ($value =~ m/x:(\d+)\t(.*)/){
				($hash{x}{$1}, $hash{y}{$1}, $value) = (1, 1, $2);	
			}
		}
		foreach $xd (keys(%{$hash{x}})){
			foreach $yd (keys(%{$hash{y}})){
					$img->setPixel(($x0 + (($del * ($xd + $yd)) * 0.5)), ($y1 - (($del * sqrt(($yd - $xd) * ($yd - $xd))) * 0.5)), $black);
			}
		}
		%hash = ();
	}
	dbmclose %TEMPFILE;
	$img->line($x0, $y1, $x1, $y1, $gray);
	$img->line($x0, $y1, (($x0 + $x1) * 0.5), ($y1 - (($x1 - $x0) * 0.5)), $gray);
	$img->line((($x0 + $x1) * 0.5), ($y1 - (($x1 - $x0) * 0.5)), $x1, $y1, $gray);
		
}

sub print_png ($$) {
    my ($img, $outfile) = @_;
    if ($outfile) {
	open(OUT, ">$outfile")
	    || die "Cannot write $outfile: $!\n";
	print OUT $img->png;
	close OUT;
    } else {
	print $img->png;
    }
}

sub revcomp ($){
	my $f = shift(@_);
	my $r = reverse($f);
	$r =~ tr/GATC/CTAG/;
	return $r;
}