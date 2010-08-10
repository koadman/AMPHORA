use strict;
use Cwd;
use Getopt::Long;

my $usage = qq~
Usage: perl -w $0 <options>
        Options:
                -AMPHORA_home: the directory where AMPHORA will be installed
                -Bioperl_home: the diretory where bioperl 1.5.2 or later is installed
~;

my ($AMPHORA_home, $bioperl_home) = ();
GetOptions(	'AMPHORA_home=s'=>\$AMPHORA_home,
			'Bioperl_home=s'=>\$bioperl_home) || die $usage;

die $usage unless ($AMPHORA_home and $bioperl_home);

die "$AMPHORA_home does not exist\n" unless (-e $AMPHORA_home);

my $working_dir = getcwd(); 
my $hmmer_path = get_program_path('hmmpfam');
my $blastp_path = get_program_path('blastp');
generate_perl_scripts();
compile_seqboot();
compile_raxml();
compile_quicktree();
deploy();

########################################################################################
sub get_program_path {
	my $program = shift;
    my $path=`which $program`;
    if($path=~/$program$/){
        $path=~s/\/*$program[^\/]*$//;
    }
    else {
    	die "Can't find $program, is it installed?\n";
    }
	return $path;
}

sub generate_perl_scripts {	
	open (IN, "find $working_dir/Scripts -type f |") || die;
	while (<IN>) {
		chop;
		my $script = $_;
		my $text = undef;
		open (FILE, $script) || die "Can't open $script\n";
		while (<FILE>) {
			$text .= $_;
		}
		close FILE;
		$text =~ s/AMPHORA_home_dir/$AMPHORA_home/g;
		$text =~ s/bioperl_home_dir/$bioperl_home/g;
		$text =~ s/HMMER_path/$hmmer_path/g;
		$text =~ s/BLASTP_path/$blastp_path/g;
		
		open (FILE, ">$script") || die "Can't write $script\n";
		print FILE $text;
		close FILE;
	}
}

sub compile_seqboot {
	print STDERR "\nCompiling seqboot ...\n";
	chdir "$working_dir/src/Seqboot";
	system("make clean");
	system("make seqboot");
	system("mv $working_dir/src/Seqboot/seqboot $working_dir/bin/.");
}

sub compile_raxml {
	print STDERR "\nCompiling RAxML ...\n";
	chdir "$working_dir/src/RAxML-7.0.4";
	system("make -f Makefile.gcc");
	system("mv $working_dir/src/RAxML-7.0.4/raxmlHPC $working_dir/bin/.");
}

sub compile_quicktree {
	print STDERR "\nCompiling quicktree ...\n";
	chdir "$working_dir/src/quicktree_1.1";
#	system("make clean");
	system("make quicktree");
	system("mv bin/quicktree $working_dir/bin/.");
}

sub deploy {
	print STDERR "\nCopying files to $AMPHORA_home ...\n";
	chdir $working_dir;
	for ('Scripts','bin','Marker','Reference','Taxonomy') {
		system("cp -r $_ $AMPHORA_home/.");
		system("rm -r $_");
	}
	system("chmod +x $AMPHORA_home/Scripts/*");
	system("chmod +x $AMPHORA_home/bin/*");
	print STDERR "AMPHORA successfully installed\n";
}	