#!/usr/bin/perl
#
# Insert one H3 antibody PDB structure into another at specified residue indices. 
#
# by Mon (17/11/2014)
#
# Modify by Ivan (12/04/2019)

use strict;
use feature 'say';
use Storable;
# use File::Basename;

if( !($#ARGV == 2 || $#ARGV == 3) )
{ 
	print "USAGE:\n\t$0 <PDB_A> <PDB_B> <Output_PDB> [noanchors]\n\n";
	print "DESCRIPTION:\n\tInsert one PDB structure (A) into another (B) at residue indices taken from N- and C-terminal ends of PDB A. Optionally, type \"noanchors\" to
 retain B structure anchors, otherwise they will be substituted by A structure ones.\n\n";
	exit 
}

# Loading things        
my $debug=0;
my $prog = "insertpdb";
my $fpdbA = $ARGV[0];
my $fpdbB = $ARGV[1];
my $fpdb_out = $ARGV[2];
my $noanchors = 0;
if($#ARGV == 3 && $ARGV[3] eq 'noanchors')
{
	print "$prog> Retaining PDB_B's anchors... (noanchors option enabled)\n";
	$noanchors = 1;
}
if($#ARGV == 3 && $ARGV[3] ne 'noanchors')
{
	print "$prog> sorry Wrong option... forcing exit\n";
	exit
}

#my $start = $ARGV[3];
#my $end = $ARGV[4];

# Loading PDB of structure A
my $completeA=readPDB($fpdbA); # Read the complete PDB (it's an array of hashes)
print "$prog> Structure A loaded from $fpdbA\n";

# Loading PDB of structure B
my $pdbB=readPDB($fpdbB); # Read the complete PDB (it's an array of hashes)
print "$prog> Structure B loaded from $fpdbB\n";

#Parte de ivan para captar loops H3 de dentro de la estructura completa en formato chothia.
#$pdbA = $completeA;
my $pdbA;
my $cur_res=0;
my $activateloop=0;
my $i=0;
foreach my $line (@{$completeA}) {
	my $res_num=$line->{resseq}; # get residue number
	my $chain =$line->{chain};
	if ($chain eq "H"){
		if (!($res_num eq $cur_res)){
			$cur_res =$res_num;
			if ($res_num eq 93){$activateloop=1;}
			if ($res_num eq 103){$activateloop=0;}
		}
	}
	
	if ($activateloop){
		$pdbA ->[$i]=$line;
		$i++;
	}
}


# Get chain id first atom of PDB A
my $chain = $pdbA->[0]->{chain};

# Get start and end indices from PDB A
my $lastatom = $#{$pdbA};
my $start = $pdbA->[0]->{resseq};
my $end = $pdbA->[$lastatom]->{resseq};
print "$prog> Loops in PDB A start at $start and end at $end of chain $chain.\n";

# writePDB($pdb,"full.pdb");
# writePDB($mpdb->[0],"model.pdb");

# Inserting cluster heads into the complete structure...
print "$prog> Inserting $fpdbA into the structure $fpdbB ... ";

# Duplicate complete structure
#my $pdb2 = Storable::dclone($pdbB);

# Trim loop anchors... if requested
if($noanchors == 1)
{
	print "$prog> Removing anchors from $fpdbA.\n";
	getchunkPDB($pdbA,$chain,$start+1,$end-1);
	#writePDB($pdbA,"removed.pdb");

	# Inserting loop structure ($mpdb) into the complete structure ($pdb)
	insertPDB($pdbA,$pdbB,$chain,$start+1,$end-1);
}
else
{
	# Inserting loop structure ($mpdb) into the complete structure ($pdb)
	insertPDB($pdbA,$pdbB,$chain,$start,$end);
}
writePDB($pdbB,$fpdb_out);
print "done! Output structure: $fpdb_out\n";
print "$prog> That's all folks!\n";
      
#################################################################
# FUNCTIONS
#################################################################


# Delete residues from PDB, only those between $start and $end will be retained.
sub getchunkPDB
{
	my $pdb = shift; # complete structure
	my $chain = shift; # chain id of the complete structure
	my $start = shift; # first residue index of the loop
	my $end = shift; # last residue index of the loop

	# Screen complete struture to get insertion/deletion indices...
	my $i=0;
	my $first=0;
	my $last=0;
	foreach my $full (@{$pdb}) # screen atoms
	{
		$i++;
		$first = $i if($full->{resseq} == $start-1 && $full->{chain} eq $chain);
		$last = $i if($full->{resseq} == $end && $full->{chain} eq $chain);
	}
#	splice @{$pdb},$first,$last-$first; # delete elements from array.
	splice @{$pdb},$last,$#{$pdb}; # delete elements from array.
	splice @{$pdb},0,$first; # delete elements from array.
}

# Inserting loop ($loop) into a complete structure ($pdb)
# Warning: - Both structures should have the same residue-numeration.
#          - Chain id of the complete structure must be provided.
sub insertPDB
{
	my $loop = shift; # loop structure
	my $pdb = shift; # complete structure
	my $chain = shift; # chain id of the complete structure
	my $start = shift; # first residue index of the loop
	my $end = shift; # last residue index of the loop

	# Screen complete struture to get insertion/deletion indices...
	my $i=0;
	my ($first,$last)=(0,0);
	foreach my $full (@{$pdb}) # screen atoms
	{
		$i++;
		$first = $i if($full->{resseq} == $start-1 && $full->{chain} eq $chain);
		# $last = $i if($full->{resseq} == $end && $full->{chain} eq $chain);
	}

	for($i = $#{$pdb} ; $i >= 0 ; $i-- ) 
	{ 
		my $full = ${$pdb}[$i];
		$last = $i if($full->{resseq} == $end+1 && $full->{chain} eq $chain);
	}

	print "$prog> start= $start  end= $end  -->  first= $first  last= $last\n";
	splice @{$pdb},$first,$last-$first,@{$loop}; # delete and insert one array into other array.
}

sub readPDBLine 
{
        my $line = shift;
        my $newAt = {};
        $newAt->{tag}=substr($line, 0,6);
        $newAt->{seq}=int(substr($line, 7,4));
        $newAt->{name}=substr($line, 12,4);
        $newAt->{altloc}=substr($line, 16,1);
        $newAt->{residueId}=substr($line,17,3);
        $newAt->{chain}=substr($line,21,1);
        $newAt->{resseq}=int(substr($_,22,4));
				$newAt->{iCode}=substr($line,26,1);
        $newAt->{x}=substr($line,30,8);
        $newAt->{y}=substr($line,38,8);
        $newAt->{z}=substr($line,46,8);
        $newAt->{occ}=substr($line, 54,6);
        $newAt->{Bfact}=substr($line, 60,6);
        $newAt->{charge}=$newAt->{occ};
        $newAt->{type}=substr($newAt->{name},1,1);
        return $newAt;
# Atomic Coordinate Entry Format Version 3.2
# http://www.wwpdb.org/documentation/format32/sect9.html
# COLUMNS        DATA  TYPE    FIELD        DEFINITION
#-------------------------------------------------------------------------------------
# 1 -  6        Record name   "ATOM  "
# 7 - 11        Integer       serial       Atom  serial number.
#13 - 16        Atom          name         Atom name.
#17             Character     altLoc       Alternate location indicator.
#18 - 20        Residue name  resName      Residue name.
#22             Character     chainID      Chain identifier.
#23 - 26        Integer       resSeq       Residue sequence number.
#27             AChar         iCode        Code for insertion of residues.
#31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
#39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
#47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
#55 - 60        Real(6.2)     occupancy    Occupancy.
#61 - 66        Real(6.2)     tempFactor   Temperature  factor.
#77 - 78        LString(2)    element      Element symbol, right-justified.
#79 - 80        LString(2)    charge       Charge  on the atom.
}

# Writes a pdb-file form a "array of hashes" format
sub writePDB 
{
  # "%6c%5c%*c%4c%c%3c%*c%c%4d%c%*3c%8c%8c%8c%6c%6c%c%3c"
  my $data=shift;
  my $file=shift;
  my $renum=shift;
  open(PDB,">$file") or die "\nFailed to open $file\n";
  if($renum eq undef) # no renumber atom index
  {
	  foreach(@{$data})
	  {
	    printf PDB "%6s%5s %4s%s%3s %s%4s%s   %8s%8s%8s%6s%6s\n",
	      $_->{tag},
	      $_->{seq},
	      $_->{name},
	      $_->{altloc},
	      $_->{residueId},
	      $_->{chain},
	      $_->{resseq},
	      $_->{iCode},
	      $_->{x},
	      $_->{y},
	      $_->{z},
	      $_->{occ},
	      $_->{Bfact},
	  }
  }
  else # renumber atom index
  {
	  my $index=1;
	  foreach(@{$data})
	  {
	    printf PDB "%6s%5s %4s%s%3s %s%4s%s   %8s%8s%8s%6s%6s\n",
	      $_->{tag},
	      $index,
	      $_->{name},
	      $_->{altloc},
	      $_->{residueId},
	      $_->{chain},
	      $_->{resseq},
				$_->{iCode},
	      $_->{x},
	      $_->{y},
	      $_->{z},
	      $_->{occ},
	      $_->{Bfact};
	    $index++;
	  }
  }
  close(PDB);
}

# Reads a pdb-file from a file
sub readPDB
{
  my $file=shift; # PDB file
  my @data=(); # The readed PDB will be stored here.
  
  my $debug=0;
  my $cont=0;
  open(PDB,"$file") or die "\nFailed to open $file\n";
  while(<PDB>)
  {
    next unless /^ATOM/; # reading only ATOM begining lines
    push(@data,readPDBLine($_));
#    printf "%5d bfact= %10.2f\n",$cont+1,$data[$cont]->{Bfact} if $debug;
    $cont++;
  }
  close PDB;
  return(\@data);
}

