#!/usr/bin/perl

#Script for the multiple analysis of the outputs of refine_loops.pl

#Iván Martín Hernández


use strict;
use warnings;
use diagnostics;
use Data::Dumper qw(Dumper);
use feature 'say';

# INPUT PARSER
if( !($#ARGV == -1) )
{
	say "USAGE:\n\t$0 <identifier> <input1> <input2>\n";
	say "\t<identifier>  --> Input some partial name that identifies only the files that you want to analyze.";
	say "\t<input1>  --> Number of decoys that RCD made.";
	say "\t<input2>  --> Number of decoys that RCD gave as output per loop (the best under Korp rank).";
	say "\nDESCRIPTION:\n\tSay something here...";
	say "\nWARNING:\n\tSomething....";
	say "\nEXAMPLE:";
	say "\tAnalyze.pl native_rcd10b_1n1Hk5k99t 100000 5000";
	say "\nHELP: some help\n";
	exit;
}

my %short = (
   "ALA" => "A",
	"CYS" => "C",
   "CYX" => "C",
   "ASP" => "D",
   "GLU" => "E",
   "PHE" => "F",
   "GLY" => "G",
   "HIS" => "H",
   "HSD" => "H",
   "HSE" => "H",
   "HSP" => "H",
   "HIP" => "H", # Rosseta
   "HID" => "H", # Rosseta
   "HIE" => "H", # Rosseta
   "ILE" => "I",
   "LYS" => "K",
   "LEU" => "L",
   "MET" => "M",
   "ASN" => "N",
   "PRO" => "P",
   "GLN" => "Q",
   "ARG" => "R",
   "SER" => "S",
   "THR" => "T",
   "VAL" => "V",
   "TRP" => "W",
   "TYR" => "Y",
   "  A" => "A",
   "  G" => "G",
   "  T" => "T",
   "  C" => "C",
   "  U" => "U",
   " DA" => "A",
   " DG" => "G",
   " DT" => "T",
   " DC" => "C" ) ;

# Loading things
my $entire_name = "*.pdb";
my @files = glob("$entire_name");
open(OUTPUT,'>',"AB_loop_RCD.txt");

foreach my $PDB(@files){
   say $PDB;
	open (PDB, $PDB) or die "Can't open file $PDB.\n";
	my $newpdb=substr($PDB,0,4)."_2.pdb";
	open(OUTPUT2,'>',$newpdb);
	
	#Set loop variable
	my $L1_start;
	my $L2_start;
	my $L3_start;
	my $H1_start;
	my $H2_start;
	my $H3_start;
	my $L1_end;
	my $L2_end;
	my $L3_end;
	my $H1_end;
	my $H2_end;
	my $H3_end;
	my $L1_seq="";
	my $L2_seq="";
	my $L3_seq="";
	my $H1_seq="";
	my $H2_seq="";
	my $H3_seq="";
	my $L1_control=0;
	my $L2_control=0;
	my $L3_control=0;
	my $H1_control=0;
	my $H2_control=0;
	my $H3_control=0;
	
		
	my $cur_res = -1;
	my $new_res_num=0;
	my $cur_chain ="";
	my $pre_alt =" ";
	foreach my $line (<PDB>) {
		chomp($line);
		#print "$line\n";	
		
		if ($line =~ /^ATOM/ || $line =~ /^HETATM/) {
			my $format = readPDBLine($line);
			#my @fields = split( / +/, $line); 
			#my $res_num = $fields[5];
			my $res_num=substr($line,22,5); # get residue number
			$res_num=~s/\s+//g;
			#my $chain =$fields[4];
			#my $chain =substr($line,21,1);
			#my $res =substr($line,17,3);
			#my $res =$fields[3];
			my $chain =$format->{chain};
			my $res =$format->{residueId};
			#print "$res_num - $chain - $res\n";
			
			my $alt = $format->{altloc};
			
			if (!($cur_chain eq $chain)){
				#say $line;
				$cur_res = -1;
				$new_res_num=0;
				$cur_chain = $chain;
				}

			
			if (!($res_num eq $cur_res)){
				#print "$cur_res - $res_num - $new_res_num - $short{$res}\n";
				$cur_res =$res_num;
				$new_res_num =$new_res_num+1;
				#say $res;
				$pre_alt =" ";
				
				if ($chain eq "L"){
					#print "$cur_res - $res_num - $new_res_num - $short{$res}\n";
					#say "aa -$res_num-";
					#L1
					if ($res_num eq 24){
						#say "$res - $short{$res}";
						$L1_start =$new_res_num;
						$L1_control=1;
					}
					if($L1_control){
						#say "$res - $short{$res}";
						my $new =$short{$res};
						$L1_seq = $L1_seq.$new;
					}
					if ($res_num eq 34){
						$L1_end =$new_res_num;
						$L1_control=0;
					}
					
					#L2
					if ($res_num eq 50){
						$L2_start =$new_res_num;
						$L2_control=1;
					}
					if($L2_control){
						#say "$res - $short{$res}";
						my $new =$short{$res};
						$L2_seq = $L2_seq.$new;
					}
					if ($res_num eq 56){
						$L2_end =$new_res_num;
						$L2_control=0;
					}
					
					#L3
					if ($res_num eq 89){
						$L3_start =$new_res_num;
						$L3_control=1;
					}
					if($L3_control){
						#say "$res - $short{$res}";
						my $new =$short{$res};
						$L3_seq = $L3_seq.$new;
					}
					if ($res_num eq 97){
						$L3_end =$new_res_num;
						$L3_control=0;
					}
					
					
					
				}
				elsif  ($chain eq "H"){
					#H1
					if ($res_num eq 26){
						$H1_start =$new_res_num;
						$H1_control=1;
					}
					if($H1_control){
						#say "$res - $short{$res}";
						my $new =$short{$res};
						$H1_seq = $H1_seq.$new;
					}
					if ($res_num eq 32){
						$H1_end =$new_res_num;
						$H1_control=0;
					}
					
					#H2
					if ($res_num eq 52){
						$H2_start =$new_res_num;
						$H2_control=1;
					}
					if($H2_control){
						#say "$res - $short{$res}";
						my $new =$short{$res};
						$H2_seq = $H2_seq.$new;
					}
					if ($res_num eq 56){
						$H2_end =$new_res_num;
						$H2_control=0;
					}
					
					#H3
					if ($res_num eq 93){
						$H3_start =$new_res_num;
						$H3_control=1;
					}
					if($H3_control){
						#say "$res - $short{$res}";
						my $new =$short{$res};
						$H3_seq = $H3_seq.$new;
					}
					if ($res_num eq 102){
						$H3_end =$new_res_num;
						$H3_control=0;
					}
				}
			}
			
			if (!($alt eq " ")){
				say "NOW - $alt";
				if ($alt eq "A"){
					$pre_alt=$alt;
				}
				else {
					if ($pre_alt eq "A"){
						$pre_alt=$alt;
					}
					else{
						#$pre_alt=$alt;
						$format->{altloc} = " ";
					}
				}
				say "ALT - $format->{altloc}";
			}
			
			
			$format->{resseq}=$new_res_num;
			
			printf OUTPUT2 "%6s%5s %4s%s%3s %s%4s%s   %8s%8s%8s%6s%6s\n",
	      $format->{tag},
	      $format->{seq},
	      $format->{name},
	      $format->{altloc},
	      $format->{residueId},
	      $format->{chain},
	      $format->{resseq},
	      " ",
	      $format->{x},
	      $format->{y},
	      $format->{z},
	      $format->{occ},
	      $format->{Bfact};	
			
			
		}
		else {
			say OUTPUT2 $line;
			if ($line =~ /TER/){
				#say $line;
				$cur_res = -1;
				$new_res_num=0;
				}
		}
		#say $line;
	}
	
	
	print "L1: $L1_start - $L1_end - seq: $L1_seq\n";
	print "L2: $L2_start - $L2_end - seq: $L2_seq\n";
	print "L3: $L3_start - $L3_end - seq: $L3_seq\n";
	print "H1: $H1_start - $H1_end - seq: $H1_seq\n";
	print "H2: $H2_start - $H2_end - seq: $H2_seq\n";
	print "H3: $H3_start - $H3_end - seq: $H3_seq\n";
	
	my $len_L1= length($L1_seq);
	my $len_L2= length($L2_seq);
	my $len_L3= length($L3_seq);
	my $len_H1= length($H1_seq);
	my $len_H2= length($H2_seq);
	my $len_H3= length($H3_seq);
	
	printf OUTPUT  "%10s  %3d  %3d   L   %20s  %2d  L1\n", $newpdb , $L1_start, $L1_end, $L1_seq, $len_L1;
	printf OUTPUT  "%10s  %3d  %3d   L   %20s  %2d  L2\n", $newpdb , $L2_start, $L2_end, $L2_seq, $len_L2;
	printf OUTPUT  "%10s  %3d  %3d   L   %20s  %2d  L3\n", $newpdb , $L3_start, $L3_end, $L3_seq, $len_L3;
	printf OUTPUT  "%10s  %3d  %3d   H   %20s  %2d  H1\n", $newpdb , $H1_start, $H1_end, $H1_seq, $len_H1;
	printf OUTPUT  "%10s  %3d  %3d   H   %20s  %2d  H2\n", $newpdb , $H2_start, $H2_end, $H2_seq, $len_H2;
	printf OUTPUT  "%10s  %3d  %3d   H   %20s  %2d  H3\n", $newpdb , $H3_start, $H3_end, $H3_seq, $len_H3;

	close (PDB);
	close (OUTPUT2);
}

close (OUTPUT);






###########################
#     FUNCTION
###########################

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
        $newAt->{resseq}=int(substr($line,22,4));
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
	      " ",
	      $_->{x},
	      $_->{y},
	      $_->{z},
	      $_->{occ},
	      $_->{Bfact};
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
	      " ",
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

