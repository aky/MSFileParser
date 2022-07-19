use strict;
use warnings;
use XML::Twig;
#use Data::Dumper;
#############################################################################################
# Parse pepXML file as ProteoStats data structure
###################################################################################


sub ReadAnyPepXML{
	#PEPXML version supported 
	my ($pepXMLfile)=@_;
	#print"reading $pepXMLfile...\n";
	my $OutFile=$pepXMLfile;
	$OutFile=~s/\.xml/\.tsv/ig;
	open(OUT,">$OutFile") or die "Cannot open $OutFile\n";
	my $t= XML::Twig->new(twig_roots=>{'msms_pipeline_analysis'=>\&Retrive_pepXML}); # give the root tag of the file 
	$t->parsefile($pepXMLfile) or print("File not found\n");
	close OUT;
	$t->purge;
	#$ArrRef=mapcol($OutFile,$ArrRef);
	return($OutFile);
}
sub Retrive_pepXML
{
	my ($twig,$ele)=@_;
	my $r=$ele;
	# my @Score;
	my $ScanID;
	my $ExpMass;
	my $Charge;
	my $Peptide;
	my $Rank;
	my $Protein;
	my $CalcMass;
	my $DeltaMass;
	my $MissedCleav;
	my $Counter=0;
	#my $RankFlag=0;
	my $headerflag=2; #2=search for algo name, 0=set headers, 1=set values
	my @ModPosition;
	my @ModMass;
	
	my $ModAA='NA';
	my $ModificationType="NA";
	my $SpectrumIndex;
	my $MassDiff;
	
	my $ionscore;
	my $identityscore;
	my $star;
	my $homologyscore;
	my $expect;
	
	my %Scores;

	my $RetentionTime;
	
	my $SearchHitFlag=0;
	
	my $No_tol_term;
	my $Algo;
	
	while($r=$r->next_elt())
	{
		if(($headerflag==2) && ($r->gi eq 'search_summary'))
		{
			$Algo=$r->{'att'}->{'search_engine'};
			#print "Fisrt \$Algo=$Algo\n";
			$headerflag=0;
			#<>;
		}
		if ($headerflag==0)
		{
			if($Algo=~/Mascot/i) # if Mascot
			{
				print OUT "ScanID\tExpMass\tCharge\tSpectrumIndex\tRetentionTime\tPeptide\tRank\tProtein\tCalcMass\tDeltaMass\tMissedCleav\tNo_tol_term\tModPosition\tModMass\tIonScore\tIdentityScore\tStar\tHomologyScore\tExpect\n";
			}
			elsif($Algo=~/Comet/i) # if Comet
			{
				print OUT "ScanID\tExpMass\tCharge\tSpectrumIndex\tRetentionTime\tPeptide\tRank\tProtein\tCalcMass\tDeltaMass\tMissedCleav\tNo_tol_term\tModPosition\tModMass\tXcorr\tDeltacn\tDeltacnStar\tSpScore\tSpRank\tExpect\n";
			}
			elsif($Algo=~/MyriMatch/i) # if MyriMatch
			{
				print OUT "ScanID\tExpMass\tCharge\tSpectrumIndex\tRetentionTime\tPeptide\tRank\tProtein\tCalcMass\tDeltaMass\tMissedCleav\tNo_tol_term\tModPosition\tModMass\tmvhmz\tFidelity\tXCorr\n";
			}
			else
			{
				#print"Nothing found\n";<>;
				print OUT "ScanID\tExpMass\tCharge\tSpectrumIndex\tRetentionTime\tPeptide\tRank\tProtein\tCalcMass\tDeltaMass\tMissedCleav\tNo_tol_term\tModPosition\tModMass\tExpect\n";
			}
			$headerflag=1;
		}		
		if($r->gi eq 'spectrum_query')
		{
			$ScanID=$r->{'att'}->{'spectrum'};
			$ExpMass=$r->{'att'}->{'precursor_neutral_mass'};
			$Charge=$r->{'att'}->{'assumed_charge'};
			$SpectrumIndex=$r->{'att'}->{'index'};
			$RetentionTime=$r->{'att'}->{'retention_time_sec'}; # Not found in Mascot pepxml file
			
			# in case if no retention time is found NA will be assigned 
			if(!$RetentionTime)
			{
				$RetentionTime='NA';
			}
		}
		if($r->gi eq 'search_hit')
		{
			$Peptide=$r->{'att'}->{'peptide'};
			$Rank=$r->{'att'}->{'hit_rank'};
			$Protein=$r->{'att'}->{'protein'};
			$CalcMass=$r->{'att'}->{'calc_neutral_pep_mass'};
			$DeltaMass=$r->{'att'}->{'massdiff'};
			$MissedCleav=$r->{'att'}->{'num_missed_cleavages'};
			$No_tol_term=$r->{'att'}->{'num_tol_term'};
			#$RankFlag=1;
		}
		if($r->gi eq 'mod_aminoacid_mass')
		{
			push(@ModPosition,$r->{'att'}->{'position'});
			push(@ModMass,$r->{'att'}->{'mass'});
		}
		if($r->gi eq 'search_score')
		{
			my $ScoresRef=\%Scores;
			if($SearchHitFlag==0)
			{
				# Return all the score in tab separated with 1 status that score is finished
				
				my $Finished;
				($Finished,$ScoresRef)=RetriveScorePepXML($Algo,$r,\$Counter,$ScoresRef);
				if($Finished == 1)
				{
					$SearchHitFlag=1;
				}
				%Scores=%$ScoresRef;
			}
			if($SearchHitFlag==1) # print every thing, reset flag, undef all variables
			{
				print OUT"$ScanID\t$ExpMass\t$Charge\t$SpectrumIndex\t$RetentionTime\t$Peptide\t$Rank\t$Protein\t$CalcMass\t$DeltaMass\t$MissedCleav\t$No_tol_term\t",join("|",@ModPosition),"\t",join("|",@ModMass),"\t";
				#print"Checking score hash....\n\n";
				if($Algo=~/^Mascot$/i) # if Mascot
				{
					print OUT $Scores{'ionscore'},"\t";
					print OUT $Scores{'identityscore'},"\t";
					print OUT $Scores{'star'},"\t";
					print OUT $Scores{'homologyscore'},"\t";
					print OUT $Scores{'expect'},"\n";
					
					#print $Scores{'ionscore'},"\t";
					#print $Scores{'identityscore'},"\t";
					#print $Scores{'star'},"\t";
					#print $Scores{'homologyscore'},"\t";
					#print $Scores{'expect'},"\n";
				}
				elsif($Algo=~/Comet/i) # if Comet
				{
					print OUT $Scores{'xcorr'},"\t";
					print OUT $Scores{'deltacn'},"\t";
					print OUT $Scores{'deltacnstar'},"\t";
					print OUT $Scores{'spscore'},"\t";
					print OUT $Scores{'sprank'},"\t";
					print OUT $Scores{'expect'},"\n";
					#
					#print $Scores{'xcorr'},"\t";
					#print $Scores{'deltacn'},"\t";
					#print $Scores{'deltacnstar'},"\t";
					#print $Scores{'spscore'},"\t";
					#print $Scores{'sprank'},"\t";
					#print $Scores{'expect'},"\n";
				}
				elsif($Algo=~/Myrimatch/i)  # if Myrimatch
				{
					print OUT $Scores{'mvh'},"\t";
					print OUT $Scores{'mzFidelity'},"\t";
					print OUT $Scores{'xcorr'},"\n";
					
					#print $Scores{'mvh'},"\t";
					#print $Scores{'mzFidelity'},"\t";
					#print $Scores{'xcorr'},"\n";
				}
				else
				{
					print OUT $Scores{'expect'},"\n";
					#print $Scores{'expect'},"\n";
				}
				# print"Done\tpress enter to continue....\n";<>;
				#,join("|",@ModMass),"\t",$score{'ionscore'},("\t",@Scores),"\n"; # \t$ionscore\t$identityscore\t$star\t$homologyscore\t$expect
				($Peptide,$Rank,$Protein,$CalcMass,$DeltaMass,$MissedCleav,$No_tol_term)=undef;
				undef(@ModPosition);
				undef(@ModMass);
				undef(%Scores);
				$SearchHitFlag=0;
			}
		}
		
		$ele->purge;		
	}
}

sub RetriveScorePepXML
{
	my($Algo,$r,$Counter,$ScoreRef)=@_;
	# my $Score;
	my $Finished=0;
	#print "Begin of RetriveScorePepXML=",scalar keys %{$ScoreRef}; <>;
	if($Algo=~/Mascot/i) # if Mascot
	{
		if ($r->{'att'}->{'name'} eq 'ionscore') 
		{
			#$Score=
			$$ScoreRef{'ionscore'}=$r->{'att'}->{'value'};
			$$Counter++;
		}
		elsif ($r->{'att'}->{'name'} eq 'identityscore') 
		{
			$$ScoreRef{'identityscore'}=$r->{'att'}->{'value'};
			$$Counter++;
		}
		elsif ($r->{'att'}->{'name'} eq 'star') 
		{
			$$ScoreRef{'star'}=$r->{'att'}->{'value'};
			$$Counter++;
		}
		elsif ($r->{'att'}->{'name'} eq 'homologyscore') 
		{
			$$ScoreRef{'homologyscore'}=$r->{'att'}->{'value'};
			$$Counter++;
		}
		elsif ($r->{'att'}->{'name'} eq 'expect') 
		{
			$$ScoreRef{'expect'}=$r->{'att'}->{'value'};
			$$Counter++;
		}
		else
		{
			print STDERR "Unknown Score type in Mascot pepXML.\n";
		}
		
		if ($$Counter == 5)
		{
			# when the mascot has no more score 
			$Finished=1;
			# reference vars value will be assigned to 0
			$$Counter=0;
		}
		
	}
	elsif($Algo=~/Comet/i) # if Comet
	{
		# Xcorr\tDeltacn\tDeltacnStar\tSpScore\tSpRank\tExpect
		if ($r->{'att'}->{'name'} eq 'xcorr')
		{
			$$ScoreRef{'xcorr'}=$r->{'att'}->{'value'};
			$$Counter++;
		}
		elsif ($r->{'att'}->{'name'} eq 'deltacn')
		{
			$$ScoreRef{'deltacn'}=$r->{'att'}->{'value'};
			$$Counter++;
		}
		elsif ($r->{'att'}->{'name'} eq 'deltacnstar')
		{
			$$ScoreRef{'deltacnstar'}=$r->{'att'}->{'value'};
			$$Counter++;
		}
		elsif ($r->{'att'}->{'name'} eq 'spscore')
		{
			$$ScoreRef{'spscore'}=$r->{'att'}->{'value'};
			$$Counter++;
		}
		elsif ($r->{'att'}->{'name'} eq 'sprank')
		{
			$$ScoreRef{'sprank'}=$r->{'att'}->{'value'};
			$$Counter++;
		}
		elsif ($r->{'att'}->{'name'} eq 'expect')
		{
			$$ScoreRef{'expect'}=$r->{'att'}->{'value'};
			$$Counter++;
		}
		else
		{
			print STDERR "Unknown Score type in Comet pepXML.\n";
		}
		
		if ($$Counter == 6)
		{
			# when the mascot has no more score 
			$Finished=1;
			# reference vars value will be assigned to 0
			$$Counter=0;
		}
	}
	elsif($Algo=~/Myrimatch/i)  # if Myrimatch
	{
		if ($r->{'att'}->{'name'} eq 'mvh')
		{
			$$ScoreRef{'mvh'}=$r->{'att'}->{'value'};
			$$Counter++;
		}
		elsif ($r->{'att'}->{'name'} eq 'mzFidelity')
		{
			$$ScoreRef{'mzFidelity'}=$r->{'att'}->{'value'};
			$$Counter++;
		}
		elsif ($r->{'att'}->{'name'} eq 'xcorr')
		{
			$$ScoreRef{'xcorr'}=$r->{'att'}->{'value'};
			$$Counter++;
		}
		else
		{
			print STDERR "Unknown Score type in Myrimatch pepXML.\n";
		}
		
		# checking whether all scores has been finished 
		if ($$Counter == 3)
		{
			# when the mascot has no more score 
			$Finished=1;
			# reference vars value will be assigned to 0
			$$Counter=0;
		}
	}
	else
	{
		if ($r->{'att'}->{'name'} eq 'expect')
		{
			$Finished=1;
			$$Counter=0;
			$$ScoreRef{'expect'}=$r->{'att'}->{'value'};
		}
		else
		{
			print STDERR "Unknown Score type in $Algo pepXML.\n";
		}
		
	}
	#print "Befor rerturning RetriveScorePepXML=",scalar keys %{$ScoreRef},"\n";
	#print "\$Finished=$Finished";<>;

	return($Finished,$ScoreRef);
}

sub mapcolumns
{
	my($pepXML)=@_;
	my $c=0;
	my $Algo;
	my @arr;
	my $scorecol1='NA'; #Score/XCorr etc
	my $scorecol2='NA'; #Evalue
	my $scancol=0; #for all
	my $pepcol=5;
	my $rankcol=6;
	my $masscol=8;
	open(XML, $pepXML) or die $!, " \"$pepXML\" cannot be opened/doesn't exist\n";#|EXPR, LIST) or die EXPR;
	while (<XML>)
	{
		$c++;
		chomp $_;
		if ($_=~m/search_engine=\"([A-Za-z]+)\"/)
		{
			$Algo=$1;
			#print "$Algo\tat=$c\n";
			last;
		}		
	}
	if($Algo=~/Mascot/i) # if Mascot
	{
		$scorecol1=14;
		$scorecol2=18;
	}
	elsif($Algo=~/Myrimatch/i)  # if Myrimatch
	{
		$scorecol1=16;
		#$scorecol2=undef;
	}
	elsif($Algo=~/Comet/i) # if Comet
	{
		$scorecol1=14;
		$scorecol2=19;
	}
	else
	{
		#$scorecol1=14;
		$scorecol2=14;
	}
	#print "Inside mapcol=",	join"\t",($scorecol1,$scorecol2,$scancol,$pepcol,$rankcol,$masscol);<>;

	push(@arr,$scorecol1,$scorecol2,$scancol,$pepcol,$rankcol,$masscol);
	return(@arr);
}

1;