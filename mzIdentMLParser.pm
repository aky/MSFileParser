####################################################################################################
#            Perl Module to Parse mzIDentML Version-1.1                                                                                                          #
#            Written by: Puneet Kumar Kadimi, Research Applicationn Developer DDRC-THSTI                                            #
#            Date : 18-January-2017                                                                                                                                               #
####################################################################################################

use strict;
use warnings;
use XML::Twig;
use Data::Dumper;

sub ParsemzID
{
	my ($Infile,$outfile)=@_;
	my ($PSMFile,$InferenceFile) = ReadAnymzID($Infile,$outfile);
	return ($PSMFile,$InferenceFile);
}

#my $ParsedFile=ReadAnymzID('./INPUT-QTOF-18Mix/InferenceonMSGFTandem/msgf/NEW_QTOF_all.MSGFPlus.FDR-0.01.Groups.mzid'); # NEW_QTOF_all_Comet.mzid Crux-QTOF-tide-search.mzid NEW_QTOF_all.MSGFPlus.mzid NEW_QTOF_all_Myrimatch.mzid NEW_QTOF_all.tandem.2016_11_21_15_15_33.t.mzid NEW_QTOF_all.OMSSA.pep.mzid
#my $ParsedFile=ReadAnymzID('./INPUT-QTOF-18Mix/InferenceonMSGFTandem/msgf/NEW_QTOF_all.MSGFPlus.FDR-0.01.Groups.mzid'); # ./OtherFiles/55merge_mascot_full.mzid
# tandem/NEW_QTOF_all.tandem.2016_11_21_15_15_33.t.FDR0.01.GroupFile.mzid
sub ReadAnymzID{

	my ($mzIDFile,$OutFile)=@_;
	
	my $OutFileB=$mzIDFile;
    #$OutFile=~s/mzid$/tsv/ig;
    $OutFileB=~s/mzid$/inference.txt/ig;
	open(OUT,">$OutFile") or die "Cannot open $OutFile\n";
    open(INFR,">$OutFileB") or die "Cannot open $OutFile\n";
	my $t= XML::Twig->new(twig_roots=>{'MzIdentML'=>\&RetriveFrommzID}); # give the root tag of the file 
	$t->parsefile($mzIDFile) or print("File not found\n");
	close OUT;
    close INFR;
	$t->purge;
	return($OutFile,$OutFileB);
}

sub RetriveFrommzID
{
    my ($twig,$ele)=@_;
	
    my $r=$ele;
    
	my %Data;
    my %InferanceData;
    my %ModHash;
    my %DBSeqHash;
	my %PeptideHash;
	my %PepEvedIDHash;
    my %PepEvedHash;
    my %SpectrumIDHash;
	
	my $Peprefid;
    my $spectrumIDindex;
    my $spectrumSIRidRef;
    my $AlgoName;
    
	my @Spectrumid;
    my @PassThreshold;
    my @Rank;
    my @Peptideref;
    my @PeptideEveRef;
    my @Calcmass;
    my @Expermass;
    my @Charge;
    my @QualityScore;
    my @Mods;
    my @SpectrumSIIRef;
    
	my $PeptideIDA;
    my $PeptideIDB;
    my $PeptideSequence;
    my $Spectrum;
    my $ConCatPepEveRef;

    my $Flag=0;
    my $ModFlag=0;
    my $ModEndFlag=0;
    my $ProtInfFlag=0;
    my $protambgrp;
    my $ProtGrpFlaf=0;
    my $PepEveRefFlag=0;
    my $ScoreCounter=0;
    my $AlgoNameFlag=0;

    my %AlgoScrndHdr=(
        'MS-GF+'=>['6',['RawScore','DeNovoScore','SpecEValue','EValue','QValue','PepQValue']],
        'OMSSA'=>['4',['P-Value','E-Value','Matached-Peaks','Unmatched-Peaks']],
        'xtandem'=>['2',['E-Value','Hyperscore']],
        'Comet'=>['8',['Matched-Peaks','Unmatched-Peaks','Xcorr','Deltacn','Deltacnstar','Spscore','Sprank','E-Value']],
        'MyriMatch'=>['4',['Matched-Peaks','Unmatched-Peaks','MVH','mzFidelity']],
        'Mascot'=>['2',['E-Value','Score']], # This name is kept as optional
		'Mascot Server'=>['2',['E-Value','Score']],
        'Crux'=>['2',['X-Corr','Deltacn']],
        'Sequest'=>['0',['None']],
    );

    while($r=$r->next_elt())
	{
        # Search engine name
        if(($r->gi eq 'AnalysisSoftware')) # this checks if the tags arrives
		{
            if($AlgoNameFlag == 0)
            {
                 # This gets Algorithm name only once; to prevent overwriting 
                 # in case pepxml to mzid conversion external software also appends own name with analysissoftware tag additionally
                $AlgoName = $r->{'att'}->{'name'};
                # chekc here if this algo exists then only proceed else die
               if($AlgoName)
               {
                    unless(exists $AlgoScrndHdr{$AlgoName})
                    {
                        die('ERROR: Unknown algo, programme exits...\n');
                    }
               }
               else
               {
                   die('ERROR: No algo name exists programme exits...\n');
               }

                $AlgoNameFlag++;
            }
        }
        
        # <DBSequence length="1668" searchDatabase_ref="SearchDB_1" accession="Rnd2psu|NC_LIV_145550" id="dbseq_Rnd2psu|NC_LIV_145550">
        # Read above line and store as id => Accession in %DBSeqHash

        if(($r->gi eq 'DBSequence')) # this checks if the tags arrives 
		{
            $DBSeqHash{$r->{'att'}->{'id'}}=$r->{'att'}->{'accession'};
        }

=start
<Peptide id="EGNTLVIVTADHAHASQIVAPDTK___">
	<PeptideSequence>EGNTLVIVTADHAHASQIVAPDTK</PeptideSequence>
</Peptide>
=cut
			
			# This section stores peptides information
			if($r->gi eq 'Peptide')
			{
					$Peprefid=$r->{'att'}->{'id'};
			}
			if($r->gi eq 'PeptideSequence')
            {
					#$PepSequence=$r->text;
					$PeptideHash{$Peprefid}=$r->text;
            }

        # Starts reading peptides and the modifications  
        if(($r->gi eq 'Peptide') or ($r->gi eq 'PeptideEvidence')) # this checks if the peptide tags arrives
		{

            # <Peptide id="PEP_1">
            # <PeptideSequence>YICDNQDTISSK</PeptideSequence>
            
            if($r->gi eq 'PeptideEvidence') 
            {
                $ModEndFlag++;
            }

            # This gets tags attributes information
            $ModFlag++;

            if($ModFlag == 2)
            {
                $PeptideIDB=$PeptideIDA;
                $PeptideIDA = $r->{'att'}->{'id'};
            }
            else
            {
                $PeptideIDA = $r->{'att'}->{'id'};
            }
        }

        # <PeptideEvidence dBSequence_ref="dbseq_Rnd3psu|NC_LIV_083320" peptide_ref="RVDSGLHCPLLPDDR" start="62" end="76" pre="K" post="A" isDecoy="true" id="PE1_2_0">
        # Read Above line and store as id => dBSequence_ref %PepEvedHash
        
        if(($r->gi eq 'PeptideEvidence')) # this checks if the tags arrives
		{
            $PepEvedHash{$r->{'att'}->{'id'}}=$r->{'att'}->{'dBSequence_ref'};
        }

        if($ModFlag == 1)
        {
            if($r->gi eq 'PeptideSequence')
            {
               $PeptideSequence=$r->text;
            }

            # Start storing modification TODO: Still working on fetching modifications not finished yet

            if($r->gi eq 'Modification')
            {
                # checking here if Amino Acid and Position is also available here then store
            
                my $VarMod= $r->{'att'}->{'monoisotopicMassDelta'};

                if($r->{'att'}->{'location'})
                {
                    $VarMod.=':'.$r->{'att'}->{'location'};
                }

                if($r->{'att'}->{'residues'})
                {
                    $VarMod.=':'.$r->{'att'}->{'residues'};
                }

                push(@Mods,$VarMod);
            }
        }
        else
        {
            # Print the stored mods , null array and set counter to 1
             if(scalar(@Mods) == 0)
             {
                 push(@Mods,'None'); # in case there is no modification
             }

             if($ModEndFlag <= 1)
             {
                 # pushing peptide here at position 0 later to Retrive

                if($PeptideIDB && scalar(@Mods) > 0)
                {
                    # HAVING A STRUCTURE LIKE Hash{pepid}['PEPTIDESEQUENCE','M1|M2|M3|M4']

                    $ModHash{$PeptideIDB}[0]=$PeptideSequence;
                    $ModHash{$PeptideIDB}[1]=join('|',@Mods);
                }
             }
             else
             {
                 # Do nothing , if i will do next , parser will skip from all related information which are coming next
				 # if it arrives store; id => peptide_ref
             }

			 if($ModEndFlag > 0)
			{
				$PepEvedIDHash{$r->{'att'}->{'id'}}=$r->{'att'}->{'peptide_ref'};
			}

             # Releasing variables data here
             ($PeptideSequence,$PeptideIDB,@Mods)=();
             $ModFlag=1;
        }

        # <SpectrumIdentificationResult id="SIR_1" name="NEW_QTOF_all.00001.00001" spectrumID="index=0" spectraData_ref="SD">
        
        if($r->gi eq 'SpectrumIdentificationResult') # this checks if the tags arrives
		{
            # When the next spectrum identification result arrives

            if($Flag  == 1) # if its is 1 means the data has been filled
            {
            
                my @PassThreshold_c=@PassThreshold;
                my @Rank_c=@Rank;
                my @Peptideref_c=@Peptideref;
                my @PeptideEveRef_c=@PeptideEveRef;
                my @Calcmass_c=@Calcmass;
                my @Expermass_c=@Expermass;
                my @Charge_c=@Charge;
                my @QualityScore_c=@QualityScore;
                my @SpectrumSIIRef_c=@SpectrumSIIRef;
                #print @SpectrumSIIRef_c,"\n"; <>;
                $Data{$Spectrum}=[\@Peptideref_c,\@PeptideEveRef_c,\@Calcmass_c,\@Expermass_c,\@Charge_c,\@PassThreshold_c,\@QualityScore_c,\@Rank_c,\@SpectrumSIIRef_c];            
                 # print Dumper(\%Data); <>;
                # Releasing data from Spectrum var, all arrays and setting flag to 0
                ($Spectrum,@Peptideref,@PeptideEveRef,@Calcmass,@Expermass,@Charge,@PassThreshold,@QualityScore,@Rank,$Flag,@SpectrumSIIRef)=();
            }


            # This section checks if the main tag has name which contains spcetrum then assigns it 
            # else for safty it will assign the id information# in case if no spectrum id is found it will use the id
            # as spectrum id
            # %SpectrumIDHash

            if($r->{'att'}->{'name'})
            {
                $Spectrum=$r->{'att'}->{'name'};
            }
            else
            {
                 $Spectrum=$r->{'att'}->{'id'};
            }

            # Storing Spectrum index
            if($r->{'att'}->{'spectrumID'})
            {
                $spectrumIDindex=$r->{'att'}->{'spectrumID'};
            }

            # Storing spectrum refrence id
            $spectrumSIRidRef=$r->{'att'}->{'id'};

            $Flag ++;  
        }

 # <SpectrumIdentificationItem id="SII_143_2" calculatedMassToCharge="997.651655" chargeState="1" experimentalMassToCharge="997.582" peptide_ref="VLRSGKPLK_00000000000" rank="2" passThreshold="false">

        if($r->gi eq 'SpectrumIdentificationItem') # this checks if the tags arrives
		{ 
            push(@PassThreshold,$r->{'att'}->{'passThreshold'});
            push(@Rank,$r->{'att'}->{'rank'});
            push(@Peptideref,$r->{'att'}->{'peptide_ref'});
            push(@Calcmass,$r->{'att'}->{'calculatedMassToCharge'});
            push(@Expermass,$r->{'att'}->{'experimentalMassToCharge'});
            push(@Charge,$r->{'att'}->{'chargeState'}); 
            push(@SpectrumSIIRef,$r->{'att'}->{'id'});
            $PepEveRefFlag=0;

            # Set here score counting flag == 0
            $ScoreCounter=0;
        }
        
        # Storing PeptideEvidenceRef information <Protein id>
        # <PeptideEvidenceRef peptideEvidence_ref="[Contaminant]SW:K1CX_HUMAN_PEP_3919"/> # @PeptideEveRef

        if($r->gi eq 'PeptideEvidenceRef') # this checks if the tags arrives
		{
            if($r->{'att'}->{'peptideEvidence_ref'})
            {
                # Concatenate here all the pepevdence refs in a variable by separator
                # when it reaches to spectrum identuification item push in array and null the variable
                # this steps would to be done because mascot gives more than one peptideevidenceref for a given spectrum
                $ConCatPepEveRef.=$r->{'att'}->{'peptideEvidence_ref'}.'*';
            }
            else
            {
                # in case there is nothing found NA will be inserted
                $ConCatPepEveRef.='NA|';
            }

        }

        if($r->gi eq 'cvParam') # this checks if the tags arrives
		{
                if ($Flag == 1)
                {
                       # Storing all scores upto here storing if the score exists
                       if(defined $r->{'att'}->{'value'}) # if value is not null
                       {
                            if ($r->{'att'}->{'name'}=~m/spectrum\s+title/ig) # Now it checks if spectrum is give some where between scores
                            {
                                $Spectrum=$r->{'att'}->{'value'}; # storing here the scan id
                            }
                            else #if not means its score pushes in array
                            {
                                # If it founds any score is zero like in then pushes character "0.000"case @ <cvParam cvRef="MS" accession="MS:1001362" name="number of unmatched peaks" value="0"/>
                                # retrive algo score counts by sending  algo name
                                
                                $ScoreCounter++;
                                if($ScoreCounter <= $AlgoScrndHdr{$AlgoName}[0]) # Here it checks that if counter reaches till the defined no of scores
                                {
                                    if($r->{'att'}->{'value'} != 0)
                                    {
                                        push(@QualityScore,$r->{'att'}->{'value'});
                                    }
                                    else
                                    {
                                        push(@QualityScore,'0.0000');
                                    }
                                }
                            }
                            
                            # Finally storing spectrum index here [0] spectrum index [1]
                            $SpectrumIDHash{$Spectrum}[0]=$spectrumIDindex;
                            $SpectrumIDHash{$Spectrum}[1]=$spectrumSIRidRef;
                       }

                       # push here the pepevedencerefs
                       if($PepEveRefFlag == 0)
                       {
                            push(@PeptideEveRef,$ConCatPepEveRef);
                            $ConCatPepEveRef=();
                            $PepEveRefFlag++;
                       }
                }

                # print Dumper(\%Data); <>;
        }

        #------------------------ Starts storing Protein Inference Information ------------------------------# #
        
         # Protein Inferance Information
        if($r->gi eq 'ProteinDetectionList' && $r->{'att'}->{'id'} eq 'PDL_1') # this checks if the tags arrives
		{
            $ProtInfFlag = 1;
        }

        if($ProtInfFlag == 1)
        {
            # If Protein Group name arrives
            if($r->gi eq 'ProteinAmbiguityGroup')
            {
                if($r->{'att'}->{'id'})
                {
                    my @Data;
                    $protambgrp=$r->{'att'}->{'id'};
                    $InferanceData{$protambgrp}=\@Data;
                }
            }

            # if Protein id arrives which forms group
            if($r->gi eq 'ProteinDetectionHypothesis')
            {
                push(@{$InferanceData{$protambgrp}},'ProteinGroup:'.$r->{'att'}->{'id'}.':'.$r->{'att'}->{'dBSequence_ref'}.':'.$r->{'att'}->{'passThreshold'});
            }

            # if peptide or spectrum arrives arrives
            if($r->gi eq 'PeptideHypothesis' or $r->gi eq 'SpectrumIdentificationItemRef')
            {
				if($r->gi eq 'PeptideHypothesis')
				{
					push(@{$InferanceData{$protambgrp}},'Peptide:'.$r->{'att'}->{'peptideEvidence_ref'});
				}
				else
				{
					push(@{$InferanceData{$protambgrp}},'Spectrum:'.$r->{'att'}->{'spectrumIdentificationItem_ref'});
				}
            }
			
            if($r->gi eq 'cvParam')
            {
				if($r->{'att'}->{'name'} && $r->{'att'}->{'value'}) # anchor protein , in case there is no value skips
				{
					  push(@{$InferanceData{$protambgrp}},'Score:'.$r->{'att'}->{'name'}.':'.$r->{'att'}->{'value'});
				}
				else
				{
					# print Dumper($r->{'att'}->{'name'}); <>;
				}
              
            }
			else
			{
				# push(@{$InferanceData{$protambgrp}},'NA');
			}
        }
    } # eof next elt
    
    #--One Final Push here, because in above loop allways last spectrum information will be left, so filling here that information--#
    if($Flag  == 1) # if its is 1 means the data has been filled
    {            
        my @PassThreshold_c=@PassThreshold;
        my @Rank_c=@Rank;
        my @Peptideref_c=@Peptideref;
        my @PeptideEveRef_c=@PeptideEveRef;
        my @Calcmass_c=@Calcmass;
        my @Expermass_c=@Expermass;
        my @Charge_c=@Charge;
        my @QualityScore_c=@QualityScore;
        my @SpectrumSIIRef_c=@SpectrumSIIRef;
        $Data{$Spectrum}=[\@Peptideref_c,\@PeptideEveRef_c,\@Calcmass_c,\@Expermass_c,\@Charge_c,\@PassThreshold_c,\@QualityScore_c,\@Rank_c,\@SpectrumSIIRef_c];
        ($Spectrum,@Peptideref,@PeptideEveRef,@Calcmass,@Expermass,@Charge,@PassThreshold,@QualityScore,@Rank,$Flag,@SpectrumSIIRef)=();
        #print Dumper(\%Data); exit;
    } # next_elt loop ends here

    #------------------------#
    # Create a Hash which stores Algo names in key and Heder for scores and no of headers at [0] [1]

    my $LoopFactor=0;

    if($AlgoName)
    {
        if($AlgoName=~m/OMSSA/ig)
        {
            print OUT "ScanID\tSpectrumIndex\tSIR-ID\tSII-ID\tPeptide\tModification\tProtein\tCalculated-Mass\tExperimental-Mass\tCharge\tPass-Threshold\t".join("\t",@{$AlgoScrndHdr{'OMSSA'}[1]})."\tRank\n";
            $LoopFactor=$AlgoScrndHdr{'OMSSA'}[0];
        }
        elsif($AlgoName=~m/crux/ig)
        {
            print OUT "ScanID\tSpectrumIndex\tSIR-ID\tSII-ID\tPeptide\tModification\tProtein\tCalculated-Mass\tExperimental-Mass\tCharge\tPass-Threshold\t".join("\t",@{$AlgoScrndHdr{'Crux'}[1]})."\tRank\n";
            $LoopFactor=$AlgoScrndHdr{'Crux'}[0];
        }
        elsif($AlgoName=~m/tandem/ig)
        {
            print OUT "ScanID\tSpectrumIndex\tSIR-ID\tSII-ID\tPeptide\tModification\tProtein\tCalculated-Mass\tExperimental-Mass\tCharge\tPass-Threshold\t".join("\t",@{$AlgoScrndHdr{'xtandem'}[1]})."\tRank\n";
            $LoopFactor=$AlgoScrndHdr{'xtandem'}[0];
        }
        elsif($AlgoName=~m/mascot/ig)
        {
            print OUT "ScanID\tSpectrumIndex\tSIR-ID\tSII-ID\tPeptide\tModification\tProtein\tCalculated-Mass\tExperimental-Mass\tCharge\tPass-Threshold\t".join("\t",@{$AlgoScrndHdr{'Mascot'}[1]})."\tRank\n";
            $LoopFactor=$AlgoScrndHdr{'Mascot'}[0];
        }
        elsif($AlgoName=~m/sequest/ig)
        {
            print OUT "ScanID\tSpectrumIndex\tSIR-ID\tSII-ID\tPeptide\tModification\tProtein\tCalculated-Mass\tExperimental-Mass\tCharge\tPass-Threshold\t".join("\t",@{$AlgoScrndHdr{'Sequest'}[1]})."\tRank\n";
            $LoopFactor=$AlgoScrndHdr{'Sequest'}[0];
        }
        elsif($AlgoName=~m/comet/ig)
        {
            print OUT "ScanID\tSpectrumIndex\tSIR-ID\tSII-ID\tPeptide\tModification\tProtein\tCalculated-Mass\tExperimental-Mass\tCharge\tPass-Threshold\t".join("\t",@{$AlgoScrndHdr{'Comet'}[1]})."\tRank\n";
            $LoopFactor=$AlgoScrndHdr{'Comet'}[0];
        }
        elsif($AlgoName=~m/myrimatch/ig)
        {
            print OUT "ScanID\tSpectrumIndex\tSIR-ID\tSII-ID\tPeptide\tModification\tProtein\tCalculated-Mass\tExperimental-Mass\tCharge\tPass-Threshold\t".join("\t",@{$AlgoScrndHdr{'MyriMatch'}[1]})."\tRank\n";
            $LoopFactor=$AlgoScrndHdr{'MyriMatch'}[0];
        }
        elsif($AlgoName=~m/ms-gf+/ig)
        {
            print OUT "ScanID\tSpectrumIndex\tSIR-ID\tSII-ID\tPeptide\tModification\tProtein\tCalculated-Mass\tExperimental-Mass\tCharge\tPass-Threshold\t".join("\t",@{$AlgoScrndHdr{'MS-GF+'}[1]})."\tRank\n";
            $LoopFactor=$AlgoScrndHdr{'MS-GF+'}[0];
        }
        else
        {
            die("Header for Algo $AlgoName is not available kindly report to developer\n");
        }
    }
    else
    {
		if(!$AlgoName)
        {
            # For now keep  Loop factor to 2 for crux as crux has no name in header Later change it to 0 
            print STDERR "WARNING: No Algo name is defined result will be printed without header and scores\n";
            print STDERR "INFO: For now temporarily loop factor is set to 2 for crux results later will be ommited\n";
            $LoopFactor=2;
        }
    }

    # Start Iterating the hash here and print the data according to the header: Brute-Force :P
    # $Data{$ScanId}=[\@Peptideref_c,\@Calcmass_c,\@Expermass_c,\@Charge_c,\@PassThreshold_c,\@QualityScore_c,\@Rank_c];
    # print Dumper(\%Data); exit;

    foreach my $Scan (keys %Data)
    {
        my $j=0;
        my $LoopTurns=0;

        for(my $i=0; $i < scalar(@{${$Data{$Scan}}[0]}); $i++)
        {
                # Split b y : and run foreach concatenate protein id and print

                my @PepEvedRefARR=split(/\*/,${${$Data{$Scan}}[1]}[$i]);
                my $PepEvedRefConc;
				
				if(scalar @PepEvedRefARR > 1)
				{
					foreach(@PepEvedRefARR)
					{
							#print $_; <>;
						 chomp($_);
						 $PepEvedRefConc.=$DBSeqHash{$PepEvedHash{$_}}.';';
					}
				}
				else
				{
					$PepEvedRefConc=$DBSeqHash{$PepEvedHash{$PepEvedRefARR[0]}};
					chomp($PepEvedRefConc);
				}

                print OUT "$Scan\t$SpectrumIDHash{$Scan}[0]"; # Scan-ID & Spectrum Index
                print OUT "\t$SpectrumIDHash{$Scan}[1]"; # SpectrumIdentificationResult ID
                print OUT "\t${${$Data{$Scan}}[8]}[$i]"; # SpectrumIdentificationItem ID
                print OUT "\t$ModHash{${${$Data{$Scan}}[0]}[$i]}[0]"; # peptide
                print OUT "\t$ModHash{${${$Data{$Scan}}[0]}[$i]}[1]"; # Modification
                print OUT "\t$PepEvedRefConc";# Protein IDs
                print OUT "\t${${$Data{$Scan}}[2]}[$i]"; # Calculates Mass
                print OUT "\t${${$Data{$Scan}}[3]}[$i]"; # Experimental-Mass
                print OUT "\t${${$Data{$Scan}}[4]}[$i]"; # Charge
                print OUT "\t${${$Data{$Scan}}[5]}[$i]\t"; # Pass Threshold

                # A loop to print all of the score variants for current spectrum here
                $LoopTurns=($LoopTurns+$LoopFactor);

                while($j <  $LoopTurns)
                {
                    if(defined ${${$Data{$Scan}}[6]}[$j])
                    {   
                        # In case no scores values is svailable then print in else part
                        print OUT ${${$Data{$Scan}}[6]}[$j]."\t";
                    }
                    else
                    {
                        print OUT "NA\t";
                    }

                    $j++;                   
                }

                $j=$LoopTurns;
                 
                # Printing Rank here
                if(defined ${${$Data{$Scan}}[7]}[$i])
                {
                    print OUT ${${$Data{$Scan}}[7]}[$i]."\n";
                }
                else
                {
                    print "NA\n";
                }
        }

    }
	
	#print Dumper(\%InferanceData);
=start
 'PAG_82' => [
             'ProteinGroup:PDH_57:dbseq_DECOY_gi|1574043|gb|AAC22672.1|:true',
             'Peptide:PE4441_2_332',
             'Spectrum:SII_4441_1',
             'Score:distinct peptide sequences:1',
             'Score:ProteoGrouper:PDH score:0.0026860059092130005',
             'Score:protein group passes threshold:true',
             'Score:ProteoGrouper:PAG score:0.0026860059092130005'
=cut

    if(scalar(keys (%InferanceData)) > 0)
    {
        foreach(keys %InferanceData)
        {
			my $ScoreFlag=0;

			if($ProtGrpFlaf == 0)
				{
					print INFR "$_\n";
					$ProtGrpFlaf++;
				}
				else
				{
					print INFR "\n$_\n";
				}

            for(my $i=0; $i < scalar @{$InferanceData{$_}}; $i++)
            {
                if(${$InferanceData{$_}}[$i]=~m/ProteinGroup/g)
                {
                    my @Line=split(':',${$InferanceData{$_}}[$i]);
					print INFR "\n$Line[1]\t$DBSeqHash{$Line[2]}\t$Line[3]\t"; 
                }
                elsif(${$InferanceData{$_}}[$i]=~m/Peptide/g)
                {
                    my @Line=split(':',${$InferanceData{$_}}[$i]);
					print INFR '*'.$PeptideHash{$PepEvedIDHash{$Line[1]}}.':';
                }
                elsif(${$InferanceData{$_}}[$i]=~m/Spectrum/g)
                {
					# Printing spectrum later to get their counts for spectral counting
					my @Line=split(':',${$InferanceData{$_}}[$i]);
					print INFR $Line[1].'|';
                }

				elsif(${$InferanceData{$_}}[$i]=~m/Score/g)
                {
					if($ScoreFlag == 0)
					{
						print INFR '*';
						$ScoreFlag++;
					}
				}
            }
        }
    }
    else
    {
        print INFR "#----No Protrein Inference Information is found----#";
    }
}
1;
#------------------------ Programme Ends Here ---------------------#