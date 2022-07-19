#!usr/bin/perl
use strict;
use warnings;
use constant proton=>1.00782504;
use Spreadsheet::ParseExcel;
use XML::Twig;
#use Data::Dumper;
################################################################################
# MODULE FOR FILE INPUT & PARSING and OUTPUT
################################################################################
# Based on Algorithm name/type , file parsing will take care of-
# 1. Ranked vs Non-Ranked target-decoy file reading
# 2. Parsing XML type file outputs
# 3. Create Outputs (and XML outputs in future)
################################################################################
my %actions = ( tab => \&split_by_tab,
                comma => \&split_by_comma,
                space => \&split_by_space, 
              );
##################################################

sub split_by_comma{
    my($var)=@_;
    my@arr=split ',',$var;
    return(@arr);
}
sub split_by_space{
    my($var)=@_;
    my@arr=split /\s/,$var;
    return(@arr);
}
sub split_by_tab{
    my($var)=@_;
    my@arr=split /\t/,$var;
    return(@arr);
}

sub ReadTEXT{
    my ($file,$scorecol,$scancol,$pepcol,$rankcol,$separator)=@_;
    #print"reading $file...\n";
    open FILE,$file or die "Cannot open file \"$file\". Check the filename\n";
    my %psm;
       
    my $lc=0; #header flag
    while(<FILE>)
    {
	chomp $_;
	if( ( $lc > 0 ) && ( $_!~/^\s+$/ ) )#no header/blank
            {
                #my @arr1=split /\s/, $_;
		
		my @arr1=$actions{$separator}->($_);
		
                if ( $arr1[$rankcol] == 1 )#rank 1 PSMs(hash removes duplicate scanIDs with same score)
                {
                    if (defined$psm{$arr1[$scancol]})
                    {
                        push ( @{$psm{$arr1[$scancol]}} , uc($arr1[$pepcol]) ); # push peptides if already defined
                    }
                    else
                    {
                        $psm{$arr1[$scancol]}=[ $arr1[$scorecol] , uc($arr1[$pepcol]) ];#scan id=>[ score , peptide]
                    }
                }
            }
        else
	    {
		$lc=1;
	    }
    }
    close FILE;
    #print"reading $file...Complete\n";
    return(\%psm);    
}
#Also saves mass(or the correlated variable for FlexiFDR calculation)
sub ReadFlexiTEXT{
    my ($file,$scorecol,$scancol,$pepcol,$rankcol,$masscol,$separator)=@_;
    open FILE,$file or die "Cannot open file \"$file\". Check the filename\n";
    my %psm;
    print"Inside ReadFlexiTEXT\n",$separator," args found\n";<>;

    #An extra value of Calc.(Theoretical) Mass
    my $lc=0; #header flag
    while(<FILE>)
    {
	chomp $_;
	if( ( $lc > 0 ) && ( $_!~/^\s+$/ ) &&($_ ne '') )#no header/blank
            {
                my @arr1=$actions{$separator}->($_);
                if ( $arr1[$rankcol] == 1 )#rank 1 PSMs(hash removes duplicate scanIDs with same score)
                {
                    if (defined$psm{$arr1[$scancol]})
                    {
                        push ( @{$psm{$arr1[$scancol]}} , uc($arr1[$pepcol]) ); # push peptides if already defined
                    }
                    else
                    {
                        #Mass at index 1(after score)
			$psm{$arr1[$scancol]}=[ $arr1[$scorecol], $arr1[$masscol], uc($arr1[$pepcol]) ];#scan id=>[ score , peptide]
                    }
                }
            }
        else
	{
            $lc=1;
	}
    }
    close FILE;
    #print "Hashes size before ReadFlexiCSV::",scalar(keys%psm),"\n";#scalar(keys%$decref),"\n";<>;

    return(\%psm);    
}

sub ReadChargeFlexiTEXT{
    my ($file,$scorecol,$scancol,$pepcol,$rankcol,$masscol,$zcol,$separator)=@_;
    open FILE,$file or die "Cannot open file \"$file\". Check the filename\n";
    my %psm; # A hash of Hashrefs, each key is a charge state and value is a hashref of %psm as defined earlier
    print "Inside ReadChargeFlexiTEXT File=$file\t";
    #An extra value of Calc.(Theoretical) Mass
    my $lc=0; #header flag
    my $zmin=undef; #need to start min z count from file
    my $zmax=1; #min z of 1 can work fine here
    my $counter=0;
    while(<FILE>)
    {
	chomp $_;
	if( ( $lc > 0 ) && ( $_!~/^\s+$/ ) &&($_ ne '') )#no header/blank
            {
		my @arr1=$actions{$separator}->($_);
		if (!defined($zmin))
		{
		    $zmin=$arr1[$zcol];
		}
                if ( $arr1[$rankcol] == 1 )#rank 1 PSMs(hash removes duplicate scanIDs with same score)
                {
		    $counter++;
		    if(!(exists$psm{$arr1[$zcol]}))
		    {
			my %localhash;
			$psm{$arr1[$zcol]}=\%localhash;
			$localhash{$arr1[$scancol]}=[ $arr1[$scorecol], $arr1[$masscol],$arr1[$zcol], uc($arr1[$pepcol]) ];#scan id=>[ score , mass, charge, peptide]
			$counter++;
		    }
		    else
		    {
			my $hashref=$psm{$arr1[$zcol]};#dereference Z Hashref value to get inner Hashref
			my %localhash=%$hashref;
			$psm{$arr1[$zcol]}=\%localhash;
			if (defined $localhash{$arr1[$scancol]})
			{
			    push ( @{$localhash{$arr1[$scancol]}} , uc($arr1[$pepcol]) ); # push peptides if already defined
			}
			else
			{
			    #Mass at index 1(after score),Charge at 2
			    $localhash{$arr1[$scancol]}=[ $arr1[$scorecol], $arr1[$masscol],$arr1[$zcol], uc($arr1[$pepcol]) ];#scan id=>[ score , mass, charge, peptide]
			    $counter++;
			    if ($arr1[$zcol]<=$zmin)
			    {
				$zmin=$arr1[$zcol];
			    }
			    elsif($arr1[$zcol]>$zmax)
			    {
				$zmax=$arr1[$zcol];
			    }
			}
			$psm{$arr1[$zcol]}=\%localhash; #reference again after assigning
		    }
		    
                }
            }
        else
	{
	    $lc=1;
	}
    }
    close FILE;
    print"Simple counter=$counter\tZHashrefs=",(scalar keys %{$psm{1}})+(scalar keys %{$psm{2}})+(scalar keys %{$psm{3}}),"\nEnter to proceed\n";
    print "Hashes size before ReadFlexiCSV::",scalar(keys%psm),"\n";#scalar(keys%$decref),"\n";<>;
    return($zmin,$zmax,\%psm);    
}

sub ReadTEXT_FWD{
    my ($tar,$scorecol,$scancol,$pepcol,$rankcol,$rank,$separator)=@_;
    open FILE,$tar or die "Cannot open file in ReadCSV_FWD \"$tar\". Check the filename\n";
    my %psm;
    my $lc=0; #header flag
    while(<FILE>)
    {
	chomp $_;
	if( ( $lc > 0 ) && ( $_!~/^\s+$/ ) )#no header/blank
            {
                my @arr1=$actions{$separator}->($_);
                if ( $arr1[$rankcol] == $rank )#$rank (=2 in FWD) PSMs(hash removes duplicate scanIDs with same score)
                {
                    if (defined$psm{$arr1[$scancol]})
                    {
                        push ( @{$psm{$arr1[$scancol]}} , uc($arr1[$pepcol]) ); # push peptides if already defined
                    }
                    else
                    {
                        $psm{$arr1[$scancol]}=[ $arr1[$scorecol] , uc($arr1[$pepcol]) ];#scan id=>[ score , peptide]
                    }
                }
            }
        else
	    {
		$lc=1;
	    }
    }
    close FILE;
    #print"reading $file...Complete\n";
    return(\%psm);   
}


################################################################################

sub Make_OMSSA_CSV{
    my ($file)=@_;
    my $out=$file;
    $out=~s/\.omssa//;
    $out=~s/\.csv/\.omssa_out\.csv/;
   
    my $lc=0; #header flag
    open FILE,$file or die "Cannot open $file\n";
    open OUT,">$out" or die "Cannot open $out\n";
    #print $out," outfile open\n";
    while(<FILE>)
    {
	chomp $_;
	#print "\$_=$_\n";
	if ($_ eq '')#skip header/blank
	{
	    next;
	}
	if( $lc == 1 )
            {
		print OUT "$_,1\n";
            }
        else
	    {
		print OUT "$_,Rank\n";
		$lc=1;
	    }
    }
    close FILE;
    close OUT; 
    #print"reading $file...Complete\n";
    return($out);
}

################################################################################

sub ReadTandemXML{
    ########################################################################################
    #Program to parse PSMs from X!Tandem output file						   #
    ########################################################################################
    # Modified By: Amit 09, August 2010 (Handling target and decoy output files)		   #
    # Written By : Dhirendra Kumar		   #
    # Modified By: Puneet Kumar on 7/6/2011 (Improved for less physical memory consumption)#
    ########################################################################################
    
    my ($file)=@_;
    my $i=0;
    #my $t1=time();
    chomp $file;
    print "Parsing file $file\n";
    my $out_file=$file;
    $out_file=~s/\.t\.xml$/\.tandem\.csv/;
    #print $out_file;
    open(OUT,">$out_file") or print "Cannot open $out_file\n";
    
    #Add calc mass,DeltaM,charge
    
    print OUT"ScanID,Missed Cleavage,Peptide,HyperScore,Expect,Modification,Protein,ExpMass,CalcMass,DeltaMass,Charge,Rank\n";
    
    my $t= XML::Twig->new(twig_roots=>{'group'=>\&retrieve_element});
    $t->parsefile($file) or print("File not found\n");
    $t->purge;
    close OUT;
    #print "DONE\t";
    #my $t2=time();
    return($out_file);
}

sub retrieve_element
{
    my ($twig,$ele)=@_;
    my $r=$ele;
    my ($hyperscore,$mc,$protein_id,$scan_id,$peptide,$expect,$ExpMass,$CalcMass,$Charge,$DeltaMass);
    my $mod='';
    my $protein_flag=0;
    my $peptide_flag=0;
    my $scan_id_flag=0;
    if(defined $r->id)
    {
	$ExpMass=$r->{'att'}->{'mh'}; #<GAML:attribute type="M+H">1346.61</GAML:attribute>

	$Charge=$r->{'att'}->{'z'};
	while($r=$r->next_elt())
	{
	    if(($r->gi eq 'group')&&($r->{'att'}->{'type'} eq 'support'))
	    {
		$scan_id_flag=1;
	    }
	    if(($r->gi eq 'GAML:attribute')&&($r->{'att'}->{'type'} eq 'M+H'))
	    {
		$ExpMass=$r->text;
		$ExpMass-= proton;
	    }
	    if(($scan_id_flag==1)&&($r->gi eq 'note'))
	    {
		$scan_id=$r->text;
	    }
	    if($r->gi eq 'protein')
	    {
		if($protein_flag==0)
		{
		    $protein_id= $r->{'att'}->{'label'};
		}
		else
		{
		    $protein_id.=';'.$r->{'att'}->{'label'};
		}
		$protein_flag++;
	    }
	    if($r->gi eq 'domain')
	    {
		if($peptide_flag==0)
		{
		    $peptide= $r->{'att'}->{'seq'};
		}
		elsif($protein_flag>1)
		{
		    #print $r->{'att'}->{'seq'},"\t";
		}
		else
		{
		    #print "\n\t\t",$r->{'att'}->{'seq'},"\t";
		}
		$CalcMass= $r->{'att'}->{'mh'};
		$CalcMass = $CalcMass - proton;
		$hyperscore=$r->{'att'}->{'hyperscore'};
		$mc=$r->{'att'}->{'missed_cleavages'};
		$expect=$r->{'att'}->{'expect'};
		$CalcMass=$r->{'att'}->{'mh'};
		$DeltaMass=$r->{'att'}->{'delta'};
		$peptide_flag=1;
	    }
	    if($r->gi eq 'aa')
	    {
		$mod.=$r->{'att'}->{'type'}.':'.$r->{'att'}->{'at'}.':'.$r->{'att'}->{'modified'}.";";
	    }
	}
	$scan_id=~s/\,//g;
	$scan_id=~s/\s+$//g;
	print OUT"$scan_id,$mc,$peptide,$hyperscore,$expect,$mod,$protein_id,$ExpMass,$CalcMass,$DeltaMass,$Charge,1\n";
	$ele->purge;
    }
}
################################################################################
sub ReadMascotDat{
    # Run Mascot2XML to create pepXML file
    # Read the pepXML file
}

################################################################################

sub ReadOMSSAXML{
    
}

################################################################################

sub ReadSequestXLS{
    my ($file)=@_;
    my $parser = Spreadsheet::ParseExcel->new();
    my $workbook = $parser->parse("$file");
    my $output_file=$file;
    $output_file=~s/.xls/.tsv/g;
    open(FILE,">$output_file") or die $!;
    # Without this message was : Wide character in print and parsing was not donne properly
    # the data being printed contained strings flagged as containing non-ASCII utf8 characters,
    # and the handle had not been declared to accommodate such data.
    binmode(FILE, ":utf8");
    if ( !defined $workbook )
    {
	die $parser->error(), ".\n";
    }
    for my $worksheet ( $workbook->worksheets())
    {
	my ( $row_min, $row_max ) = $worksheet->row_range();
	my ( $col_min, $col_max ) = $worksheet->col_range();
	for my $row ( $row_min .. $row_max )
	{
	    my @concate;
	    for my $col ( $col_min .. ($col_max-1) )
	    {
		my $cell = $worksheet->get_cell( $row, $col );
		next unless $cell;
		my $val=$cell->value();#for preventing ions matched from turning to dates
		#my $val="\"".$cell->value()."\"";#for preventing ions matched from turning to dates
		push(@concate,$val);
	    }
	    #shift @concate;
	    print FILE join("\t",@concate)."\n";
	}
    }
    close FILE;
    #print "Returning $output_file\n";
    return ($output_file);    
}

################################################################################

sub ReadPepXML{
    #UseReadCSV for Any CSV file format
    my($pepxml)=@_;
    my $out=ReadAnyPepXML($pepxml);
    return ($out);
}

################################################################################
# Filtering FDR results for score based Algos
sub filterFDRPSMsDesc{
    my ($tar,$ROCref,$scorecol,$rankcol,$FDRoutfile,$separator)=@_;
    #print"Inside filterFDRPSMsDesc...Cutoff=$cutoff\nCreating FDRout as =$FDRoutfile\n";<>;
    open TAR,$tar or die;
    my $f=0; #flag for header
    my $i=0; #result counter
    my @tar2d;
    my $tar2dheader;

    while(<TAR>)
    {
	chomp $_;
	
	if($f==1)
	{
	    my @arr=$actions{$separator}->($_);
	    my $j=0;
	    if ($arr[$rankcol]==1)
	    {
		foreach my $val(@arr)
		{
		    $tar2d[$i][$j]=$val;
		    $j++;
		}
		$i++;
	    }
	}
	elsif($f==0)
	{
	    $tar2dheader=$_.",q-value";
	    $tar2dheader=~s/\t/,/g;
	    $f=1;
	    next;
	}
    }
    close TAR;
    @tar2d=sort{$b->[$scorecol]<=>$a->[$scorecol]} @tar2d;
    #print "Size of Tar 2d=",scalar @tar2d,"\n";
    my @ROC=@$ROCref; #2darray with ROC data
    #print "Size of ROC 2d=",scalar @ROC,"\n";<>;

    my $ROC_header=shift@ROC;
    open (FDR,">$FDRoutfile") or die $!;
    print FDR $tar2dheader,"\n";
    my $start=0;
    #my $flag=0;#flag for q-values below FDR thr
    #print"\$tar2d[43][7]=$tar2d[43][7] is culprit\n";<>;
    for (my $r=0;$r<@ROC;$r++)
    {
	#print"Outerloop counter=$r\t\$start=$start\n";<>;
	for(my $i=$start;$i<@tar2d;$i++)
	{
	    #print"\tInner loop counter=$i\n";

	    #print "\t\$ROC[$r][0]=$ROC[$r][0]\t\$tar2d[$i][$scorecol]=$tar2d[$i][$scorecol]\n";
	    if ($tar2d[$i][$scorecol]>$ROC[$r][0])
	    {
		#print"\$tar2d[$i][-1] is assigned a value of $ROC[$r][1]\n";<>;
		$tar2d[$i][scalar(@{$tar2d[$i]})]=$ROC[$r][1];
		
		print FDR shift@{$tar2d[$i]};#first element(scan id)
		foreach(@{$tar2d[$i]})#from next element, put a comma first
		{
		    print FDR ",$_";
		    #print "$_,";
		}
		print FDR "\n";
		$start++;
		#print "\n";<>;
	    }
	    else
	    {
		#print"\$r=$r\t\$i=$i for \$tar2d[$i][$scorecol]=$tar2d[$i][$scorecol]\t\$ROC[$r][0]=$ROC[$r][0]\n";
		$start=$i;
		last;
	    }		    
	}
    }
    #<>;
    
    print FDR "\nROC\n",join",",@$ROC_header,"\n";
    foreach (@ROC)
    {
	print FDR join ",",@{$_},"\n";
    }
    close FDR;
    return $FDRoutfile; 
}

################################################################################
sub filterFDRPSMsAsc{
    my ($tar,$ROCref,$scorecol,$rankcol,$FDRoutfile,$separator)=@_;
    open TAR,$tar or die;
    my $f=0; #flag for header
    my $i=0; #result counter
    my @tar2d;
    my $tar2dheader;
    while(<TAR>)
    {
	chomp $_;
	if($f==1)
	{
	    #my @arr=split",",$_;
	    my @arr=$actions{$separator}->($_);
	    my $j=0;
	    if ($arr[$rankcol]==1)
	    {
		foreach my $val(@arr)
		{
		    $tar2d[$i][$j]=$val;
		    $j++;
		}
		$i++;
	    }
	}
	elsif($f==0)
	{
	    $tar2dheader=$_.",q-value";
	    $tar2dheader=~s/\t/,/g;
	    $f=1;
	    next;
	}
    }
    close TAR;
    @tar2d=sort{$a->[$scorecol]<=>$b->[$scorecol]} @tar2d;
    #print "Size of Tar 2d=",scalar @tar2d,"\n";
    my @ROC=@$ROCref; #2darray with ROC data
    #print "Size of ROC 2d=",scalar @ROC,"\n";<>;

    my $ROC_header=shift@ROC;
    open (FDR,">$FDRoutfile") or die $!;
    print FDR $tar2dheader,"\n";
    my $start=0;
    #my $flag=0;#flag for q-values below FDR thr

    for (my $r=0;$r<@ROC;$r++)
    {
	for(my $i=$start;$i<@tar2d;$i++)
	{
	    if ($tar2d[$i][$scorecol]<$ROC[$r][0])
	    {
		#print"\$tar2d[$i][-1] is assigned a value of $ROC[$r][1]\n";<>;
		$tar2d[$i][scalar(@{$tar2d[$i]})]=$ROC[$r][1];
		
		print FDR shift@{$tar2d[$i]};#first element
		foreach(@{$tar2d[$i]})#from next element, put a comma first
		{
		    print FDR ",$_";
		    #print "$_,";
		}
		print FDR "\n";
		$start++;

		#print "\n";<>;
	    }
	    
	    else
	    {
		#print"last for $tar2d[$i][$scorecol] .... $ROC[$r][0]";<>;
		$start=$i;
		last;
	    }		    
	}
    }
    
    print FDR "\nROC\n",join",",@$ROC_header,"\n";
    foreach (@ROC)
    {
	print FDR join ",",@{$_},"\n";
    }
    close FDR;
    return $FDRoutfile;
}

################################################################################
# Filtering Concat FDR results for score based Algos
sub filterConFDRPSMsDesc{
    my ($tar,$ROCref,$scorecol,$rankcol,$FDRoutfile,$XTar,$scancol,$separator)=@_;
    my %XTars=%$XTar;
    open TAR,$tar or die;
    my $f=0; #flag for header
    my $i=0; #result counter
    my @tar2d;
    my $tar2dheader;
    while(<TAR>)
    {
	chomp $_;
	if($f==1)
	{
	    my @arr=$actions{$separator}->($_);
	    if (exists$XTars{$arr[$scancol]})
	    {
		#print "$arr[$scancol] score is poorer than decoy\n";
		next;
	    }
	    my $j=0;
	    if ($arr[$rankcol]==1)
	    {
		foreach my $val(@arr)
		{
		    $tar2d[$i][$j]=$val;
		    $j++;
		}
		$i++;
	    }
	}
	elsif($f==0)
	{
	    $tar2dheader=$_.",q-value";
	    $tar2dheader=~s/\t/,/g;
	    $f=1;
	    next;
	}
    }
    close TAR;
    @tar2d=sort{$b->[$scorecol]<=>$a->[$scorecol]} @tar2d;
    #print "Size of Tar 2d=",scalar @tar2d,"\n";
    my @ROC=@$ROCref; #2darray with ROC data
    
    my $ROC_header=shift@ROC;
    open (FDR,">$FDRoutfile") or die $!;
    print FDR $tar2dheader,"\n";
    my $start=0;
    #my $flag=0;#flag for q-values below FDR thr

    for (my $r=0;$r<@ROC;$r++)
    {
	for(my $i=$start;$i<@tar2d;$i++)
	{
	    if ($tar2d[$i][$scorecol]>$ROC[$r][0])
	    {
		#print"\$tar2d[$i][-1] is assigned a value of $ROC[$r][1]\n";<>;
		$tar2d[$i][scalar(@{$tar2d[$i]})]=$ROC[$r][1];
		
		print FDR shift@{$tar2d[$i]};#first element
		foreach(@{$tar2d[$i]})#from next element, put a comma first
		{
		    print FDR ",$_";
		    #print "$_,";
		}
		print FDR "\n";
		$start++;

		#print "\n";<>;
	    }
	    else
	    {
		#print"last for $tar2d[$i][$scorecol] .... $ROC[$r][0]";<>;
		$start=$i;
		last;
	    }		    
	}
    }
    
    print FDR "\nROC\n",join",",@$ROC_header,"\n";
    foreach (@ROC)
    {
	print FDR join ",",@{$_},"\n";
    }
    close FDR;
    return $FDRoutfile;   
}

################################################################################

sub filterConFDRPSMsAsc{
    my ($tar,$ROCref,$scorecol,$rankcol,$FDRoutfile,$XTar,$scancol,$separator)=@_;
    my %XTars=%$XTar;
    open TAR,$tar or die;
    my $f=0; #flag for header
    my $i=0; #result counter
    my @tar2d;
    my $tar2dheader;
    while(<TAR>)
    {
	chomp $_;
	if($f==1)
	{
	    my @arr=$actions{$separator}->($_);
	    if (exists$XTars{$arr[$scancol]})
	    {
		#print "$arr[$scancol] score is poorer than decoy\n";
		next;
	    }
	    my $j=0;
	    if ($arr[$rankcol]==1)
	    {
		foreach my $val(@arr)
		{
		    $tar2d[$i][$j]=$val;
		    $j++;
		}
		$i++;
	    }
	}
	elsif($f==0)
	{
	    $tar2dheader=$_.",q-value";
	    $tar2dheader=~s/\t/,/g;
	    $f=1;
	    next;
	}
    }
    close TAR;
    @tar2d=sort{$a->[$scorecol]<=>$b->[$scorecol]} @tar2d;
    #print "Size of Tar 2d=",scalar @tar2d,"\n";
    my @ROC=@$ROCref; #2darray with ROC data
    
    my $ROC_header=shift@ROC;
    open (FDR,">$FDRoutfile") or die $!;
    print FDR $tar2dheader,"\n";
    my $start=0;
    #my $flag=0;#flag for q-values below FDR thr

    for (my $r=0;$r<@ROC;$r++)
    {
	for(my $i=$start;$i<@tar2d;$i++)
	{
	    if ($tar2d[$i][$scorecol]<$ROC[$r][0])
	    {
		#print"\$tar2d[$i][-1] is assigned a value of $ROC[$r][1]\n";<>;
		$tar2d[$i][scalar(@{$tar2d[$i]})]=$ROC[$r][1];
		
		print FDR shift@{$tar2d[$i]};#first element
		foreach(@{$tar2d[$i]})#from next element, put a comma first
		{
		    print FDR ",$_";
		    #print "$_,";
		}
		print FDR "\n";
		$start++;
		#print "\n";<>;
	    }
	    
	    else
	    {
		#print"last for $tar2d[$i][$scorecol] .... $ROC[$r][0]";<>;
		$start=$i;
		last;
	    }		    
	}
    }
    
    print FDR "\nROC\n",join",",@$ROC_header,"\n";
    foreach (@ROC)
    {
	print FDR join ",",@{$_},"\n";
    }
    close FDR;
    return $FDRoutfile;  
}


1;