use strict;
use warnings;
use constant H=>1.007825;
use Data::Dump qw(dump);
########################################################################
#Author: Amit Kumar Yadav
#Suruchi Aggarwal
#Puneet KUmar Kadimi

#amit.yadav@thsti.res.in

############################################################################
# PROGRAM TO READ MS SPECTRA FILES 
#########################################################################

sub spectra_parser{
	my ($RelInt,$MinPeaks,$spectra_file)=@_;
	#print $spectra_file,"\n";
	$spectra_file=~m/\.(...)$/; #get extension
	my $ext=lc($1);
	my @ExpSpectra;#2D array structure [0]=ScanID/SpectrumID; [1]=precursor_mz; [2]=intensity; [3]=charge [4]=RT; [5]=spectrum hashref; [6]=Mr experimental;
	if ($ext eq "mgf")
	{
		@ExpSpectra=MGFReader($RelInt,$MinPeaks,$spectra_file);
	}
	elsif($ext eq "pkl")
	{
		@ExpSpectra=parse_pkl($RelInt,$MinPeaks,$spectra_file);
	}
	elsif($ext eq "dta")
	{
		@ExpSpectra=parse_dta($RelInt,$MinPeaks,$spectra_file);
	}
	else
	{
		die "Unsupported spectra file format\nOnly mgf,PKL and DTA are allowed\n";
	}
	return @ExpSpectra;
}

###########################################################################################
# MGF Parser, to allow skipping spectra without charges
###########################################################################################
sub MGFReader{
	my ($RelInt,$MinPeaks,$spectra_file)=@_;
	my $flag=0;
	my $i=0; #for diff spectra scans
	my $j=0; # spectrum 2-D for each scan
	my @ExpSpec; #2D array of scans
	my @spectrum; #mz-int pairs 2d
	my $flag_last=0;
	open MGF,$spectra_file or die "$spectra_file not found\n";
	#open (SELECTED,">".$spectra_file."_SELECTED") or die("cannot create output file");
	my @mgfdata=<MGF>;
	close MGF;
	my $con=join"", @mgfdata;
	@mgfdata=split ('BEGIN IONS',$con);
	shift @mgfdata;#remove 1st empty element
	
	print scalar @mgfdata," spectra present(total) in MGF file\n";
	my $k=0;
	for (my$i =0;$i<@mgfdata;$i++)
	{
		my @arr=GetSpectrum_MGF($RelInt,$MinPeaks,$mgfdata[$i]);
		if ($arr[0]eq 'NA') #if spectrum is skipped
		{
			next;
		}
		else
		{
			for ($j=0;$j<@arr;$j++)
			{
				$ExpSpec[$k][$j]=$arr[$j];
			}
			#print "printing final spectrum...\n";
			#print($ExpSpec[$k][5]);<>;
			
			$k++;
		}
	}
	print scalar @ExpSpec," spectra returned from MGF file after filtering\n";

	return @ExpSpec;
}

sub GetSpectrum_MGF{
	my($RelInt,$MinPeaks,$sp)=@_;
	my @arr=split"\n",$sp; #Every line of the spectrum
	my $endtag=pop @arr; #remove "END IONS" tag
	#print "\$endtag=$endtag\n";
	my @spectrum;#Array for ALL Scan info
	my @Spec2D; #2d array for holding fragment mz->intensity pairs
	my %attrib;
	my $j=0;
	foreach(@arr)
	{
		chomp $_;
		if ($_=~m/^\s*$/)
		{
			#print "\$_=\"$_\"\n";
			next;
		}		
		if ($_=~m/^\d+/)
		{
			my @sp=split /\s+/,$_;
			$Spec2D[$j][0]=$sp[0];
			$Spec2D[$j][1]=$sp[1];
			$j++;
		}
		else #Save headers
		{
			#print $_;
			my @temp=split '=',$_;
			$attrib{$temp[0]}=join"=",@temp[1..$#temp];#if title has = sign
		}		
	}
	foreach(keys %attrib)
	{
		#print $_.'='.$attrib{$_}."\n";<>;
		if ($_ eq 'TITLE')
		{
			$spectrum[0]=$attrib{$_};
			$spectrum[0]=~s/,/;/g;
		}
		elsif ($_ eq 'PEPMASS')
		{
			#print $_ ,"->$attrib{$_} matches PEPMASS\n";
			my @arr2=split/\s+/,$attrib{$_};
			$spectrum[1]=$arr2[0]; #mz
			$spectrum[2]=$arr2[1]; #intensity
			#print "Precursor in READ loop 2=$spectrum[1]  ...\n";
			#print "Intensity in READ loop 2=$spectrum[2]  ...\n";<>;
			
		}
		elsif ($_ eq 'CHARGE')
		{
			$spectrum[3]=$attrib{$_};
			#print "charge in calc loop=$ExpSpec[$i][3]  ...\t";
			#print "Precursor in calc loop=$ExpSpec[$i][1]  ...\t";
			$spectrum[3]=~s/\+//g;
			$spectrum[3]=~s/\-//g;
		}
		elsif ($_ eq 'RTINSECONDS')
		{
			$spectrum[4]=$attrib{$_};
		}
		else
		{
			#print"$_ was missed \n";
		}
	}
	#Check if charge is defined, else throw the spectrum away
	if (defined($spectrum[3]))
	{
		#print "MZ=$spectrum[1]\tZ=$spectrum[3]\n";
		$spectrum[6]=($spectrum[1]*$spectrum[3])-($spectrum[3]*H); #Exp Mr (Deconvolution step)
	}
	else #delete all data for current scan and begin again
	{
		return 'NA'; #Code to throw away this spectrum
	}
	
	#Check peak count critria and filter spectrum
	@Spec2D=intensity_filter($RelInt,$MinPeaks,\@Spec2D,$spectrum[6]);
	if ((scalar @Spec2D)>=$MinPeaks)
	{
		my $spec_string='';
		foreach my $value (@Spec2D)
		{
			$spec_string.="$value->[0] $value->[1]\n";
		}
		$spectrum[5]=$spec_string;
	}
	else
	{
		return 'NA';
	}
	#END OF NEW CODE Saturday, April 17, 2010
	return @spectrum;
}

###########################################################################################
# MGF Parser, OLDER
#########################################################################################
sub parse_mgf{
	my ($RelInt,$MinPeaks,$spectra_file)=@_;
	my $flag=0;
	my $i=0; #for diff spectra scans
	my $j=0; # spectrum 2-D for each scan
	my @ExpSpec;
	my @spectrum;
	my $flag_last=0;
	open MGF,$spectra_file or die "$spectra_file not found\n";
	#open (SELECTED,">".$spectra_file."_SELECTED") or die("cannot create output file");
	while(<MGF>)
	{
		chomp $_;
		if ($_=~m/^\s+$/) 
		{
			next;
		}
		if ($_=~m/^BEGIN IONS/)
		{
			$flag=1;
			next;
		}
		elsif($_=~m/^END IONS$/)
		{
			$flag=0;
			my @temp=@spectrum;
			#5 AUG 2013 AKY
			# This chunk checks if charge is defined and then deconvolutes
			if (defined($ExpSpec[$i][3]))
			{
				$ExpSpec[$i][6]=($ExpSpec[$i][1]*$ExpSpec[$i][3])-($ExpSpec[$i][3]*H); #Exp Mr (Deconvolution step)
			}
			else #delete all data for current scan and begin again
			{
				
				#ADD Charge estimation here
				#if ($ExpSpec[$i][1]>=($temp[-1][0]+5))#within 5 Da
				#{
				#	$ExpSpec[$i][3]=1;
				#}
				#elsif ( ($ExpSpec[$i][1]<$temp[-1][0]-5) &&(2*$ExpSpec[$i][1]>=($temp[-1][0]+5))) #B/W 1x and 2x of mz
				#{
				#	$ExpSpec[$i][3]=2;
				#}
				#elsif ( ($ExpSpec[$i][1]<$temp[-1][0]-5) &&(3*$ExpSpec[$i][1]>=($temp[-1][0]+5))) #B/W 1x and 2x of mz
				#{
				#	$ExpSpec[$i][3]=3;
				#}
				#elsif ( ($ExpSpec[$i][1]<$temp[-1][0]-5) &&(4*$ExpSpec[$i][1]>=($temp[-1][0]+5))) #B/W 1x and 2x of mz
				#{
				#	$ExpSpec[$i][3]=4;
				#}
				#else {$ExpSpec[$i][3]=2;}
				
			}
			
			#NEW CODE Saturday, April 17, 2010
			#print "Precursor in main loop=$ExpSpec[$i][0] with ", scalar @temp," peaks\n";
			@temp=intensity_filter($RelInt,$MinPeaks,\@temp,$ExpSpec[$i][6]);
			if ((scalar@temp)>=$MinPeaks) 
			{
				my $spec_string='';
				foreach my $value (@temp)
				{
					$spec_string.="$value->[0] $value->[1]\n";
				}
				#	print $ExpSpec[$i][0]."\n";
				#print SELECTED "BEGIN IONS\nTITLE=$ExpSpec[$i][0],\nPEPMASS=$ExpSpec[$i][1]\nCHARGE=$ExpSpec[$i][3]".'+'."\n";
				#	print SELECTED $spec_string;
				#	print SELECTED "END IONS\n";
				$ExpSpec[$i][5]=$spec_string;
				@temp=();
				$spec_string='';
				@spectrum=();
				$i++;
				$flag_last=0;
				$j=0;
				next;
			}
			else
			{
				$j=0;
				@spectrum=();
				$flag_last=1;
				next;
			}
			#END OF NEW CODE Saturday, April 17, 2010
		}
		######################
		elsif ($flag==1)#add spectral info
			{
				if ($_=~m/^\d+/) 
				{
					my @sp=split /\s+/,$_;
					$spectrum[$j][0]=$sp[0];
					$spectrum[$j][1]=$sp[1];
					$j++;
				}
				elsif ($_=~m/^TITLE/) 
				{
					my @arr=split "=",$_;
					shift @arr;
					my $title=join "=",@arr; #scan id with = sign dont change now (15 April 2010)
					$title=~s/,/;/g;
					$ExpSpec[$i][0]=$title;
				}
				elsif ($_=~m/^PEPMASS/) 
				{
					my @arr=split /=/,$_;
					if ($arr[1]=~m/\s+/g) 
						{
						my @arr2=split/\s+/,$arr[1];
						$ExpSpec[$i][1]=$arr2[0]; #mz
						#print "Precursor in READ loop=$ExpSpec[$i][1]...\t";
						$ExpSpec[$i][2]=$arr2[1]; #intensity
						}
					else
						{
						$ExpSpec[$i][1]=$arr[1]; #mz
						$ExpSpec[$i][2]=undef; #intensity
						}
					#print "Precursor in READ loop 2=$ExpSpec[$i][1]  ...\n";

				}
				elsif ($_=~m/^CHARGE/) 
				{
					my @arr=split /=/,$_;
					$ExpSpec[$i][3]=$arr[1];
					#print "charge in calc loop=$ExpSpec[$i][3]  ...\t";
					#print "Precursor in calc loop=$ExpSpec[$i][1]  ...\t";
					$ExpSpec[$i][3]=~s/\+//g;
					$ExpSpec[$i][3]=~s/\-//g;
					
				}
				elsif ($_=~m/^RTINSECONDS/) 
				{
					my @arr=split /=/,$_;
					$ExpSpec[$i][4]=$arr[1];
				}
			}
		elsif($flag==0)
		{
			next;
		}

	}
	if ($flag_last==1)
	{
		pop(@ExpSpec);
	}
	close MGF;
	#close SELECTED;
	#dump(@ExpSpec);
	return@ExpSpec;
}
###########################################################################################
###########################################################################################

sub parse_dta {
	my ($RelInt,$MinPeaks,$spectra_file)=@_;
	my $flag=0;
	my $i=0;#Counter for spectrum/scan ID
	my $j=0; #counter for spectrum 2D
	my $lc=0;#line counter to catch errors in input
	my @ExpSpec;
	my @spectrum;
	open DTA, $spectra_file or die $!;
	while(<DTA>)
	{
		chomp $_;
			$lc++;
			if (($_!~m/^\s+$/) && ($flag==0))
				{
					my @row=split /\s|\t/,$_;
					$ExpSpec[0][0]=$i+1;#assign ID =1 to precursor
					$ExpSpec[0][1]=$row[0]; # precursor mz [MH+]
					$ExpSpec[0][2]=undef; # intensity not in DTA file
					$ExpSpec[0][3]=$row[1]; #charge
					$ExpSpec[0][4]=undef; #no RT in DTA file
					$ExpSpec[0][6]=($ExpSpec[$i][1]*$ExpSpec[$i][3])-(H*$ExpSpec[$i][3]); #Exp Mr (Deconvolution step)
					$flag=1;
				}
			elsif ($flag==1)
				{
					my @row=split /\s|\t/,$_;
					$spectrum[$j][0]=$row[0];
					$spectrum[$j][1]=$row[1];
					$j++;
				}
			else
				{
					die "Error in DTA format at line $lc\n";
				}
		}
	close DTA;
	@spectrum=intensity_filter($RelInt,$MinPeaks,\@spectrum,$ExpSpec[0][6]);
	if ((scalar @spectrum)>=$MinPeaks) 
		{
					my $spec_string='';
					foreach my $value (@spectrum)
						{
							$spec_string.="$value->[0] $value->[1]\n";
						}
					$ExpSpec[$i][5]=$spec_string;
					@spectrum=();
					$spec_string='';
			#$ExpSpec[0][5]=\@spectrum;
					return @ExpSpec;
		}
	else
		{
			my @nullarray=();
			return @nullarray;
		}
}

###########################################################################################
###########################################################################################

sub parse_pkl {
	my ($RelInt,$MinPeaks,$spectra_file)=@_;
	my $spflag=0;#spflag=1 if spectrum header found(keep filling mass int pairs); if spflag=0, find new spectrum header;
	my $flag=0;#flag to detect 1st spectrum
	my $i=0;#Counter for spectrum/scan ID
	my $j=0; #Counter for spectrum 2D
	my $lc=0;#line counter to catch errors in input
	my $mz=0;#flag for frag mass intensity pairs
	my @ExpSpec;
	my @spectrum;
	open PKL,$spectra_file or die $!;
	while(<PKL>)
		{
			chomp $_;
			$lc++;
			#print"\n\nline no.\t$lc\n";
			if ($lc==1) # capture 1st spectrum
				{
					my @row=split /\s|\t/,$_;
					if (scalar@row==3) 
						{
							$ExpSpec[$i][0]=$i; #scan ID
							$ExpSpec[$i][1]=$row[0]; #Precursor mz
							$ExpSpec[$i][2]=$row[1]; # precursor Intensity
							$ExpSpec[$i][3]=$row[2]; # charge
							$ExpSpec[$i][4]=undef; #RT
							$ExpSpec[$i][6]=($ExpSpec[$i][1]*$ExpSpec[$i][3])-(H*$ExpSpec[$i][3]); #Exp Mr (Deconvolution step)
							$spflag=1; #to record mz-int pairs
							next;
						}
					else
						{
							die"error in PKL file format at line $lc\n";
						}
				}
			#################################################
			if (($_ eq '') or (eof PKL))
				{
					if (scalar(@spectrum)==0)
						{
							$spflag=0;#start looking for new spectrum
							next;
						}
					else
						{
							my @temp=@spectrum;
							@temp=intensity_filter($RelInt,$MinPeaks,\@temp,$ExpSpec[$i][6]);
							if (scalar@temp>=$MinPeaks) 
								{
									my $spec_string='';
									foreach my $value (@temp)
										{
											$spec_string.="$value->[0] $value->[1]\n";
										}
									$ExpSpec[$i][5]=$spec_string;
									@temp=();
									$spec_string='';
									@spectrum=();
									$i++;
									$j=0;
									$spflag=0;#start looking for new spectrum
								}
							else
								{
									@spectrum=();
									$j=0;
									$spflag=0;
									if (eof(PKL)) 
										{
											pop @ExpSpec;
										}
								}
						}
				}
			else # Spectral assignment of mass-int pairs
				{
					my @row=split /\s|\t/,$_;
					if ( ( $spflag==0) && ( scalar @row==3) ) #new spectrum header
						{
							$ExpSpec[$i][0]=$i; #scan ID
							$ExpSpec[$i][1]=$row[0]; #Precursor mz
							$ExpSpec[$i][2]=$row[1]; # precursor Intensity
							$ExpSpec[$i][3]=$row[2]; # charge
							$ExpSpec[$i][4]=undef; #RT
							$ExpSpec[$i][6]=($ExpSpec[$i][1]*$ExpSpec[$i][3])-($ExpSpec[$i][3]*H); #Exp Mr (Deconvolution step)
							$spflag=1; #start looking for mz-int pairs
							next;
						}
					elsif ( ( ($spflag==0) && (scalar @row!=3) ) )
						{
							die "Error in PKL format at line no. $lc\n";
						}
					elsif ( ( $spflag==1) && ( scalar @row==2) )
						{
							$spectrum[$j][0]=$row[0];
							$spectrum[$j][1]=$row[1];
							$j++;
						}
					elsif ( ($spflag==1) &&(scalar @row!=2) )
						{
							die "Error in finding mz - intensity pair in PKL file format at line $lc\n";
						}
					else
						{
							die "Some uncaught error in PKL file format at line $lc\n";
						}
				}
		}
	return @ExpSpec;
}

######################################################
# Added Monday, May 31, 2010
# Removes base peak ± 5 Da and peak selection from zones
######################################################
sub intensity_filter 
	{
		my ($RelInt,$minpeaks,$ref,$precursor_mass)=@_;
		my %sp;
		foreach (@{$ref}) 
			{
				$sp{$_->[0]}=$_->[1];
			}
		my@intensity=sort{$a<=>$b}values%sp;
		my ($min,$max) = (sort {$a<=>$b}@intensity)[0,-1];
		#print"\$precursor_mass=$precursor_mass\n";
		foreach my $key (keys(%sp))
				{
					if($key>$precursor_mass-5)
					{
						delete $sp{$key};
					}
				}
	if (scalar keys(%sp)<=1)
	{
		my @spec;
		return @spec;
	}
	
	my $total_ions_required=int($precursor_mass/110)*4;
	#print "\$precursor_mass=$precursor_mass\t\$total_ions_required 1 =$total_ions_required\n";
	if(($total_ions_required%5)>0)
		{
			$total_ions_required=(int($total_ions_required/5)+1)*5;
		}
	#print $total_ions_required, " total ion reqd\n";<>;
	
	my $bin_number=$total_ions_required/5;
	#print "\$bin_number=$bin_number\n";
	#print "\$total_ions_required=$total_ions_required\t\$bin_number=$bin_number\n\n";
	
	my ($min_mass,$max_mass)=(sort{$a<=>$b} keys(%sp))[0,-1];
	my $increment=($max_mass-$min_mass)/$bin_number;
	my %sp_mod;
	for(my $bin_value=$min_mass;$bin_value<$max_mass;$bin_value+=$increment)
	{
		my @sub_array=();
		my $i=0;
		foreach my $key(sort{$a<=>$b} keys(%sp)) 
			{
				if(($key>=$bin_value)&&($key<$bin_value+$increment))
				{
					$sub_array[$i][0]=$key;
					
					$sub_array[$i][1]=$sp{$key};
					$i++;
				}
				else
				{
					next;
				}
			}
		my @sorted_sub_array=sort{$b->[1]<=>$a->[1]}@sub_array;
		my $max_limit=scalar(@sorted_sub_array);
		if($max_limit>5)
			{
				$max_limit=5;
			}
		for(my $j=0;$j<$max_limit;$j++)
			{
				$sp_mod{$sorted_sub_array[$j][0]}=$sorted_sub_array[$j][1];
			}
	}
	foreach my $in (keys%sp_mod) 
		{
			if (($sp_mod{$in}*100/$max)>= $RelInt)
				{
					next;
				}
			else{
					delete $sp_mod{$in};
				}
		}		
		
	if (scalar(keys %sp_mod)>=$minpeaks)
		{
			my @spec;
			my $i=0;
			foreach (sort{$a<=>$b}keys%sp_mod) {
				$spec[$i][0]= $_;
				$spec[$i][1]= $sp_mod{$_};
				$i++;
			}
			return @spec;
		}
	else
		{
			my @spec;
			return @spec;
		}
}
1;