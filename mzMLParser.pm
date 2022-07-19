use strict;
use warnings;
use XML::Twig;
use MIME::Base64;
use IQ;
my  $prev_id="NA";
sub mzML_quant{
    my($spectra_folder,$par,$head,$plexref,$are)=@_;
    my %params=%$par;
    my %header=%$head;
    my %plex=%$plexref;
    my %area_name=%$are;
    ###################################################################
    opendir(DIR,$spectra_folder) || die "\nCannot open directory \"$spectra_folder\"\n";
    my @raw_files = readdir(DIR);
    closedir DIR;
    @raw_files=grep/\.mzML$/,@raw_files;
    print "Number of files to process:: ",scalar @raw_files,"\n";
    print "Type of area method:: \"$area_name{$params{Spectrum_area}}\"\n";
    ########################################################################
    foreach my $raw_input (@raw_files)
    {
        chomp $raw_input;
        print "Processing file :\t$raw_input\n";
        my $out=$spectra_folder.'/'.$raw_input.'_'.$params{Spectrum_area}.'_'.$params{Tol_min}.'.tsv';
        $out=~s/\.mzML//;
	my $out1=$out;
	$out1=~s/\.tsv/_noquant\.tsv/;
        open(my $fh,">$out") or die $!;##without my the file handle $fh will fall into strict refs clause and will not open file    #code will die
        print "Output file opened as: $out\n";
        print $fh "Title\tRetention_time\tBase_Intensity\tPEPMASS\tCharge\tMSLevel\tPrevious Scan";
	open(my $fh1,">$out1") or die $!;##without my the file handle $fh will fall into strict refs clause and will not open file    #code will die
        print "Output file no-Quant opened as: $out1\n";
        print $fh1 "Title\tRetention_time\tBase_Intensity\tPEPMASS\tCharge\tMSLevel\tPrevious Scan";
        foreach my $tag(@{$header{$params{Label_type}}{$params{Plex}}})
        {
            #print $fh "\t$tag-Area\t$tag-Error\t$tag-Peak_Max\t$tag-Ratio\t$tag-Ratio_error";
            print $fh "\t$tag/$params{Base_reporter}";
        }
        foreach my $tag(@{$header{$params{Label_type}}{$params{Plex}}})
        {
            print $fh "\t$tag/$params{Base_reporter}-Err";
        }
        foreach my $tag(@{$header{$params{Label_type}}{$params{Plex}}})
        {
            print $fh "\t$tag-Peak_Max";
        }
        foreach my $tag(@{$header{$params{Label_type}}{$params{Plex}}})
        {
            print $fh "\t$tag-Area";
        }
        foreach my $tag(@{$header{$params{Label_type}}{$params{Plex}}})
        {
            print $fh "\t$tag-Error";
        }
        print $fh "\n";
	print $fh1 "\n";
        #####################################################################################################
        #print "Before twig...\n";<>;
        #print "\nPARAM KEYS==",keys %params,"\n";<>;
        #print "\nPLEX KEYS==",keys %plex,"\n";<>;
	
    ##########################################################################TODOTODO\
    #for ms3 funday
    #read file in while and concat
    #MS1, MS2, Ms3
    #then give to twig one by one and qunatitate
    #send MS2 level info side by side
    #once ms1 recomes
    #flush all hashes
    #
    #########################################################################
    #open(MZML,"$raw_input") or die $!;
    #my %mzmlspec;
    #my $flag=0;
    #while (<MZML>)
    #{
    #    chomp $_;
    #    
    #}
    #
    
    
	#http://www.perlmonks.org/?node_id=901332
        my $t= XML::Twig->new(twig_roots=>{'spectrum'=> sub{Parse_mzML_spec(@_,\%plex,\%params,$fh,$fh1)}}); # give the root tag of the file 
        #print Dumper($t);<>;
        $t->parsefile($spectra_folder."/".$raw_input) or print("File not found\n");#check TODO
        #print "\nENDED\n\n";<>;
        $t->purge;
    }
}
sub Parse_mzML_spec{
    my ($twig,$ele,$plexref,$par,$fh,$fh1)=@_;
    my %plex=%$plexref;
    my %params=%$par;
    my $r=$ele;
    my $flag=0;
    my %spectra;##MS2=[array]
    ##my %scan;
    my $id;
    #my %res;
    #my $label=$plex{$params{Label}};
    my($charge,$prec_intent,$activationMethod)=('NA','NA','NA');
    my @scan_arr;#pushing all hashes to one array
    my ($mz_arr,$int_arr)=(0,0);
    $id=$r->{'att'}->{'id'};
    #$id=~m/(.+) scan=(.+)$/;
    #$id=$2;
    $spectra{TITLE}=$id;
    #print "id\t$id\n";<>;
    $spectra{PrevScan}="NA";  
    while($r=$r->next_elt())
    {
        #if ($r->gi eq 'spectrum')
        #{
        #    $id=$r->{'att'}->{'id'};
        #    $id=~m/(.+) scan=(.+)$/;
        #    $id=$2;
        #    $spectra{TITLE}=$id;
        #    #code
        #}
        #
        #
        
        if($r->gi eq 'precursor')# only present if mslevel2/3
	{
            my  $prev_scan=$r->{'att'}->{'spectrumRef'};
	   #print "$r\t$prev_scan";<>;
       
        #$prev_scan=~m/(.+) scan=(.+)$/;
	    #$prev_scan=$2;
        if (!defined $prev_scan || $prev_scan eq "")
        {
           $prev_scan=$prev_id;
        }
        
	    $spectra{PrevScan}=$prev_scan;
	    #$res{ $id}{Prev_scan}=$prev_scan;
	    #print "id\t$id\tprev_scan\t$prev_scan\n";<>;
            #print OUT "\n$id\t";
        }
        
        my $acc=$r->{'att'}->{'accession'} || 0 ;#gets accesions in all cvparams
	#print $acc,"\n";<>;
        if ($acc eq "MS:1000511")#mslevel
        {
            my $mslevel=$r->{'att'}->{'value'};
	    $spectra{MSLevel}=$mslevel;
            #$res{$id}{mslevel}=$mslevel;
            #print OUT "$mslevel\t";
        }
	elsif($acc eq "MS:1000041")#charge
        {
            my $z=$r->{'att'}->{'value'};
	    $spectra{Z}=$z;
	    #print OUT "$charge\t";    
        }
        elsif($acc eq "MS:1000016")#rt
        {
            my $rtinsec=$r->{'att'}->{'value'};
	    $spectra{RT}=$rtinsec;
	    #$res{$id}{RT}=$rtinsec;
            #print OUT "$rtinsec\t";
        }
        elsif($acc eq "MS:1000574")#zlib
        {
            my $zlib="zlib";
	    $spectra{ZLIB}=$zlib;
            #print "COMP=$zlib\t";    
        }
        elsif($acc eq "MS:1000744")#mz
        {
            my $mz=$r->{'att'}->{'value'};
            $spectra{MZ}=$mz;
	    #print OUT "$mz\t";
	    #$res{$id}{MZ}=$mz;
        }    
        elsif($acc eq "MS:1000042")#int
        {
            my $Int=$r->{'att'}->{'value'};
            $spectra{INT}=$Int;
	    #print OUT "$Int\t";
	    #$res{$id}{Int}=$Int;
        }
        elsif($acc eq "MS:1000576")#no compression
        {
            my $compression="none";#TODO
            $spectra{compression}=$compression;
	    #print "COMP=$compression\t";    
        }
        elsif($acc eq "MS:1000514")#m/z array
        {
            $flag=1;
            next;
            #if ($r->gi eq "binaryDataArray")
            #print "m/zarr=present\t";
	}
        elsif ($r->gi eq "binary" && $flag==1)
        {
            my $val=$r->text;
	    $spectra{MZ_frag}=$val;
	    #$res{$id}{MZ_bin}=$val;
	    #my $val=defined($r->first_child_text('binary'))?($r->first_child_text('binary')):'NA';#peaks
            #print "\n::m/z:: $val\t";<>;
        }
        elsif($acc eq "MS:1000515")#int array
        {
            $flag=2;
            next;
	}
        #if ($r->gi eq "binaryDataArray")
        #{
        #    my $val=defined($r->first_child_text('binary'))?($r->first_child_text('binary')):'NA';#peaks
        #    print "\n::INT:: $val\t";
        #
        #}
        #code
        #print "int_arr=present\t";
         elsif ($r->gi eq "binary"&& $flag==2)    
        {
            my $val=$r->text;
	    $spectra{Int_frag}=$val;
	    #$res{$id}{INT_bin}=$val;
            #my $val=defined($r->first_child_text('binary'))?($r->first_child_text('binary')):'NA';#peaks
            #print "\n::int:: $val\t";<>;
        }    
    }
    if ($spectra{MSLevel}>1)
    {	
        #print "MSLEVEL=>$spectra{TITLE}\n";<>;
	if (exists($spectra{compression}))
	{
            #print "Compression=>$spectra{TITLE}\n";<>;
	    my $Base64DecodedMZ=decode_base64($spectra{MZ_frag});
	    my $Base64DecodedIntensity=decode_base64($spectra{Int_frag});
	    # unpacking directly unzipped base64 decoded data 
            # Refernce: http://www.tutorialspoint.com/perl/perl_unpack.htm
	    # d - A double-precision floating-point number
	    my @DataMZ=unpack("d*",$Base64DecodedMZ);
	    # f - A single-precision floating-point number
	    my @DataIntensity=unpack("f*",$Base64DecodedIntensity);
	    for(my $i=0;$i<@DataMZ;$i++)
	    {
#            {
#		    if ($params{Tol_unit} eq "Da")
#		    {
#                        #print "$DataMZ[$i]";<>;
#			if ($DataMZ[$i]>=($plex{$params{Label_type}}{$params{Plex}}[0] - $params{Tol_min}))
#			{
#                            #print "$DataMZ[0]";<>;
#                            #print "\nhere min\n";
#			    if($DataMZ[$i]<=($plex{$params{Label_type}}{$params{Plex}}[-1] + $params{Tol_max}))
#			    {
#                                #print "$DataMZ[0]";<>;
#                                #print "\nhere max\n";
#			        if ($params{Label_type} eq "iTRAQ" && $DataMZ[$i]=~m/^120/)
#			        {
#			            #code
#			        }
#			        else
#			        {
#                                    #print "$DataMZ[0]";<>;
#                                    #print "\nhere string\n";
			            my $string="$DataMZ[$i] $DataIntensity[$i]";
			            push(@{$spectra{MS2}},$string);#keys->MS2 value->[(m/z int),(m/z int)....]
			#	}
			#    }
			#}
		    }
#		    else#ppm
#		    {
#			if ($DataMZ[$i]>=($plex{$params{Label_type}}{$params{Plex}}[0] - (($plex{$params{Label_type}}{$params{Plex}}[0]*$params{Tol_min})/1000000)))
#			{
#			    if($DataMZ[$i]<=($plex{$params{Label_type}}{$params{Plex}}[-1] + (($plex{$params{Label_type}}{$params{Plex}}[-1]*$params{Tol_max})/1000000)))
#			    {
#			        if ($params{Label_type} eq "iTRAQ" && $DataMZ[$i]=~m/^120/)
#			        {
#			            #code
#			        }
#			        else
#			        {
#			            my $string="$DataMZ[$i] $DataIntensity[$i]";
#			            push(@{$spectra{MS2}},$string);#keys->MS2 value->[(m/z int),(m/z int)....]
#				}
#			    }
#			}
#		    }
#		}
#		if (!exists($spectra{Z}))
#                {
#                    $spectra{Z}='NA';
#                }
#                if (!exists($spectra{RT}))
#                {
#                    $spectra{RT}='NA';
#                }
#                if (!exists($spectra{INT}))
#                {
#                    $spectra{INT}='NA';
#                }
                if (exists $spectra{MS2}) #if no peaks are available the subroutine will not be called
                {
		    #print "$spectra{TITLE}\n";<>;
                    my($spectra_hashref,$area_hashref,$rep_count)= spectrum_quantitator(\%spectra,\%params,\%plex);
                    if ($rep_count>=$params{min_reporters})
                    {
			print_spec($spectra_hashref,$area_hashref,$fh);
                    }                     
                        #}
                }
    
                else
                {
		    print $fh1 "$spectra{TITLE}\t$spectra{RT}\t$spectra{INT}\t$spectra{MZ}\t$spectra{Z}\n";
                }
                #undef %spectra;##TODO remove one statement after testing
                %spectra=();
		#$spectra{MS2}=@mz_arr;
	}
	
	else#zlib TODO
	{
	    print "zlib=>$spectra{TITLE}\n";
	    print "zlib compression not supported\n";
	    exit;
	    #<>;
	}
	
    }
    $prev_id=$id;
}
1;