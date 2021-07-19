use strict;
my %uniref_string;
my %uniref_map;
open(IN,$ARGV[0]);
while(<IN>)
{
	my $FF = $_;
	chomp $FF;
	my @a = split(/\t/,$FF);
	my @b = split(/\s/,$a[1]);
	my $i;
	for($i=0;$i<=$#b;$i++)
	{
		my $cazy = $b[$i];
		if(not exists $uniref_string{$cazy})
		{
			$uniref_string{$cazy} = "$a[0]";
			$uniref_map{$cazy}{$a[0]} = 1;
		}
		else
		{
			if(not exists $uniref_map{$cazy}{$a[0]})
			{
				$uniref_string{$cazy} .= "\t$a[0]";
			}
		}
	}
}

foreach my $ckeys (keys %uniref_string)
{
	print "$ckeys\t$uniref_string{$ckeys}\n";
}
			
