use strict;

my %cazy_ec_mapping;

open(MAP,"Combined_EC_CAZy_Mappings.txt");
while(<MAP>)
{

	my $FF = $_;
	chomp $FF;
	my @a = split(/\t/,$FF);
	$cazy_ec_mapping{$a[0]} = $a[1];
	#print "$a[0]\n";
}

my %uniref_cazy_mapping;

open(FILE,"bzcat $ARGV[0] |");
while(<FILE>)
{
	my $FF = $_;
	chomp $FF;
	my @a = split(/\t/,$FF);
	#print "$a[1]\n";
	if(exists $cazy_ec_mapping{$a[1]})
	{
		my $i;
		for($i=2;$i<=$#a;$i++)
		{
			print "$a[$i]\t$cazy_ec_mapping{$a[1]}\n";
		}
	}
}
