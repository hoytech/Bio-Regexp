use strict;

use Test::More;
use Data::Dumper;

use Bio::Regexp;


verify(Bio::Regexp->new->add('AA')->single_stranded,
       'AAA',
       'basic overlap',
       matches => [[0,2],[1,3]]);


verify(Bio::Regexp->new->add('GATA')->single_stranded,
       'GATTGATC',
       'basic no match',
       matches => []);


verify(Bio::Regexp->new->add('GAT{2,3}C[AT]')->single_stranded,
       'GGGATTCAAGATTTCTA',
       'simple regexp operations',
       matches => [[2,8],[9,16]]);


verify(Bio::Regexp->new->add('AA')->add('AC')->single_stranded,
       'AAC',
       'basic multi-regexp',
       matches => [[0,2],[1,3]]);


verify(Bio::Regexp->new->add('ATG'),
       'AAAGACATCC',
       'basic reverse complement',
       matches => [[5,8]]);


verify(Bio::Regexp->new->add('ATG'),
       'AAAGACATCC',
       'basic reverse complement',
       matches => [[5,8]]);


verify(Bio::Regexp->new->rna->add('AUG'),
       'GGCCGGCATAA',
       'RNA pattern, DNA string',
       matches => [[6,9]]);


verify(Bio::Regexp->new->add('TAT'),
       'AUGUAUAA',
       'DNA pattern, RNA string',
       matches => [[3,6],[4,7]]);


verify(Bio::Regexp->new->add('GAATTC'),
       'AGACTGAGAATTCGGG',
       'palindrome matches twice same place',
       matches => [[7,13],[7,13]]);


done_testing();


sub verify {
  my ($obj, $input, $desc, %checks) = @_;

  my @matches = $obj->match($input);

  print Dumper(\@matches) if $checks{dumper};

  if (exists $checks{matches}) {
    is(scalar @matches, scalar @{ $checks{matches}}, "$desc: length check");

    foreach my $i (0 .. $#matches) {
      is($matches[$i]->{start}, $checks{matches}->[$i]->[0], "$desc: $i start");
      is($matches[$i]->{end}, $checks{matches}->[$i]->[1], "$desc: $i end");
    }
  }
}
