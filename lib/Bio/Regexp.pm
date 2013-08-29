package Bio::Regexp;

our $VERSION = '0.100';

use v5.10;
use common::sense;

use Regexp::Exhaustive;

use Bio::Regexp::AST;


sub new {
  my ($class, @args) = @_;

  my $self = {};
  bless $self, $class;

  return $self;
};


sub add {
  my ($self, $regexp) = @_;

  die "Can't add new regexp because regexp has already been compiled" if $self->{compiled_regexp};

  push @{ $self->{regexps} }, $regexp;

  return $self;
}


sub dna { _arg($_[0], 'type', 'dna') }
sub rna { _arg($_[0], 'type', 'rna') }
sub protein { _arg($_[0], 'type', 'protein') }

sub circular { _arg($_[0], 'circular', 1) }
sub linear { _arg($_[0], 'circular', 0) }

sub single_stranded { _arg($_[0], 'strands', 1) }
sub double_stranded { _arg($_[0], 'strands', 2) }

sub strict_thymine_uracil { _arg($_[0], 'strict_thymine_uracil', 1) }
sub strict_case { _arg($_[0], 'strict_case', 1) }


sub _arg {
  my ($self, $arg, $val) = @_;
  die "Can't set $arg to $val because it was already set to $val" if exists $self->{$arg};
  die "Can't set $arg to $val because regexp has already been compiled" if $self->{compiled_regexp};
  $self->{arg}->{$arg} = $val;
  return $self;
}



sub _process_args {
  my ($self) = @_;

  $self->{type} //= 'dna';

  if ($self->{type} eq 'dna') {
    $self->{arg}->{strands} //= 2;
  } elsif ($self->{type} eq 'rna') {
    $self->{arg}->{strands} //= 1;
  } elsif ($self->{type} eq 'protein') {
    die "protein search not implemented";
  }
}


sub compile {
  my ($self) = @_;

  return if $self->{compiled_regexp};

  $self->_process_args;

  my $regexp_index = 0;
  my @regexp_fragments;

  foreach my $regexp (@{ $self->{regexps} }) {
    ## Parse

    $regexp =~ $Bio::Regexp::AST::parser || die "Couldn't parse regexp: $regexp";
    my $ast = \%/;

    ## Compute meta data

    my ($min, $max) = $ast->{regexp}->compute_min_max;

    $self->{min} = $min if !defined $self->{min} || $min < $self->{min};
    $self->{max} = $max if !defined $self->{max} || $max > $self->{max};

    ## Main "sense" strand

    my $rendered = $ast->{regexp}->render;

    push @regexp_fragments, "$rendered(?{ $regexp_index })";
    $regexp_index++;

    my $component = { regexp => $regexp, };

    $component->{strand} = 1 if $self->{arg}->{strands} == 2;

    push @{ $self->{components} }, $component;

    ## Reverse complement strand

    if ($self->{arg}->{strands} == 2) {
      $rendered = $ast->{regexp}->reverse_complement->render;

      push @regexp_fragments, "$rendered(?{ $regexp_index })";
      $regexp_index++;

      my $component = { regexp => $regexp, strand => 2, };

      push @{ $self->{components} }, $component;
    }
  }

  my $compiled_regexp = ($self->{arg}->{strict_case} ? '' : '(?i)') .
                        '(' .
                        join('|', @regexp_fragments) .
                        ')';

  {
    use re 'eval';
    $self->{compiled_regexp} = qr{$compiled_regexp};
  }

  return $self;
}




sub match {
  my ($self, $input, $callback) = @_;

  $self->compile;

  my @matches = Regexp::Exhaustive::exhaustive($input => $self->{compiled_regexp},
                                               qw[ $1 @- @+ $^R ]);

  my @output;

  foreach my $match (@matches) {
    my $element = {
                    match => $match->[0],
                    start => $match->[1]->[0],
                    end => $match->[2]->[0],
                    %{ $self->{components}->[$match->[3]] },
                  };

    push @output, $element;
  }

  return @output;
}




1;



__END__

=head1 NAME

Bio::Regexp - Exhaustive DNA/RNA/protein regexp searches

=head1 SYNOPSIS

    my @matches = Bio::Regexp->new
                             ->add('AGC[2]YY[GTZ]{3,5}GCGC')
                             ->dna
                             ->circular
                             ->match($input);

    ## Example match:
    ##   {
    ##     match => 'AGCCGGGTGCC',
    ##     start => 203,
    ##     end => 214,
    ##     strand => 1,
    ##   }


=head1 DESCRIPTION

This module is for searching inside DNA or RNA or protein sequences. The sequence to be found is specified by a regular expression. Actually the search language is a restricted sub-set of perl regular expressions. IUPAC short form expansions are also in effect for DNA/RNA (they are kind of like character classes).

One goal of this module is to provide a complete search. Given the particulars of a sequence (DNA/RNA/protein, linear molecule/circular plasmid, single/double stranded) it attempts to figure out all of the possible matches with no missing, false-positive, or duplicate matches.

It handles cases where matches overlap in the sequence or your regular expression can match in multiple ways. For circular DNA plasmids it will find matches even if they span the "end" followed by the "beginning" of the arbitrary location in the circular sequence selected to be interbase coordinate 0. For double-stranded DNA it will find matches on the reverse complement strand as well.

The typical use case of this module is to search for few or many small patterns in large amounts of input data. Although it is optimised for that workflow it is efficient at a variety of tasks. Usually none of the input sequence data is copied at all except to extract matches (though this can be disabled too leaving only indices).



=head1 INPUT FORMAT

The input string passed to C<match> must be a nucleotide sequence for now (protein sequences will be supported better soon). There must be no line breaks or other whitespace, or any other kind of FASTA-like header/data.

Unless C<strict_case> is specified, the case of your patterns and the case of your input doesn't matter. I suggest using uppercase everywhere.






=head1 INTERBASE COORDINATES

All offsets returned by this module are in "interbase coordinates". Rather than the first base in a sequence being considered index 1 as most biologists might think of it, in interbase coordinates the first base is the sequence spanning coordinates 0 through 1.

One of the reasons this is useful is because it allows us to unambiguously specify 0-width sequences like for example endonuclease cut sites. If biology-style base coordinates are used it is ambiguous whether the cut is before or after.

Unlike with most string operations in perl, the start coordinate can actually be greater than the end coordinate if the pattern was found on the reverse complement strand and C<double_stranded> is set. For DNA, C<double_stranded> is the default so you should use C<single_stranded> if you don't want reverse complement matches. For RNA the default is C<single_stranded>.

For circular inputs, interbase coordinates can also be greater than the length of the input. This is interpreted as wrapping back around to the beginning in a modular arithmetic fashion. Similarly for negative coordinates. However, "Out-of-range" interbase coordinates are only defined for circular inputs and referencing them on linear inputs will throw errors.



=head1 EXHAUSTIVE SEARCH

Most methods of searching nucleotide sequences will only find non-overlapping matches in the input. For example, when searching for the sequence C<AA> in the input C<AAAA>, perl's C<m/AA/g> searches will only return 2 matches:

    AAAA
    --
      --

With this module you get all three matches that occur anywhere:

    AAAA
    --
     --
      --

For DNA data this can be useful for finding the comprehensive set of possible molecules that could exist after a restriction enzyme cleaving.



=head1 IUPAC SHORT FORMS

For DNA and RNA, IUPAC incompletely specified amino acid sequences can be used. These are analogous to regular expression character classes. Just like perl's C<\s> is short for C<[ \r\n\t]>, in IUPAC form C<V> is short for C<[ACG]>, or C<[^T]>. Unless C<strict_thymine_uracil> is in effect this will actually be like C<[^TU]> for both DNA and RNA inputs.




=head1 ADDING MULTIPLE SEARCH PATTERNS

One requirement for this module is that any number of regular expressions can be combined into one so that many patterns can be searched for simultaneously while doing a single pass over the data.

Doing a single pass is generally more efficient because of memory locality and has a few other positive side-effects. For instance, we can scan the reverse complement strand using a single pass over the sense strand and therefore avoid copying the input.

This module should be able to support quite a large number of simultaneous search patterns although I have some ideas for future optimisations. Large numbers of patterns may come in handy when building a list of all restriction enzymes that don't cut a target sequence, or finding all PCR primer sites accounting for IUPAC expanded primers.

Multiple patterns can be added at once simply by calling C<add()> multiple times before attempting a C<match> (or a C<compile>):

    my $re = Bio::Regexp->new;

    $re->add($_) for ('GAATTC', 'CCWGG');

    my @matches = $re->match($input);

Which pattern matched is returned as the C<match> key in the returned match results. You should probably have a hash of all your patterns so that you can look them up while processing matches. The way this is implemented is similar to the very useful L<Regexp::Assemble> except that it doesn't implement the hacks needed for ancient perl versions.

When matching, only a single pass will be made over the data in order to find all possible locations that either of the added sequences could have matched. Large numbers of patterns should be fairly efficient because the perl 5.10+ regular expression engine uses a trie data structure for such patterns (and 5.10 is the minimum required perl for other reasons).





=head1 CIRCULAR INPUTS

The search sequence C<GAATTC> will match the following input:

    ATTCGGGGGGGGGGGGGGGGGGA

The C<start> index in the match will be 22 and the end will be 26. Since the input's length is only 23, we know that it must have wrapped around. Note that you must call the C<circular> method on the object before attempting a C<match> if you wish to find such patterns.

In order to make this efficient even with really long input sequences, this module copies only the maximum length your search pattern could possibly be. Being able to figure out the minimum and maximum sequence lengths is one of the reasons why the types of regular expressions you can use with this module are limited.



=head1 SEE ALSO

L<Bio-Regexp github repo|https://github.com/hoytech/Bio-Regexp>

L<Bio::Tools::SeqPattern> from the BioPerl distribution also allows the manipulation of patterns but is less advanced than this module. Also, the way L<Bio::Tools::SeqPattern> reverses a regular expression in order to match the reverse complement is... wow. Just wow. :)

L<Bio::Grep> is an interface to various programs that search biological sequences. L<Bio::Grep::Backend::RE> is probably the most comparable to this module.

L<Bio::DNA::Incomplete>


=head1 AUTHOR

Doug Hoyte, C<< <doug@hcsw.org> >>


=head1 COPYRIGHT & LICENSE

Copyright 2013 Doug Hoyte.

This module is licensed under the same terms as perl itself.

=cut



TODO:

Implement back references in language? http://www.bioperl.org/wiki/Regular_expressions_and_Repeats
