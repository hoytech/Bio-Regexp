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

Bio::Regexp - Exhaustive DNA/RNA/Protein regexp searches

=head1 SYNOPSIS

    my @matches = Bio::Regexp->new
                             ->add(q{AGC[2]YY[GTZ]{3,5}GCGC})
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

=head1 SEE ALSO

L<Bio-Regexp github repo|https://github.com/hoytech/Bio-Regexp>

=head1 AUTHOR

Doug Hoyte, C<< <doug@hcsw.org> >>

=head1 COPYRIGHT & LICENSE

Copyright 2013 Doug Hoyte.

This module is licensed under the same terms as perl itself.
