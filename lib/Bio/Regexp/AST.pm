package Bio::Regexp::AST;

use common::sense;

use Data::Dumper;

{
  use Regexp::Grammars;

  $Bio::Regexp::AST::parser = qr{
    ## MAIN
    <regexp>

    ## GRAMMAR

    <objtoken: Bio::Regexp::AST::regexp>
    ^
      <[element]>*
    $

    <objtoken: Bio::Regexp::AST::element>
      <literal> <repeat>? |
      <charclass> <repeat>?

    <objtoken: Bio::Regexp::AST::literal>
      [a-zA-Z]

    <objtoken: Bio::Regexp::AST::charclass>
      \[ <negate_charclass>? <[literal]>+ \]

    <objtoken: Bio::Regexp::AST::negate_charclass>
      \^

    <objtoken: Bio::Regexp::AST::repeat>
      \{
         (?:
            <max=(?: \d+ )> <min=(?{ $MATCH{max} })> |
            <min=(?: \d+ )> , <max=(?: \d+ )>
         )
      \} |
      \? <min=(?{ 0 })> <max=(?{ 1 })>
  }xs;
}


sub Bio::Regexp::AST::regexp::render {
  my ($self, $level) = @_;

  my $output = '';

  for my $element (@{ $self->{element} }) {
    $output .= $element->render;
  }

  return $output;
}


sub Bio::Regexp::AST::regexp::reverse_complement {
  my ($self, $level) = @_;

  print Dumper($self);
  die;
}


sub Bio::Regexp::AST::regexp::compute_min_max {
  my ($self, $level) = @_;

  my $min = my $max = 0;

  for my $element (@{ $self->{element} }) {
     if (my $repeat = $element->{repeat}) {
       $min += $repeat->{min};
       $max += $repeat->{max};
     } else {
       $min++;
       $max++;
     }
  }

  return ($min, $max);
}



sub Bio::Regexp::AST::element::render {
  my ($self, $level) = @_;

  my $output = '';

  if ($self->{literal}) {
    $output .= $self->{literal}->render;
  } elsif ($self->{charclass}) {
    $output .= $self->{charclass}->render;
  } else {
    die "unknown element type";
  }

  $output .= $self->{repeat}->render if $self->{repeat};

  return $output;
}

sub Bio::Regexp::AST::literal::render {
  my ($self, $level) = @_;

  return $self->{''};
}

sub Bio::Regexp::AST::charclass::render {
  my ($self, $level) = @_;

  my $output = '';

  $output .= '[';

  $output .= '^' if $self->{negate_charclass};

  foreach my $literal (@{ $self->{literal} }) {
    $output .= $literal->render;
  }

  $output .= ']';

  return $output;
}

sub Bio::Regexp::AST::repeat::render {
  my ($self, $level) = @_;

  if ($self->{min} == 0 && $self->{max} == 1) {
    return '?';
  } elsif ($self->{min} == $self->{max}) {
    return "{$self->{min}}";
  } else {
    return "{$self->{min},$self->{max}}";
  }

  return $self->{''};
}



1;
