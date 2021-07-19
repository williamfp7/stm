#!/usr/bin/perl
package HydrophobicMoment;
use constant PI => 3.14159;
my %KyteAndDoolittle=(
	A => 1.8,
	R => -4.5,
	N => -3.5,
	D => -3.5,
	C => 2.5,
	Q => -3.5,
	E => -3.5,
	G => -0.4,
	H => -3.2,
	I => 4.5,
	L => 3.8,
	K => -3.9,
	M => 1.9,
	F => 2.8,
	P => -1.6,
	S => -0.8,
	T => -0.7,
	W => -0.9,
	Y => -1.3,
	V => 4.2
);

my %Eisenberg=(
	A => 0.25,
	R => -1.8,
	N => -0.64,
	D => -0.72,
	C => 0.04,
	Q => -0.69,
	E => -0.62,
	G => 0.16,
	H => -0.4,
	I => 0.73,
	L => 0.53,
	K => -1.1,
	M => 0.26,
	F => 0.61,
	P => -0.07,
	S => -0.26,
	T => -0.18,
	W => 0.37,
	Y => 0.02,
	V => 0.54
);

my %RadzickaWolfenden=(
	A =>   1.81,
	C =>   1.27,
	D =>  -8.72,
	E =>  -6.81,
	F =>   2.97,
	G =>   0.93,
	H =>  -4.66,
	I =>   4.92,
	K =>  -5.55,
	L =>   4.92,
	M =>   2.35,
	N =>  -6.64,
	P =>   0.00,
	Q =>  -5.54,
	R => -14.92,
	S =>  -3.40,
	T =>  -2.75,
	V =>   4.04,
	W =>   2.32,
	Y =>  -0.14
);

sub hydrophobicMoment{
	my($seq,$scale)=@_;
	my $tamanho=length($seq);
	my %hydrophobicity = ($scale eq "e") ? %Eisenberg : ($scale eq "rw") ? %RadzickaWolfenden : %KyteAndDoolittle;
	my $sen=0;
	my $cos=0;
	for(my $i=0;$i<$tamanho;$i++){
		$sen+=$hydrophobicity{substr($seq,$i,1)}*sin(($i+2)*100*PI/180);
		$cos+=$hydrophobicity{substr($seq,$i,1)}*cos(($i+2)*100*PI/180);
	}
	return (sqrt $sen*$sen+$cos*$cos)/$tamanho;
}

sub geoMeanHM{
	my ($seq)=@_;
	return (hydrophobicMoment($seq,"kt")*hydrophobicMoment($seq,"rw")*hydrophobicMoment($seq,"e"))**(1/3);
}

1;
