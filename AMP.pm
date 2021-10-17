#!/usr/bin/perl
package AMP;
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

sub calculaMM{#calcula e retorna o peso molecular da proteina em Daltons
	my($seq)=@_;#sequencia de aas
	chomp $seq;
	my $peso=0;#será incrementada com o peso de cada aa e será o retorno da função
	my %massa=(#massa dos aas em Daltons
		A => 89,
		R => 174,
		N => 132,
		D => 133,
		C => 121,
		Q => 146,
		E => 147,
		G => 75,
		H => 155,
		I => 131,
		L => 131,
		K => 146,
		M => 149,
		F => 165,
		P => 115,
		S => 105,
		T => 119,
		W => 204,
		Y => 181,
		V => 117
	);
	for(my $cont=0;$cont<length($seq);$cont++){	
		if(exists($massa{substr($seq,$cont,1)})){ 
			$peso+=$massa{substr($seq,$cont,1)};
		}else{
			print "sequencia com aminoácido desconhecido!!!";
		}
	}
	return $peso;
}
sub volumeVanDerWaals{#calcula o volume de Van der Waals
	my($seq)=@_;
	chomp $seq;
	my $vol=0;
	my $desconto=(length($seq)-1)*7;#Desconto para cada ligação peptídica
	my %volumes=(#volumes dos aas em angstrons³
		A => 67,
		R => 148,
		N => 96,
		D => 91,
		C => 86,
		Q => 114,
		E => 109,
		G => 48,
		H => 118,
		I => 124,
		L => 124,
		K => 135,
		M => 124,
		F => 135,
		P => 90,
		S => 73,
		T => 93,
		W => 163,
		Y => 141,
		V => 105
	);
	for(my $cont=0;$cont<length($seq);$cont++){
		$vol+=$volumes{substr($seq,$cont,1)};
	}
	return $vol-$desconto;#falta descontar as contantes da intersecção C-N
}
sub hydrophobicity{
	my($seq,$scale)=@_;
	my $tamanho=length($seq);
	my $hydro=0;
	my %hydrophobicity = ($scale eq "e") ? %Eisenberg : %KyteAndDoolittle;
	for(my $i=0;$i<$tamanho;$i++){
		$hydro+=$hydrophobicity{substr($seq,$i,1)};
	}
	$hydro/=$tamanho;
	return $hydro;
}
sub hydrophobicitymoment{
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
sub hydropercent{
	my($seq)=@_;
	my $tamanho=length($seq);
	my $hidrofobico=0;
	for(my $i=0;$i<$tamanho;$i++){ 
		my $aa=substr($seq,$i,1);
		if($aa=~/[AVLIFYW]/i){
			$hidrofobico++;
		}
	}
	return $hidrofobico/$tamanho;
}
sub ratio{
	my($seq)=@_;
	my $tamanho=length($seq);
	my $charged=0;
	my $hidrofobico=0;
	for(my $i=0;$i<$tamanho;$i++){
		my $aa=substr($seq,$i,1);
		if($aa=~/[DEHKR]/i){
			$charged++;
		}
		if($aa=~/[AVLIFYW]/i){
			$hidrofobico++;
		}
	}
	if($charged!=0){
		return $hidrofobico/$charged;
	}else{
		return "NaN";
	}
}
sub carga{
	my($seq)=@_;
	my $tamanho=length($seq);
	my $charge=0;
	for(my $i=0;$i<$tamanho;$i++){
		my $aa=substr($seq,$i,1);
		if($aa=~/[DE]/i){
			$charge--;
		}
		if($aa=~/[KRH]/i){
			$charge++;
		}
	}
	return $charge;
}
sub BomanIndex{#calcula o índice de Boman, dado pela média das energias livres dos aminoácidos
	my($seq)=@_;#sequencia de aas
	chomp $seq;
	my $index=0;#será incrementada com o index de cada aa e será o retorno da função
	for(my $cont=0;$cont<length($seq);$cont++){	
		if(exists($RadzickaWolfenden{substr($seq,$cont,1)})){ 
			$index+=$RadzickaWolfenden{substr($seq,$cont,1)};
		}
	}
	return -1*($index/length $seq);
}
sub Helix{#calcula a média de propensão de hélice,
	my($seq)=@_;#sequencia de aas
	chomp $seq;
	my $index=0;#será incrementada com o index de cada aa e será o retorno da função
	my %freeenergy=(
		A => 0,
		C => 0.68,
		D => 0.69,
		E => 0.40,
		F => 0.54,
		G => 1,
		H => 0.66,
		I => 0.41,
		K => 0.26,
		L => 0.21,
		M => 0.24,
		N => 0.65,
		P => 3.16,
		Q => 0.39,
		R => 0.21,
		S => 0.5,
		T => 0.66,
		V => 0.61,
		W => 0.49,
		Y => 0.53
	);
	for(my $cont=0;$cont<length($seq);$cont++){	
		if(exists($freeenergy{substr($seq,$cont,1)})){ 
			$index+=$freeenergy{substr($seq,$cont,1)};
		}else{
			print "sequencia com aminoácido desconhecido!!! -> ".substr($seq,$cont,1)."\n\n";
		}
	}
	return $index/length $seq;
}

sub geoMeanHM{
	my ($seq)=@_;
	return (hydrophobicitymoment($seq,"kt")*hydrophobicitymoment($seq,"rw")*hydrophobicitymoment($seq,"e"))**(1/3);
}

1;
