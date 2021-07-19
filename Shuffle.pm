#!/usr/bin/perl
# Created on: 21/06/2011
#Shuffe de amion�cidos
package Shuffe;
use warnings;
use strict;
use HydrophobicMoment;

sub verify{
	my($number,@sorteados)=@_;
	for(my $i=0;$i<@sorteados;$i++){
		if($number==$sorteados[$i]){
			return 1;
		}
	}
	return 0;
}

sub shuffle{
	my ($seq)=@_;
	my @tset=();
	my $shuffle="";
	for(my $i=0;$i<length $seq;$i++){
		my $valor=int rand length $seq;
		if(verify($valor, @tset)==0){
			$tset[$i]=$valor;
			$shuffle.=substr($seq,$valor,1);
		}else{
			$i--;
		}
	}
	return $shuffle;
}

sub averageHM{
	my($sequence,$repetitions)=@_;
	my $average=0;
	for(my $i=0;$i<$repetitions;$i++){
		$average+=AMP::geoMeanHM(shuffle($sequence));#hydrophobicitymoment(shuffle($sequence),"kt");#
	}
	return $average/$repetitions;
}

1;
