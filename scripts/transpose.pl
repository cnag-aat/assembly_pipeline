#!/usr/bin/env perl

while(<>){
        $_=~ tr/\n//d;
        @arr=split("\t",$_);
        $col=$#arr;
        for($i=0;$i<=$#arr;$i++){
                $index=sprintf("%s%s",$.,$i);
                $hash{$index}=$arr[$i];
        }
        $row=$.;
}

for($a=0;$a<=$col;$a++){
        for($b=1;$b<=$row;$b++){
                $t=sprintf("%s%s",$b,$a);
                print $hash{$t},"\t";
        }
        print "\n";
}
