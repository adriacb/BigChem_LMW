#!/usr/bin/perl
#
# Execute as:
#   filter_dist.pl ref_file sd_file > output
#
   if ($#ARGV < 0) {
     printf "\n";
     printf "To Filter a SD file based on distances, execute as:\n";
     printf "   filter_dist.pl REF_FILE SD_FILE > OUTPUT\n";
     printf "\n";
     printf "To get a PDB file of the cube(s) to be matched, execute as:\n";
     printf "   filter_dist.pl REF_FILE > OUTPUT\n";
     printf "\n";
     die;
   }
# Step 1: Read the file with centers to be matched ($ARGV[0])
# Format: x y z tol
#
   open(F1,"$ARGV[0]") || die "Cannot open file: $!";
   while (<F1>) {
     if (/^#/) {next}
     push (@center,$_);
   }
   close(F1);
   for $j (0..$#center) {
     @line = split(' ',$center[$j]);
     $tol[$j] = $line[3];
     $limx1[$j] = $line[0]-$tol[$j];
     $limy1[$j] = $line[1]-$tol[$j];
     $limz1[$j] = $line[2]-$tol[$j];
     $limx2[$j] = $line[0]+$tol[$j];
     $limy2[$j] = $line[1]+$tol[$j];
     $limz2[$j] = $line[2]+$tol[$j];
   }
#
   if ($#ARGV == 0) {
     printf "HEADER  Any molecule containing at least 1 atom\n";
     printf "HEADER  within this Cube(s) will be accepted\n";
     for $j (0..$#center) {
      printf "HEADER  Tolerance = %4.1fA\n",$tol[$j];
      printf "ATOM   %4i  W   XXX    %2i    %8.3f%8.3f%8.3f\n",($j*8)+1,$j+1,$limx1[$j],$limy1[$j],$limz1[$j];
      printf "ATOM   %4i  W   XXX    %2i    %8.3f%8.3f%8.3f\n",($j*8)+2,$j+1,$limx1[$j],$limy1[$j],$limz2[$j];
      printf "ATOM   %4i  W   XXX    %2i    %8.3f%8.3f%8.3f\n",($j*8)+3,$j+1,$limx1[$j],$limy2[$j],$limz1[$j];
      printf "ATOM   %4i  W   XXX    %2i    %8.3f%8.3f%8.3f\n",($j*8)+4,$j+1,$limx1[$j],$limy2[$j],$limz2[$j];
      printf "ATOM   %4i  W   XXX    %2i    %8.3f%8.3f%8.3f\n",($j*8)+5,$j+1,$limx2[$j],$limy2[$j],$limz2[$j];
      printf "ATOM   %4i  W   XXX    %2i    %8.3f%8.3f%8.3f\n",($j*8)+6,$j+1,$limx2[$j],$limy2[$j],$limz1[$j];
      printf "ATOM   %4i  W   XXX    %2i    %8.3f%8.3f%8.3f\n",($j*8)+7,$j+1,$limx2[$j],$limy1[$j],$limz2[$j];
      printf "ATOM   %4i  W   XXX    %2i    %8.3f%8.3f%8.3f\n",($j*8)+8,$j+1,$limx2[$j],$limy1[$j],$limz1[$j];
     }
     for $j (0..$#center) {
      printf "CONECT %4i %4i %4i %4i\n",($j*8)+1,($j*8)+2,($j*8)+3,($j*8)+8;
      printf "CONECT %4i %4i %4i\n",($j*8)+2,($j*8)+4,($j*8)+7;
      printf "CONECT %4i %4i %4i\n",($j*8)+3,($j*8)+4,($j*8)+6;
      printf "CONECT %4i %4i %4i %4i\n",($j*8)+5,($j*8)+4,($j*8)+6,($j*8)+7;
      printf "CONECT %4i %4i %4i\n",($j*8)+8,($j*8)+6,($j*8)+7;
     }
     printf "END\n";
     die;
   }
#
#
# Read SD file
#
   open(F0,"$ARGV[1]") || die "Cannot open file: $!";
   while (<F0>) {
    if (/^\$\$/) {
     push (@ENTRY,$_);
     @line4 = substr ($ENTRY[3],0,3);
#    @line4 = split(' ',$ENTRY[3]); # FAILS FOR N > 99
     $natoms = $line4[0];
#
     BLOCK1:
     for $i (4..$natoms+3) {
# skip non_H
      $ele = substr ($ENTRY[$i],31,1);
      if ($ele eq 'H' || $ele eq 'O' || $ele eq 'N') {next BLOCK1}
#
      @line = split(' ',$ENTRY[$i]);
      $x = $line[0];
      $y = $line[1];
      $z = $line[2];
#
      BLOCK2:
      for $j (0..$#center) {
      if ($limx1[$j] <= $x && $x <= $limx2[$j] &&
          $limy1[$j] <= $y && $y <= $limy2[$j] &&
          $limz1[$j] <= $z && $z <= $limz2[$j]) { $save[$j] = 1 }
      } # End of BLOCK2
     } # End of BLOCK1
#
     $matched = 0;
     for $j (0..$#center) {
       $matched += $save[$j];
       $save[$j] = 0;
     }
     if ($matched == $#center+1) {
      while (@ENTRY) {
       $outline = shift (@ENTRY);
       print $outline;
      }
     }
     undef(@ENTRY);
    }else{
      push (@ENTRY,$_);
    }
   }
   close(F0);
