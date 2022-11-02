#!/bin/bash
#
#  *################################## DFTD3 INFO ###############################################*
#  #                                                                                             #
#  # Use this script when VASP assigns the wrong functional for dispersion corrections.          #
#  # This happens at least with hse06 and pbe0, where pbe is chosen by VASP as default.          #
#  #                                                                                             #
#  #   -v    for verbose (<=1 errors; <=2 warnings; <=3 info; <=4 debugging)                     #
#  #                                                                                             #
#  #   -func for functional type (pbe, hse06,...)                                                #
#  #                                                                                             #
#  #   -type for type of d3 corrections (-zero, -bj, ...)                                        #
#  #                                                                                             #
#  #   -traj for the trajectory (e.g. POSCAR)                                                    #
#  #                                                                                             #
#  #   -pbc  to read periodic structures such as VASP POSCARS - This has to be the last argument #
#  *#############################################################################################*
#

dbg() {              if test "$DFTD3_VERB" -ge 4; then echo "DFTD3: dbg: $*" 1>&2; fi; }
inf() {              if test "$DFTD3_VERB" -ge 3; then echo "DFTD3: inf: $*" 1>&2; fi; }
wrn() {              if test "$DFTD3_VERB" -ge 2; then echo "DFTD3: WARNING: $*" 1>&2; fi; }
kill_DFTD3() { local ex=$?; if test "$DFTD3_VERB" -ge 1; then echo "DFTD3: ERROR: $* (exit code $ex)" 1>&2; fi; if test $ex -eq 0; then ex=199; fi; exit $ex; }

write_to_file() {        echo "$*" > final_corrected.e; }

test -z "$DFTD3_VERB"    && DFTD3_VERB="2"
test -z "$DFTD3_FUNC"    && DFTD3_FUNC="hse06"
test -z "$D3TYPE"        && D3TYPE="-zero"
test -z "$TRAJECTORY"    && TRAJECTORY="POSCAR"
test -z "$PBC"           && PBC=""

while test $# -gt 0; do
  case "$1" in
    (-v)      DFTD3_VERB="$2"; shift;;
    (-func)   DFTD3_FUNC="$2"; shift;;
    (-type)       D3TYPE="$2"; shift;;  
    (-traj)   TRAJECTORY="$2"; shift;;
    (-pbc)           PBC="-pbc"; shift;;
    (-[Hh?]*) awk 'BEGIN{f=0} /^[^#]/{f=1} /^ *$/{f=1} /END COMMENTS/{f=1} {if(f==0&&NR>1){print}}' "$0" | sed 's/^#/ /g'; exit 0;;
    (*)       kill_DFTD3 "Unkown command line argument '$1'. Please call 'dftd3.sh -help' for help.";;
  esac
  shift
done


test -z $PBC  && periodicity='no periodicity' || periodicity='under periodic boundary conditions'


inf "Verbosity $DFTD3_VERB selected."
inf "Functional $DFTD3_FUNC selected with type $D3TYPE ."
inf "DFTD3 will read $TRAJECTORY $periodicity ."



execute_pbe="dftd3 $TRAJECTORY $PBC $D3TYPE -func pbe"
execute_dftd3="dftd3 $TRAJECTORY $PBC $D3TYPE -func $DFTD3_FUNC"

if test -z $PBC; then
  wrn "Command will be created without periodic boudary conditions. Make sure this is really what you want. For a VASP calculation (e.g. POSCAR) this is needed."
fi

if test ! -e "final.e"; then
  kill_DFTD3 "final.e not found. WAIT until the calculation is over, this file is needed to obtain the correct energy."
fi

if [ "$DFTD3_FUNC" = "pbe" ]; then
  kill_DFTD3 "Your chosen functional is pbe, you do not need this script."
fi

dbg "execute_pbe command   : $(echo $execute_pbe)"
dbg "execute_dftd3 command : $(echo $execute_dftd3)"


$execute_pbe >dftd3_pbe.log        ## for PBE-D3
$execute_dftd3 >dftd3_$DFTD3_FUNC.log  ## for chosen functional

pbedisp=$(grep "Edisp" dftd3_pbe.log | awk '{print $5}')
funcdisp=$(grep "Edisp" dftd3_$DFTD3_FUNC.log | awk '{print $5}')

final_energy=$(echo "$(cat final.e) - $pbedisp + $funcdisp" | bc -l)
nodisp_energy=$(echo "$(cat final.e) - $pbedisp" | bc -l)
pbedisp_energy=$(cat final.e)

functional=$(echo $DFTD3_FUNC | tr '[:lower:]' '[:upper:]') 

write_to_file "final energy without dispersion : $nodisp_energy eV
final energy with PBE dispersion : $pbedisp_energy eV
final energy with $functional dispersion : $final_energy eV
"
