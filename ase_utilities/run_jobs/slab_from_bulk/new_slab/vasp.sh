#!/bin/bash
#SBATCH --job-name=vasp_job
#SBATCH --partition razi
#SBATCH --nodes=1
#SBATCH --mem=20G   
#SBATCH --ntasks-per-node=8
#SBATCH --time=2:00:00



####SBATCH --output $(JOBNAME)
##SBATCH --mem-per-cpu=3200mb
###SBATHC -j oe
echo " "
echo "### Setting up shell environment ..."
echo " "
if test -e "/etc/profile"; then source "/etc/profile"; fi;
if test -e "$HOME/.bash_profile"; then source "$HOME/.bash_profile"; fi;
unset LANG; export LC_ALL="C"; export MKL_NUM_THREADS=1; export OMP_NUM_THREADS=1
export USER=${USER:=`logname`}
export MOAB_JOBID=${SLURM_JOB_ID:=`date +%s`}
export MOAB_SUBMITDIR=${SLURM_SUBMIT_DIR:=`pwd`}
export MOAB_JOBNAME=${SLURM_JOB_NAME:=`basename "$0"`}
export MOAB_JOBNAME=$(echo "${SLURM_JOB_NAME}" | sed 's/[^a-zA-Z0-9._-]/_/g')
export MOAB_NODECOUNT=${SLURM_JOB_NUM_NODES:=1}
export MOAB_PROCCOUNT=${SLURM_NTASKS:=1}
ulimit -s unlimited


echo " "
echo "### Printing basic job infos to stdout ..."
echo " "
echo "START_TIME           = `date +'%y-%m-%d %H:%M:%S %s'`"
echo "HOSTNAME             = ${HOSTNAME}"
echo "USER                 = ${USER}"
echo "MOAB_JOBNAME         = ${SLURM_JOB_NAME}"
echo "MOAB_JOBID           = ${SLURM_JOB_ID}"
echo "MOAB_SUBMITDIR       = ${SLURM_SUBMIT_DIR}"
echo "MOAB_NODECOUNT       = ${SLURM_JOB_NUM_NODES}"
echo "MOAB_PROCCOUNT       = ${SLURM_NTASKS}"
echo "SLURM_NODELIST       = ${SLURM_NODELIST}"
echo "PBS_NODEFILE         = ${SLURM_JOB_NODELIST}"
if test -f "${SLURM_JOB_NODELIST}"; then
  echo "PBS_NODEFILE (begin) ---------------------------------"
  cat "${SLURM_JOB_NODELIST}"
  echo "PBS_NODEFILE (end) -----------------------------------"
fi
echo "---------------- ulimit -a -S ----------------"
ulimit -a -S
echo "---------------- ulimit -a -H ----------------"
ulimit -a -H
echo "----------------------------------------------"

module unload gcc

module load gcc/8.4.0
module load intel/2019.3
module load openmpi/4.1.1-gnu
module load vasp/5.4.4-ompi411

echo " "
echo "### Creating TMP_WORK_DIR directory and changing to it ..."
echo " "
#JOB_WORK_DIR="${SLURM_JOB_NAME}.uc1.${SLURM_JOB_ID##*.}.$(date +%y%m%d_%H%M%S)"
#if test -z "$SLURM_JOB_NUM_NODES" -o "$SLURM_JOB_NUM_NODES" = "1"; then
# # TMP_BASE_DIR="/tmp/${USER}"
#  ws_allocate "$SLURM_JOB_ID" 10
#  TMP_BASE_DIR="$(ws_find $SLURM_JOB_ID)"
#else
#  # In case of 2 or more nodes, use a common scratch dir available on all nodes...
#  ws_allocate "$SLURM_JOB_ID" 10
#  TMP_BASE_DIR="$(ws_find $SLURM_JOB_ID)"
#  #TMP_BASE_DIR="$(ws_list | grep $JOBID|tail -1)"
#fi
#TMP_WORK_DIR="${TMP_BASE_DIR}/${JOB_WORK_DIR}"
#ln -s $TMP_WORK_DIR $SLURM_SUBMIT_DIR/scratch_dir
#RUN_DIR='/scratch'
TMP_WORK_DIR=/scratch/${SLURM_JOB_ID}
JOB_WORK_DIR="$(basename $TMP_WORK_DIR)"
TMP_BASE_DIR="$(dirname $TMP_WORK_DIR)"



echo "JOB_WORK_DIR         = ${JOB_WORK_DIR}"
echo "TMP_BASE_DIR         = ${TMP_BASE_DIR}"
echo "TMP_WORK_DIR         = ${TMP_WORK_DIR}"
#mkdir -vp "${TMP_WORK_DIR}"
#echo "Copying queueing system script to TMP_WORK_DIR ..."
#cp -v "$0" "${TMP_WORK_DIR}/"
ln -s "${TMP_WORK_DIR}" scratch_dir
cd "${TMP_WORK_DIR}"


#echo " "
#echo "### Loading software module:"
#echo " "
## Load a specific software version:
#if test -z "$VASP_VERSION"; then
#  echo "ERROR: Failed to load module 'chem/vasp/5.3.3.4'."
#  exit 101
#fi
#echo "VASP_VERSION         = $VASP_VERSION"
#module list



echo " "
echo "### Copying input files for job (if required):"
echo " "
#cp -v $SLURM_SUBMIT_DIR/{INCAR,KPOINTS,POSCAR,POTCAR,vdw_kernel.bindat}        $TMP_WORK_DIR/
cp -v $SLURM_SUBMIT_DIR/{*py,*cif,POSCAR,coord,*xyz,in.traj,CONTCAR,*.pkl}        $TMP_WORK_DIR/

echo " "



script_name=$BASH_SOURCE

var_loop=0
echo ""
echo ""
while IFS= read -r line
do
  if [[ "$line" = *"--time"* ]]; then
     timeslurm=$(echo $line | sed 's/#SBATCH --time=//g')
     slurm_minutes=$(echo $timeslurm | awk -F ":" '{print $(NF-1)}')
     slurm_seconds=$(echo $timeslurm | awk -F ":" '{print $(NF)}')
     if [[ $(echo $timeslurm | awk -F ":" '{print (NF)}') -ne 3 ]]; then
       slurm_hours=0
     else
       slurm_hours=$(echo $timeslurm | awk -F ":" '{print $(NF-2)}')
     fi
     echo "Running for $(echo "$slurm_hours*1" |bc)h $(echo "$slurm_minutes*1" |bc)m and $(echo "$slurm_seconds*1" |bc)s."
     timeslurm=$(echo "$slurm_hours*3600 + $slurm_minutes*60 + $slurm_seconds" | bc)
     echo "This means $timeslurm seconds."
     timeslurm=$(echo "$timeslurm *0.95" |bc)
     echo "Will terminate at ${timeslurm}s to copy back necessary files from scratch"
  fi
  var_loop=$((var_loop+1))
  if [[ $var_loop = 10 ]]; then
    break
  fi
done < "$script_name"

echo ""
echo ""


time timeout ${timeslurm}s ./run.py $MOAB_SUBMITDIR  > $MOAB_SUBMITDIR/out.txt

exit_status=$?
if [ $exit_status -eq 124 ]; then
  echo " "
  echo "Cancelled due to timelimit."
  echo " "
fi


#mpirun -bootstrap slurm -n ${SLURM_NTASKS} /home/fh2-project-catfut1/xl0762/VASP/VASP_BIN/vasp_std
#time ./run.py $MOAB_SUBMITDIR  > $MOAB_SUBMITDIR/out.txt
#time mpirun -np ${SLURM_NTASKS} /home/fh2-project-catfut1/xl0762/vasp/vasp_std
exit_code=$?
echo " "
echo "### Cleaning up files ... removing unnecessary scratch files ..."
echo " "
rm -vf EIGENVAL* IBZKPT* PCDAT* PROCAR* vasprun.xml \
       ELFCAR* LOCPOT* OUTPAR* PROOUT* WAVECAR \
       TMPCAR* CHG* POTCAR vasp.dipcor
sleep 10 # Sleep some time so potential stale nfs handles can disappear.


echo " "
echo "### Compressing results and copying back result archive ..."
echo " "
cd "${TMP_BASE_DIR}"
mkdir -vp "${SLURM_SUBMIT_DIR}" # if user has deleted or moved the submit dir
echo "Creating result tgz-file '${SLURM_SUBMIT_DIR}/${JOB_WORK_DIR}.tgz' ..."
tar -zcvf "${SLURM_SUBMIT_DIR}/${JOB_WORK_DIR}.tgz" "${JOB_WORK_DIR}" \
  || { echo "ERROR: Failed to create tgz-file. Please cleanup TMP_WORK_DIR '$TMP_WORK_DIR' on host '$HOSTNAME' manually (if not done automatically by queueing system)."; exit 102; }


echo " "
echo "### Final cleanup: Remove TMP_WORK_DIR ..."
echo " "
rm -rvf "${TMP_WORK_DIR}"
cd ${SLURM_SUBMIT_DIR}
tar -xzf ${JOB_WORK_DIR}.tgz
mv ${JOB_WORK_DIR}/* .
rm -r ${JOB_WORK_DIR}.tgz ${JOB_WORK_DIR}

ws_release "$SLURM_JOB_ID"
rm $SLURM_SUBMIT_DIR/scratch_dir
echo "END_TIME             = `date +'%y-%m-%d %H:%M:%S %s'`"

echo "${SLURM_JOB_ID} is complete: on `date +'%y.%m.%d %H:%M:%S'` ${SLURM_SUBMIT_DIR}" >> ~/job.log

echo " "
echo "### Exiting with exit code '$exit_code' ..."
echo " "
exit $exit_code


