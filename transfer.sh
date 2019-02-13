SIM_DIR=$1
SIM_ID="MIS_BG"

mkdir $SIM_DIR
scp 	poincare:/gpfs1l/gpfshome/mds/grptraining/training12/sim/$SIM_ID.{err,log}\
	poincare:/gpfs1l/gpfshome/mds/grptraining/training12/sim/job.sh\
	poincare:/gpfs1l/gpfshome/mds/grptraining/training12/sim/output.csv\
	poincare:/gpfs1l/gpfshome/mds/grptraining/training12/sim/.git/refs/heads/master\
	$SIM_DIR
