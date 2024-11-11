# ARGS1: nametag

nametag=$1

mkdir -p configs

cat config.4th.yaml | sed "s#PLACEHOLDER#$nametag#g" > configs/config.$nametag.yaml

for i in `seq 1 22`
do
  if [[ -f logs/run_2nd_${nametag}_${i}.log ]]
  then
    e=`cat logs/run_2nd_${nametag}_${i}.log | tail -n 1 | grep 'FileExistsError\|FileNotFoundError\|unlock'`
    if [[ ! -z $e ]]
    then
      sbatch --export=CHR=$i,NAME=$nametag --job-name=${i}_${nametag} run_2nd.sbatch 
      #qsub -v CHR=$i,NAME=$nametag -N ${i}_${nametag}_gwas run_2nd.qsub
    fi
  else
    # :
    sbatch --export=CHR=$i,NAME=$nametag --job-name=${i}_${nametag} run_2nd.sbatch
    #qsub -v CHR=$i,NAME=$nametag -N ${i}_${nametag}_gwas run_2nd.qsub
  fi
  sleep 30
done


# bash 03-run-gwas.sh 10000

