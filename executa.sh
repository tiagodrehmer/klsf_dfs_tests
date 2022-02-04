MULTIPLICADOR_ITERATION="$1"
MAX_IN_SOLUTION="$2"
FRAC_FIX="$3"


SEED=1
while [ $SEED -lt 31 ];
do
   echo "$INSTANCE"
   for INSTANCE in klsf_instances/*;
   do 
      PART=0
         echo "./p $INSTANCE $SEED $MULTIPLICADOR_ITERATION $MAX_IN_SOLUTION $FRAC_FIX $FRAC_SUB_FICA $MAX_POND 1"
         ./p $INSTANCE $SEED $MULTIPLICADOR_ITERATION $MAX_IN_SOLUTION $FRAC_FIX 1
         let PART=PART+1;
      echo "$SEED"
   done
   let SEED=SEED+1;
done
      
