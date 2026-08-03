#!/bin/bash
#BSUB -o /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/work_dir/%J.o
#BSUB -e /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/work_dir/%J.e
#BSUB -M 4000 -R "select[mem>4000] rusage[mem=4000]" -M 4000
#BSUB -q normal
#BSUB -n 2
 
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms8G -Xmx8G -Dnxf.pool.maxThreads=2000"
export NXF_VER=22.04.0-5697

# Load Nextflow module
module load nextflow/25.04.1-5946
  
# Run the Nextflow script with the specified parameters
nextflow run \
/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/Talleh_code/singleR_scripts_and_refs/submsn_singleR.nf \
-w /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/work_dir/work_singleR \
-profile sanger \
-c /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/work_dir/profiles.config \
-qs 1000 \
-resume \
--softpack_module "HGI/softpack/groups/team222/Pf_scRNAseq/34" \
--scrpt /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/Talleh_code/singleR_scripts_and_refs/singleR_mult.R \
--input_seu "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/data/processed/Pf/MSC*/seu_obj/dense_CbendrCr.RDS" \
--ref_obj_p "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/Talleh_code/singleR_scripts_and_refs/annot_refs/msc_lp_sce.RDS" \
--marker_obj_p "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/Talleh_code/singleR_scripts_and_refs/annot_refs/com_markers_localMac.RDS" \
--ref_stage "stageHL" \
--o_sufx "pred_sR" \
--ncores "5" \
--mem "70 GB" \
--run_time "3.hour" \
--resume \

# ## clean up on exit 0 - delete this if you want to keep the work dir
# status=$?
# if [[ $status -eq 0 ]]; then
#   rm -r /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/work_dir/work_singleR
# fi

