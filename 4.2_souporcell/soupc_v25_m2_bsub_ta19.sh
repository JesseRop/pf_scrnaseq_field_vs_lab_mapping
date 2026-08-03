#!/bin/bash
#BSUB -o /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/work_dir/%J.o
#BSUB -e /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/work_dir/%J.e
#BSUB -M 4000 -R "select[mem>4000] rusage[mem=4000]" -M 4000
#BSUB -q oversubscribed
#BSUB -n 3
 
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms8G -Xmx8G -Dnxf.pool.maxThreads=2000"
export NXF_VER=22.04.0-5697
  
# Load Nextflow module
module load nextflow/25.04.1-5946
  
# Run the Nextflow script with the specified parameters
nextflow run \
/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/scripts_nbooks/souporcell/soupc_v25_m2_ta19.nf \
-w /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/work_dir/work_soupc \
-profile sanger \
-c /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/work_dir/profiles.config \
-qs 1000 \
-resume \
--id_decode "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/data/raw/pf_all_id_decode.csv" \
--bam "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf_Talleh/5736STDY*/outs/possorted_genome_bam.bam" \
--bcodes "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/data/processed/Pf/MSC*/soupc_resc_hspf/barcodes.tsv.gz" \
--soup_dir "soupc_resc_hspf" \
--o_dir "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/data/processed/Pf/" \
--scrpt "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/scripts_nbooks/souporcell/soupc_v25_mnmap_hsat_310_ta19.sh" \
--ref_file "PlasmoDB-68_Pfalciparum3D7_Genome_VAR_masked.fasta" \
--hsat_ref_file "Pf68" \
--ncores "15" \
--mem "100 GB" \
--run_time "12.hour"
  
# ## clean up on exit 0 - delete this if you want to keep the work dir
status=$?
if [[ $status -eq 0 ]]; then
  rm -r /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Talleh/work_dir/work_soupc
fi

