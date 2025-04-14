source /home/rdu230/.bashrc;
conda activate /home/rdu230/miniconda3/envs/d3tales_fw;
export DB_INFO_FILE=/scratch/rdu230/d3tales/high_throughput/config/db_infos.json
export PYTHONPATH=$PYTHONPATH:/scratch/rdu230/d3tales/high_throughput
export FW_CONFIG_FILE=/scratch/rdu230/d3tales/high_throughput/config/FW_config.yaml
cd /scratch/rdu230/d3tales/high_throughput