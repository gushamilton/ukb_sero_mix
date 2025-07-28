# This script runs the R script to create the GWAS input files and then uploads
# them to the specified location on the DNAnexus platform. THIS SHOULD BE RUN HERE

dx run app-cloud_workstation --ssh \
  -isnapshot=r_env \
    -imax_session_length="24h" \
  --instance-type mem1_ssd1_v2_x4 \
   -y --allow-ssh 0.0.0.0/16

DX_PROJECT_CONTEXT_ID=project-J1GYbG0JQz0gxQ6yZf1GqYkb
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:



dx extract_dataset app81499_20250623113809.dataset --fields 'participant.eid,participant.p23000_i0,participant.p23001_i0,participant.p23049_i0,participant.p23048_i0,participant.p23026_i0,participant.p23039_i0,participant.p23043_i0,participant.p23018_i0,participant.p23030_i0,participant.p23031_i0,participant.p23006_i0,participant.p23004_i0,participant.p23042_i0,participant.p23016_i0,participant.p23017_i0,participant.p23025_i0,participant.p23024_i0,participant.p23023_i0,participant.p23022_i0,participant.p23010_i0,participant.p23011_i0,participant.p23027_i0,participant.p23015_i0,participant.p23029_i0,participant.p23032_i0,participant.p23014_i0,participant.p23028_i0,participant.p23019_i0,participant.p23041_i0,participant.p23037_i0,participant.p23013_i0,participant.p23044_i0,participant.p23003_i0,participant.p23040_i0,participant.p23005_i0,participant.p23002_i0,participant.p23034_i0,participant.p23033_i0,participant.p23012_i0,participant.p23020_i0,participant.p23038_i0,participant.p23009_i0,participant.p23008_i0,participant.p23007_i0,participant.p23021_i0,participant.p23035_i0,participant.p23036_i0,participant.p23050_i0,participant.p23051_i0,participant.p23052_i0,participant.p23053_i0,participant.p23054_i0,participant.p23055_i0,participant.p23056_i0,participant.p23057_i0,participant.p23058_i0,participant.p23059_i0,participant.p23060_i0,participant.p23061_i0,participant.p23062_i0,participant.p23063_i0,participant.p23064_i0,participant.p23065_i0,participant.p23066_i0,participant.p23067_i0,participant.p23068_i0,participant.p23075_i0,participant.p23069_i0,participant.p23070_i0,participant.p23071_i0,participant.p23073_i0,participant.p23074_i0,participant.p50_i0,participant.p21001_i0,participant.p4080_i0_a0,participant.p20258,participant.p49_i0' --delimiter ' ' --output phenotypes_serology.txt

dx extract_dataset app81499_20250623113809.dataset --fields 'participant.eid,participant.p23012_i0,participant.p4080_i0_a0,participant.p20258,participant.p49_i0' --delimiter ' ' --output phenotypes_serology.txt


dx run table-exporter \
  -idataset_or_cohort_or_dashboard=app81499_20250623113809.dataset \
  -ifield_names_file_txt=serology_fields.tsv \
  -ioutput_format=TSV \
  -ientity=participant \
  -ioutput=serology_export \
  --yes