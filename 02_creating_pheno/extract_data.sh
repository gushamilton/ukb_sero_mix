
dx run table-exporter \
  -idataset_or_cohort_or_dashboard=app81499_20250623113809.dataset \
  -ifield_names_file_txt=serology_fields.tsv \
  -ioutput_format=TSV \
  -ientity=participant \
  -ioutput=serology_export_title \
  -iheader_style=FIELD-TITLE


dx run app-cloud_workstation --ssh \
  -isnapshot=r_env \
    -imax_session_length="24h" \
  --instance-type mem1_ssd1_v2_x4 \
   -y --allow-ssh 0.0.0.0/16


DX_PROJECT_CONTEXT_ID=project-J1GYbG0JQz0gxQ6yZf1GqYkb
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:


awk '{$1=$1}1' covariates.txt > covariates_clean.txt
