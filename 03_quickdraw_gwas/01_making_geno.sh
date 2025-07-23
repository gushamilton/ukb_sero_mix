dx mkdir -p QuickdrawsGWAS/

for CHR in {1..22}; do
dx run app-swiss-army-knife \
  --folder="QuickdrawsGWAS/" \
  -iin="Bulk/Genotype Results/Genotype calls/ukb22418_c${CHR}_b0_v2.bed" \
  -iin="Bulk/Genotype Results/Genotype calls/ukb22418_c${CHR}_b0_v2.bim" \
  -iin="Bulk/Genotype Results/Genotype calls/ukb22418_c${CHR}_b0_v2.fam" \
  -icmd="cp ukb22418_c${CHR}_b0_v2.bed copy_ukb22418_c${CHR}_b0_v2.bed && cp ukb22418_c${CHR}_b0_v2.bim copy_ukb22418_c${CHR}_b0_v2.bim && cp ukb22418_c${CHR}_b0_v2.fam copy_ukb22418_c${CHR}_b0_v2.fam" \
  --instance-type mem1_ssd1_v2_x8 \
  --name move_genotypes \
  --priority normal \
  -y
done

echo -n '' > merge_list.txt
for CHR in {1..22}; do
  echo /mnt/project/QuickdrawsGWAS/copy_ukb22418_c${CHR}_b0_v2 >> merge_list.txt
done
dx upload merge_list.txt --path /QuickdrawsGWAS/

dx run app-swiss-army-knife \
  --folder="QuickdrawsGWAS/" \
  -iin="QuickdrawsGWAS/merge_list.txt" \
  -icmd="plink --bfile /mnt/project/QuickdrawsGWAS/copy_ukb22418_c1_b0_v2 --merge-list merge_list.txt --make-bed --out genotype_500k" \
  --instance-type mem1_ssd1_v2_x8 \
  --name merge_chromosomes \
  --priority normal \
  -y


  for CHR in {1..22}; do
  dx rm QuickdrawsGWAS/copy_ukb22418_c${CHR}_b0_v2.bed 
  dx rm QuickdrawsGWAS/copy_ukb22418_c${CHR}_b0_v2.bim 
  dx rm QuickdrawsGWAS/copy_ukb22418_c${CHR}_b0_v2.fam 
done