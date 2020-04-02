name="simulation1"
region=1000000 # simulate 1 Mbase
nref=100 # number of sample in the reference panel
ntarget=10 # number of imputation target
sample_size=500 # This is the total number msprime simulates, which should be more than nref + ntarget
seed=123

mkdir -p $name
cd $name

# simulate
python3 ../simulate_tree_standalone.py -demo ../CEU.demo -sample_size=$sample_size -out sim -n $region -seed $seed

if [ -f sim.vcf.gz.pos ]; then
  rm sim.vcf.gz.pos
fi
bcftools +fill-AN-AC sim.vcf.gz | bcftools query -f '%POS\t%AN\t%AC\n' > sim.vcf.gz.pos

# subset samples
bcftools query -l sim.vcf.gz > sim.vcf.gz.sample
head -n $nref sim.vcf.gz.sample > ref.sample
tail -n $ntarget sim.vcf.gz.sample > target.sample
bcftools view -S ref.sample sim.vcf.gz -Oz -o ref.vcf.gz
tabix -f ref.vcf.gz
bcftools view -S target.sample sim.vcf.gz -Oz -o target.vcf.gz
tabix -f target.vcf.gz

# create chip data
Rscript ../chip_variants.r sim.vcf.gz 1
bcftools view -R sim.vcf.gz.chip.variants target.vcf.gz -Oz -o target_chip.vcf.gz
tabix target_chip.vcf.gz

java -jar ../beagle.25Nov19.28d.jar ref=ref.vcf.gz gt=target_chip.vcf.gz ap=true out=imputed
tabix -f imputed.vcf.gz

bcftools +fill-tags imputed.vcf.gz | bcftools query -f '%POS\t%AF\n' > imputed.vcf.gz.pos &
bcftools +fill-tags target.vcf.gz | bcftools query -f '%POS\t%AF\n' > target.vcf.gz.pos &
