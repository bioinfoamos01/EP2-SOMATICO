```bash
          RANCHO DO AMÓS
BEM VINDO AO PIPE DE VARIANTE SOMATICA !
  --------------------------------------
         \   ^__^ 
          \  (oo)\_______
             (__)\       )\/\\
                 ||----w |
                 ||     ||
```
    

# EP2-SOMATICO
O desafio é mapear, anotar, filtrar e comparar arquivos da amostra tumoral WP312, com foco no uso do Mutect2 e PoN (panel of normals).

# Instalação .FASTQ

Instalar (sratoolskit) e fazer Download do arquivo WP312.

```bash
brew install sratoolkit

pip install parallel-fastq-dump

# instalar vdb-config

wget -c https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz
tar -zxvf sratoolkit.3.0.0-ubuntu64.tar.gz
export PATH=$PATH://workspace/somaticoEP1/sratoolkit.3.0.0-ubuntu64/bin/
echo "Aexyo" | sratoolkit.3.0.0-ubuntu64/bin/vdb-config

# baixar arquivo WP312

time parallel-fastq-dump --sra-id SRR8856724 \
--threads 4 \
--outdir ./ \
--split-files \
--gzip
```
# Mapeamento - bwa
Mapear a amostra com o chr9 do hg19
```bash
# bwa para mapear FASTQ
brew install bwa

# baixar chr9-hg19 no UCSC
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr9.fa.gz

# BWA index do arquivo chr9.fa.gz
gunzip chr9.fa.gz
bwa index chr9.fa
brew install samtools
samtools faidx chr9.fa

# combinar com pipes: bwa + samtools view e sort

NOME=WP312; Biblioteca=Nextera; Plataforma=illumina;
bwa mem -t 10 -M -R "@RG\tID:$NOME\tSM:$NOME\tLB:$Biblioteca\tPL:$Plataforma" chr9.fa SRR8856724_1.fastq.gz SRR8856724_2.fastq.gz | samtools view -F4 -Sbu -@2 - | samtools sort -m4G -@2 -o WP312_sorted.bam

# remover duplicatas de PCR
samtools rmdup WP312_sorted.bam WP312_sorted_rmdup.bam

# indexar
samtools index WP312_sorted_rmdup.bam 
```
# Arquivos de referência hg19 de normais
Obter os arquivos de Panel of Normal (PoN), Gnomad AF e adicionar a palavra 'chr' neles.

```bash
# download das sequências
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf.idx
wget -c  https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf
wget -c  https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx

# adicionando 'chr' aos arquivos Gnomad AF para evitar conflitos com GATK
grep "\#" af-only-gnomad.raw.sites.vcf > af-only-gnomad.raw.sites.chr.vcf
grep  "^9" af-only-gnomad.raw.sites.vcf |  awk '{print("chr"$0)}' >> af-only-gnomad.raw.sites.chr.vcf

# indexação
bgzip af-only-gnomad.raw.sites.chr.vcf
tabix -p vcf af-only-gnomad.raw.sites.chr.vcf.gz

# adicionando 'chr' aos arquivos PoN para evitar conflitos com GATK
grep "\#" Mutect2-WGS-panel-b37.vcf > Mutect2-WGS-panel-b37.chr.vcf
grep  "^9" Mutect2-WGS-panel-b37.vcf |  awk '{print("chr"$0)}' >> Mutect2-WGS-panel-b37.chr.vcf

#indexação
bgzip Mutect2-WGS-panel-b37.chr.vcf 
tabix -p vcf Mutect2-WGS-panel-b37.chr.vcf.gz
```
# Arquivos BED
Obter os arquivos BEDs e coberturas
```bash
# instalação bedtools
brew install bedtools

# geração BED a partir de BAM
bedtools bamtobed -i WP312_sorted_rmdup.bam > WP312_sorted_rmdup.bed
bedtools merge -i WP312_sorted_rmdup.bed > WP312_sorted_rmdup_merged.bed
bedtools sort -i WP312_sorted_rmdup_merged.bed > WP312_sorted_rmdup_merged_sorted.bed

# BED com média de cobertura
bedtools coverage -a WP312_sorted_rmdup_merged_sorted.bed \
-b WP312_sorted_rmdup.bam -mean \
> WP312_coverageBed.bed

# Filtragem de cobertura >=20x
cat WP312_coverageBed.bed | \
awk -F "\t" '{if($4>=20){print}}' \
> WP312_coverageBed20x.bed
```
# Chamada de variantes - GATK Mutect Call+ PoN 

Instalação GATK4
```bash
# download GATK4
wget -c https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip
# descompactação
unzip gatk-4.2.2.0.zip
# geração arquivo dict
./gatk-4.2.2.0/gatk CreateSequenceDictionary -R chr9.fa -O chr9.dict
# geração interval_list chr9
./gatk-4.2.2.0/gatk ScatterIntervalsByNs -R chr9.fa -O chr9.interval_list -OT ACGT
# conversão bed para interval_list
./gatk-4.2.2.0/gatk BedToIntervalList -I WP312_coverageBed20x.bed \
-O WP312_coverageBed20x.interval_list -SD chr9.dict
```
Avaliar contaminação e chamar variantes
```bash
# Cálculo de contaminação
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I WP312_sorted_rmdup.bam  \
	-V af-only-gnomad.raw.sites.chr.vcf.gz \
	-L WP312_coverageBed20x.interval_list \
	-O WP312.table

./gatk-4.2.2.0/gatk CalculateContamination \
-I WP312.table \
-O WP312.contamination.table

# Chamada de variantes Mutect2
./gatk-4.2.2.0/gatk Mutect2 \
  -R chr9.fa \
  -I WP312_sorted_rmdup.bam \
  --germline-resource af-only-gnomad.raw.sites.chr.vcf.gz  \
  --panel-of-normals Mutect2-WGS-panel-b37.chr.vcf.gz \
  --disable-sequence-dictionary-validation \
  -L WP312_coverageBed20x.interval_list \
  -O WP312.somatic.pon.vcf.gz

# Filtragem MutectCalls
./gatk-4.2.2.0/gatk FilterMutectCalls \
-R chr9.fa \
-V WP312.somatic.pon.vcf.gz \
--contamination-table WP312.contamination.table \
-O WP312.filtered.pon.vcf.gz

#Comparação vcf-compare
brew install vcftools

# Download dos arquivos VCFs da versão hg19 da análise antiga do Projeto LMA Brasil:
# https://drive.google.com/drive/folders/1m2qmd0ca2Nwb7qcK58ER0zC8-1_9uAiE?usp=sharing
# WP312.filtered.vcf.gz
# WP312.filtered.vcf.gz.tbi

vcf-compare WP312.filtered.pon.vcf.gz WP312.filtered.vcf.gz

```
# Comparar apenas as variantes do chr9 (novo vcf vs download LMA Brasil)
```bash
#adicionar chr no arquivo vcf LMA Brasil
zgrep "\#" WP312.filtered.vcf.gz > header.txt
zgrep -v "\#" WP312.filtered.vcf.gz | awk '{print("chr"$0)}' > variants.txt
cat header.txt variants.txt > WP312.filtered.chr.vcf

#zgrep para separar cabeçalho # e variantes do cromossomo 9:
grep "^\#\|chr9"  WP312.filtered.chr.vcf > WP312.filtered.chr9.vcf

#Compactar o arquivo VCF com o comando bgzip:
bgzip WP312.filtered.chr9.vcf

#Rodar o comando tabix para indexar o arquivo VCF compactado:
tabix -p vcf WP312.filtered.chr9.vcf.gz

#Fazer nova comparação com o arquivo vcf contendo apenas das variantes do chr9:
vcf-compare WP312.filtered.pon.vcf.gz WP312.filtered.chr9.vcf.gz 
```
# RESULTADO:
```bash
#This file was generated by vcf-compare.
#The command line was: vcf-compare(v0.1.14-12-gcdb80b8) WP312.filtered.pon.vcf.gz WP312.filtered.chr9.vcf.gz
#
#VN 'Venn-Diagram Numbers'. Use `grep ^VN | cut -f 2-` to extract this part.
#VN The columns are: 
#VN        1  .. number of sites unique to this particular combination of files
#VN        2- .. combination of files and space-separated number, a fraction of sites in the file
#VN      169     WP312.filtered.chr9.vcf.gz (25.8%)      WP312.filtered.pon.vcf.gz (0.2%)
#VN      487     WP312.filtered.chr9.vcf.gz (74.2%)
#VN      77879   WP312.filtered.pon.vcf.gz (99.8%)
#SN Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
#SN      Number of REF matches:  168
#SN      Number of ALT matches:  166
#SN      Number of REF mismatches:       1
#SN      Number of ALT mismatches:       2
#SN      Number of samples in GT comparison:     0^C
```
