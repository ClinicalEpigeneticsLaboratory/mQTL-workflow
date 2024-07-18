echo Installed tools:

echo Python: $(python3.10 --version)
echo R: $(R --version)
echo bgzip: $(bgzip --version)
echo tabix: $(tabix --version)
echo BCFtools: $(bcftools --version)
echo Illumina CLI: $(array-analysis-cli --version)
echo VEP: $(vep -v)

echo R environment:
echo $(Rscript -e 'sesame::sesame_checkVersion()')

echo Python environment:
echo $(python3.10 -m pip freeze)
