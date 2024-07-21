## Run netmhcstab

cat hlas |while read id
do
cut -f 1 ${id}|uniq > hla_type
cut -f 2 ${id} > pep.seq
netMHCstabpan -a $(cat hla_type) -p pep.seq -s -1 -xls -l 8 -xlsfile out/netmhcstab_${id}.xsl
done

sed -i '1,2d' *txt.xsl