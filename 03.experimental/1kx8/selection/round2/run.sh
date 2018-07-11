for i in {1..7}
do
    name=mota_1kx8_d${i}
    echo $name
    cd $name/strfrags
    bash run.sh ../${name}.pdb
    cd -
done
