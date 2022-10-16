for i in data/fastq/*/*.fastq.gz; do
    dir=$(dirname $i)
    echo $i
    ncftpput -F -R -z -u geoftp -p "rebUzyi1" ftp-private.ncbi.nlm.nih.gov ./uploads/yuanxue@stanford.edu_6F4lOCKG/${dir}/ $i
done

for i in data/h5ad/*.h5ad; do
    dir=$(dirname $i)
    echo $i
    ncftpput -F -R -z -u geoftp -p "rebUzyi1" ftp-private.ncbi.nlm.nih.gov ./uploads/yuanxue@stanford.edu_6F4lOCKG/${dir}/ $i
done
