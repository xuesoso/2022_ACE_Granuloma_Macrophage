for i in $(find ../data/*); do
    ext=${i##*.}
    if [ -f ${i} ] && [ ! -e ${i}.md5 ] && [ "${ext}" != 'md5' ]; then
        echo $i
        md5=($(md5sum $i))
        echo $md5 >> $i.md5
    fi
done
