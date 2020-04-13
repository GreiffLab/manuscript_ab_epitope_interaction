for f in $(find . -type f -size +99M | grep dataset)
	do git lfs track $f
done
