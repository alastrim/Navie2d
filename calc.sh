if [ -d H ]; then
	rm -rf H
	rm -rf V1
	rm -rf V2
	rm -rf realH
	rm -rf realV1
	rm -rf realV2
fi

mkdir H
mkdir V1
mkdir V2
mkdir realH
mkdir realV1
mkdir realV2
make
./navie $1 $2 $3 $4 $5 $6 $7

