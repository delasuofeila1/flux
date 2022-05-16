theta=1000; phi=120.0;b=30.0;c=0;d=0.25
while (( $theta <=3000 ))
do
	phi=120.0
	c=0
	b=30.0
	while (( $c<=360 ))
	do      
		
		fluxgo -p fluxzbl testall1
		sed -i -e '/too/d' testall1.txt
		sed -i -e '/!/d' testall1.txt
		python shabi.py 
		((c=$c+1))
		a=`echo $phi+$d | bc`
		sed -i 's/'"$theta"'\,'"$phi"'.*$/'"$theta"'\,'"$a"'/g' testall1.inp
		sed -i 's/'"$theta"'_'"$phi".txt\',heatmap\)'.*$/'"$theta"'_'"$a".txt\',heatmap\)'/g' shabi.py
		sed -i 's/'"$theta"'_'"$phi".png\'\)'.*$/'"$theta"'_'"$a".png\'\)'/g' shabi.py
		a=`echo $b+$d | bc`
		sed -i 's/'0.03' '"-$b"'.*$/'0.03'\ '"-$a"'/g' testall1.inp
	        phi=`echo $phi+$d | bc`
		b=`echo $b+$d | bc`

	done
	((a=$theta+50))
	sed -i 's/'0.03' '"$b"'.*$/'0.03'\ '"210.0"'/g' testall1.inp
	sed -i 's/'"$theta"'\,'"$phi"'.*$/'"$a"'\,'"120.0"'/g' testall1.inp
	sed -i 's/'"$theta"'_'"$phi".png\'\)'.*$/'"$a"'_'"120.0".png\'\)'/g' shabi.py
	sed -i 's/'"$theta"'_'"$phi".txt\',heatmap\)'.*$/'"$a"'_'"120.0".txt\',heatmap\)'/g' shabi.py
	((theta=$theta+50))
done

sed -i 's/'"$theta"'\,'"$phi"'.*$/'"1500"'\,0/g' testall1.inp
sed -i 's/'"$theta"'_'"$phi".png\'\)'.*$/'"1500"'_'"0".png\'\)'/g' shabi.py
