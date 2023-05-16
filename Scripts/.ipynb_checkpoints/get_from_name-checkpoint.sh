wget --no-check-certificate http://neuromorpho.org/api/neuron/name/$1 -O $1_temp.json 
cat $1_temp.json | python -m json.tool >$2/$1.json
rm $1_temp.json
