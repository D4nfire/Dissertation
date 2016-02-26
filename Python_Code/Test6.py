d = []
with open('C:\ThirdYear\Dissertation\Other helpful documents\Gene_ontology_annotations.txt') as source:
    for line in source:
        fields = line.split('\t')
        d.append(fields)
print (d[0][0])